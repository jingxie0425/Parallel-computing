#define true 1
#define false 0

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <mpi.h>
#include <time.h>

void perror_exit(const char *message)
{
    perror(message);
    exit(EXIT_FAILURE);
}

void printFile(char **outGrid, int width, int height)
{
    FILE *fout = fopen("./life.2d.csv", "w"); // printing the result to a file
    // with
                                                 // 1 or 0 (1 being an alive cell and 0 a dead cell)
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            if (outGrid[i][j] == '1')
                fprintf(fout, "%d,%d\n", i, j);
        }
    }

    fflush(fout);
    fclose(fout);
}

void update(char **local, char **new, int width_local, int height_local)
{
    // Access the cells of the actual grid, leaving out the auxiliary cells
    // around the grid, yet, taking them into account for the calculations
    for (int y = 1; y <= height_local; y++)
    {
        for (int x = 1; x <= width_local; x++)
        {
            int neighbors = 0;

            // Add the value of each cell to neighbor's variable
            // Adds the ASCII value of each cell, '1' = 49, '0' = 48
            neighbors = local[y - 1][x - 1] + local[y - 1][x] +
                        local[y - 1][x + 1] + local[y][x - 1] +
                        local[y][x + 1] + local[y + 1][x - 1] +
                        local[y + 1][x] + local[y + 1][x + 1];

            // Determine if the current cell is going to be alive or not
            // 387 means that it has 3 neighbors ( (3 * 49) + (5 * 48))
            // 386 means that it has 2 neighbors ( (2 * 49) + (6 * 48))
            if (neighbors == 387 || (neighbors == 386 && (local[y][x] == '1')))
                new[y][x] = '1';
            else
                new[y][x] = '0';
        }
    }
}

int empty(char **local, int width_local, int height_local)
{
    // Checks if local is empty or not (a.k.a. all the cells are dead)
    for (int y = 1; y <= height_local; y++)
    {
        for (int x = 1; x <= width_local; x++)
        {
            if (local[y][x])
                return false;
        }
    }

    return true;
}

int empty_all(char **local, int width_local, int height_local, MPI_Comm *new_comm, int comm_sz)
{
    // Calculates if all subgrids are empty
    int local_flag = empty(local, width_local, height_local),
        global_sum;

    MPI_Allreduce(&local_flag, &global_sum, 1, MPI_INT, MPI_SUM, *new_comm);

    // Compare the number of instances that have an empty grid
    // with the total amount of instances
    return (global_sum == comm_sz);
}

int main(int argc, char *argv[])
{
    // read configurations
    if (argc != 5) {
        perror_exit("argc count ");
        exit(1);
    }
    int nRows = 0, nCols = 0, nTime = 0;

    nTime = atoi(argv[2]);
    nRows = atoi(argv[3]);
    nCols = atoi(argv[4]);
    char *fileName;
    fileName = argv[1];

    // Allocate space for the universal array, the main grid
    char **srcGrid = malloc(nRows * sizeof(char *));
    char *a = malloc(nRows * nCols * sizeof(char));
    if (srcGrid == NULL || a == NULL)
        perror_exit("malloc: ");
    for (int i = 0; i < nRows; i++)
        srcGrid[i] = &a[i * nCols];

    int mpirank, mpisize;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpisize);

    MPI_Comm nextComm;
    int ndims, reorder, periods[2], dim_size[2];

    ndims = 2;
    int blockSize = (int)sqrt(mpisize);

    int nRowsLocal, nColsLocal;

    nRowsLocal = nColsLocal = nRows / blockSize;

    dim_size[0] = blockSize;
    dim_size[1] = blockSize;
    periods[0] = 1;
    periods[1] = 1;
    reorder = 1;

    MPI_Cart_create(MPI_COMM_WORLD, ndims, dim_size, periods, reorder, &nextComm);
    int thisRank, coords[2];

    MPI_Comm_rank(nextComm, &thisRank);
    MPI_Cart_coords(nextComm, thisRank, ndims, coords);

    double t_start = MPI_Wtime();
    double secs;
    FILE *filePtr = NULL;

    // Allocate space
    char **currGrid = malloc((nRowsLocal + 2) * sizeof(char *));
    char *b = malloc((nRowsLocal + 2) * (nColsLocal + 2) * sizeof(char));
    if (currGrid == NULL || b == NULL)
        perror_exit("malloc: ");
    for (int i = 0; i < (nRowsLocal + 2); i++)
        currGrid[i] = &b[i * (nColsLocal + 2)];

    MPI_Status status;

    if (thisRank == 0)
    {
        filePtr = fopen(fileName, "r");
        if (filePtr == NULL)
            perror_exit("fopen: ");

        for (int y = 0; y < nCols; y++)
        {
            for (int x = 0; x < nRows;)
            {
                char c = fgetc(filePtr);
                if ((c != EOF) && (c != '\n'))
                {
                    srcGrid[y][x] = c;
                    x++;
                }
            }
        }

        if (mpisize > 1)
        {
            for (int i = 1; i < mpisize; i++)
            {
                MPI_Cart_coords(nextComm, i, ndims, coords);

                for (int i = 0; i < nRowsLocal; i++)
                {
                    for (int j = 0; j < nColsLocal; j++)
                    {
                        currGrid[i + 1][j + 1] = srcGrid[i + coords[0] * nColsLocal][j + coords[1] * nRowsLocal];
                    }
                }

                MPI_Send(&currGrid[0][0], (nRowsLocal + 2) * (nColsLocal + 2), MPI_CHAR, i, 0, nextComm);
            }
        }

        for (int i = 0; i < nRowsLocal; i++)
        {
            for (int j = 0; j < nColsLocal; j++)
            {
                currGrid[i + 1][j + 1] = srcGrid[i][j];
            }
        }
    }
    else
    {
        MPI_Recv(&currGrid[0][0], (nRowsLocal + 2) * (nColsLocal + 2), MPI_CHAR, 0, MPI_ANY_TAG, nextComm, &status);
    }

    if (thisRank == 0)
    {
        fclose(filePtr);
        filePtr = NULL;
    }

    secs = (MPI_Wtime() - t_start) * 1000;

    // Allocate space for the new array which holds the next generation of the local grid
    char **nextGrid = malloc((nRowsLocal + 2) * sizeof(char *));
    char *c = malloc((nRowsLocal + 2) * (nColsLocal + 2) * sizeof(char));
    if (nextGrid == NULL || c == NULL)
        perror_exit("malloc: ");
    for (int i = 0; i < (nRowsLocal + 2); i++)
        nextGrid[i] = &c[i * (nColsLocal + 2)];

    int iTime = 1;

    MPI_Cart_coords(nextComm, thisRank, ndims, coords);

    // Calculating the coordinates of the neighbours, relative to the current process ones
    int north;
    int south;
    int east;
    int west;

    int north_coords[2];
    int south_coords[2];
    int west_coords[2];
    int east_coords[2];

    north_coords[0] = coords[0] + 1;
    south_coords[0] = coords[0] - 1;
    west_coords[0] = coords[0];
    east_coords[0] = coords[0];

    north_coords[1] = coords[1];
    south_coords[1] = coords[1];
    west_coords[1] = coords[1] - 1;
    east_coords[1] = coords[1] + 1;

    int north_west;
    int north_east;
    int south_west;
    int south_east;

    int north_west_coords[2];
    int north_east_coords[2];
    int south_west_coords[2];
    int south_east_coords[2];

    north_west_coords[0] = coords[0] - 1;
    north_east_coords[0] = coords[0] - 1;
    south_west_coords[0] = coords[0] + 1;
    south_east_coords[0] = coords[0] + 1;

    north_west_coords[1] = coords[1] - 1;
    north_east_coords[1] = coords[1] + 1;
    south_west_coords[1] = coords[1] - 1;
    south_east_coords[1] = coords[1] + 1;

    MPI_Cart_rank(nextComm, north_coords, &north);
    MPI_Cart_rank(nextComm, south_coords, &south);
    MPI_Cart_rank(nextComm, west_coords, &west);
    MPI_Cart_rank(nextComm, east_coords, &east);

    MPI_Cart_rank(nextComm, north_west_coords, &north_west);
    MPI_Cart_rank(nextComm, north_east_coords, &north_east);
    MPI_Cart_rank(nextComm, south_west_coords, &south_west);
    MPI_Cart_rank(nextComm, south_east_coords, &south_east);

    MPI_Datatype vertical_type;

    MPI_Type_vector(nColsLocal, 1, nRowsLocal + 2, MPI_CHAR, &vertical_type);
    MPI_Type_commit(&vertical_type);

    MPI_Request requests_odd[16];
    MPI_Request requests_even[16];

    MPI_Recv_init(&currGrid[0][1], nRowsLocal, MPI_CHAR, north, 1, nextComm, &requests_odd[0]);
    MPI_Send_init(&currGrid[1][1], nRowsLocal, MPI_CHAR, north, 2, nextComm, &requests_odd[1]);
    MPI_Recv_init(&currGrid[nColsLocal + 1][1], nRowsLocal, MPI_CHAR, south, 2, nextComm, &requests_odd[2]);
    MPI_Send_init(&currGrid[nColsLocal][1], nRowsLocal, MPI_CHAR, south, 1, nextComm, &requests_odd[3]);

    MPI_Recv_init(&currGrid[1][nRowsLocal + 1], 1, vertical_type, east, 3, nextComm, &requests_odd[4]);
    MPI_Send_init(&currGrid[1][nRowsLocal], 1, vertical_type, east, 4, nextComm, &requests_odd[5]);
    MPI_Recv_init(&currGrid[1][0], 1, vertical_type, west, 4, nextComm, &requests_odd[6]);
    MPI_Send_init(&currGrid[1][1], 1, vertical_type, west, 3, nextComm, &requests_odd[7]);

    MPI_Recv_init(&currGrid[0][0], 1, MPI_CHAR, north_west, 5, nextComm, &requests_odd[8]);
    MPI_Send_init(&currGrid[1][1], 1, MPI_CHAR, north_west, 6, nextComm, &requests_odd[9]);
    MPI_Recv_init(&currGrid[0][nRowsLocal + 1], 1, MPI_CHAR, north_east, 7, nextComm, &requests_odd[10]);
    MPI_Send_init(&currGrid[1][nRowsLocal], 1, MPI_CHAR, north_east, 8, nextComm, &requests_odd[11]);

    MPI_Recv_init(&currGrid[nColsLocal + 1][0], 1, MPI_CHAR, south_west, 8, nextComm, &requests_odd[12]);
    MPI_Send_init(&currGrid[nColsLocal][1], 1, MPI_CHAR, south_west, 7, nextComm, &requests_odd[13]);
    MPI_Recv_init(&currGrid[nColsLocal + 1][nRowsLocal + 1], 1, MPI_CHAR, south_east, 6, nextComm, &requests_odd[14]);
    MPI_Send_init(&currGrid[nColsLocal][nRowsLocal], 1, MPI_CHAR, south_east, 5, nextComm, &requests_odd[15]);

    MPI_Recv_init(&nextGrid[0][1], nRowsLocal, MPI_CHAR, north, 1, nextComm, &requests_even[0]);
    MPI_Send_init(&nextGrid[1][1], nRowsLocal, MPI_CHAR, north, 2, nextComm, &requests_even[1]);
    MPI_Recv_init(&nextGrid[nColsLocal + 1][1], nRowsLocal, MPI_CHAR, south, 2, nextComm, &requests_even[2]);
    MPI_Send_init(&nextGrid[nColsLocal][1], nRowsLocal, MPI_CHAR, south, 1, nextComm, &requests_even[3]);

    MPI_Recv_init(&nextGrid[1][nRowsLocal + 1], 1, vertical_type, east, 3, nextComm, &requests_even[4]);
    MPI_Send_init(&nextGrid[1][nRowsLocal], 1, vertical_type, east, 4, nextComm, &requests_even[5]);
    MPI_Recv_init(&nextGrid[1][0], 1, vertical_type, west, 4, nextComm, &requests_even[6]);
    MPI_Send_init(&nextGrid[1][1], 1, vertical_type, west, 3, nextComm, &requests_even[7]);

    MPI_Recv_init(&nextGrid[0][0], 1, MPI_CHAR, north_west, 5, nextComm, &requests_even[8]);
    MPI_Send_init(&nextGrid[1][1], 1, MPI_CHAR, north_west, 6, nextComm, &requests_even[9]);
    MPI_Recv_init(&nextGrid[0][nRowsLocal + 1], 1, MPI_CHAR, north_east, 7, nextComm, &requests_even[10]);
    MPI_Send_init(&nextGrid[1][nRowsLocal], 1, MPI_CHAR, north_east, 8, nextComm, &requests_even[11]);

    MPI_Recv_init(&nextGrid[nColsLocal + 1][0], 1, MPI_CHAR, south_west, 8, nextComm, &requests_even[12]);
    MPI_Send_init(&nextGrid[nColsLocal][1], 1, MPI_CHAR, south_west, 7, nextComm, &requests_even[13]);
    MPI_Recv_init(&nextGrid[nColsLocal + 1][nRowsLocal + 1], 1, MPI_CHAR, south_east, 6, nextComm, &requests_even[14]);
    MPI_Send_init(&nextGrid[nColsLocal][nRowsLocal], 1, MPI_CHAR, south_east, 5, nextComm, &requests_even[15]);

    t_start = MPI_Wtime();

    while ((!empty_all(currGrid, nRowsLocal, nColsLocal, &nextComm, mpisize))
    && (iTime <= nTime))
    {

        if ((iTime % 2) == 1)
        {
            MPI_Startall(16, requests_odd);
            MPI_Waitall(16, requests_odd, MPI_STATUSES_IGNORE);
        }
        else
        {
            MPI_Startall(16, requests_even);
            MPI_Waitall(16, requests_even, MPI_STATUSES_IGNORE);
        }

        update(currGrid, nextGrid, nRowsLocal, nColsLocal);

        // The pointer swap
        char **tempGrid = currGrid;
        currGrid = nextGrid;
        nextGrid = tempGrid;

        iTime++;

    } // End of while loop

    // Each MPI process sends its rank to reduction, root MPI process collects the result
    double time_sum = 0;
    double time_max = 0;
    double time_min = 0;
    secs = MPI_Wtime() - t_start;
    MPI_Reduce(&secs, &time_sum, 1, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(&secs, &time_max, 1, MPI_DOUBLE, MPI_MAX, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(&secs, &time_min, 1, MPI_DOUBLE, MPI_MIN, 0,
               MPI_COMM_WORLD);

    if (thisRank == 0) { // If I am not the master instance
        double avg_time = time_sum / (double)mpisize;
        printf("TIME: Min: %f s Avg: %f s Max: %f s\n", time_min, avg_time,
                time_max);
    }
    if (thisRank == 0) // If I am the master instance
    {
        // Copy the dat from my array to univ
        for (int i = 0; i < nRowsLocal; i++)
        {
            for (int j = 0; j < nColsLocal; j++)
            {
                srcGrid[i][j] = currGrid[i + 1][j + 1];
            }
        }

        // For every other instance, receive the local array
        if (mpisize > 1)
        {
            for (int i = 1; i < mpisize; i++)
            {
                MPI_Cart_coords(nextComm, i, ndims, coords);
                MPI_Recv(&currGrid[0][0], (nRowsLocal + 2) * (nColsLocal + 2), MPI_CHAR, i, 0, nextComm, &status);

                // And copy the data to univ
                for (int i = 0; i < nRowsLocal; i++)
                {
                    for (int j = 0; j < nColsLocal; j++)
                    {
                        srcGrid[i + coords[0] * nColsLocal][j + coords[1] * nRowsLocal] = currGrid[i + 1][j + 1];
                    }
                }
            }
        }

        printFile(srcGrid, nRows, nCols);

    }
    else // If I am not the master instance
    {
        // Send the local array to the master instance
        MPI_Send(&currGrid[0][0], (nRowsLocal + 2) * (nColsLocal + 2), MPI_CHAR, 0, 0, nextComm);
    }

    // Deallocate space no longoer needed
    free(b);
    free(currGrid);
    b = NULL;
    currGrid = NULL;

    free(c);
    free(nextGrid);
    c = NULL;
    nextGrid = NULL;

    free(a);
    free(srcGrid);
    a = NULL;
    srcGrid = NULL;

    MPI_Type_free(&vertical_type);

    MPI_Finalize();

    fflush(stdout);

    return 0;

}
