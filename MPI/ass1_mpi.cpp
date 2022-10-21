#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

int main (int argc, char* argv[]){
    int mpirank, mpisize;

    // Initialize the MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpisize);

    int mpiroot = 0;
    int nRows, nCols, nTime;

    // read configurations
    if (mpirank == mpiroot) {
        if (argc != 5) {
            cerr << "incorrect number of input, expected 4" << endl;
            exit(1);
        }
        nRows = atoi(argv[3]);
        nCols = atoi(argv[4]);
        nTime = atoi(argv[2]);
    }

    MPI_Bcast(&nRows, 1, MPI_INT, mpiroot, MPI_COMM_WORLD);
    MPI_Bcast(&nCols, 1, MPI_INT, mpiroot, MPI_COMM_WORLD);
    MPI_Bcast(&nTime, 1, MPI_INT, mpiroot, MPI_COMM_WORLD);

    int nRowsLocal = nRows / mpisize;
    if (mpirank == mpisize - 1){
        nRowsLocal += nRows % mpisize;
    }

    int nRowsLocalWithGhost = nRowsLocal + 2;
    int nColsWithGhost = nCols + 2;

    vector<vector<int> > currGrid(nRowsLocalWithGhost, vector<int>
            (nColsWithGhost, 0));
    vector<vector<int> > nextGrid(nRowsLocalWithGhost, vector<int>
            (nColsWithGhost, 0));


    if (mpirank == mpiroot) {
        // read input file
//        char* fileArg;
//        FILE *filePtr = NULL;
//        fileArg = argv[4];
//        if (nRows < 1 || nCols < 1 || nTime < 1) {
//            cerr << "incorrect input range" << endl;
//            exit(1);
//        }
//
//        filePtr = fopen(fileArg, "r");
//        if (filePtr == NULL) {
//            cerr << "fopen" << endl;
//            exit(1);
//        }
//        cout << "read success: " << endl;
//
//        for (int iRow = 0; iRow < nRows; ++iRow){
//            for (int iCol = 0; iCol < nCols;){
//                char c = fgetc(filePtr);
//                if ((c != EOF) && (c != '\n')) {
//                    srcGrid[iRow][iCol] = c - '0';
//                    iCol++;
//                }
//            }
//        }
        vector<vector<int> > srcGrid(nRows, vector<int> (nCols, 0));

        ifstream infile( argv[1] );
        if (!infile.is_open()){
            cerr << "fopen wrong" << endl;
            exit(1);
        }
        while (infile)
        {
            string s;
            if (!getline( infile, s )) break;

            istringstream ss( s );
            vector <string> record;

            while (ss)
            {
                string s;
                if (!getline( ss, s, ',' )) break;
                record.push_back( s );
            }
            srcGrid[stoi(record[0])][stoi(record[1])] = 1;
        }
//        if (!infile.eof())
//        {
//            cerr << "Fooey!\n";
//        }



//        cout << "read from src: " << endl;
//        for (int iRow = 0; iRow < nRows; ++iRow){
//            for (int iCol = 0; iCol < nCols; ++iCol){
//                cout << srcGrid[iRow][iCol] << " ";
//            }
//            cout << endl;
//        }
//        cout << endl;

        // formulate its own grid
        for (int iRow = 0; iRow < nRowsLocal; ++iRow){
            for (int iCol = 0; iCol < nCols; ++iCol){
                currGrid[iRow + 1][iCol + 1] = srcGrid[iRow][iCol];
            }
        }

        // send to other ranks
        for (int destRank = 1; destRank < mpisize; ++destRank){
            int nSend = nRows / mpisize;
            int nRowInPrevRank = nSend;
            // last rank
            if (destRank == mpisize - 1){
                nSend += nRows % mpisize;
            }

            vector<int> buff(nCols, 0);
            for (int iRow = destRank * nRowInPrevRank; iRow < (destRank *
            nRowInPrevRank + nSend); ++iRow){
                MPI_Send(&srcGrid[iRow][0], nCols, MPI_INT, destRank, 0,
                         MPI_COMM_WORLD);
            }
        }
    }

    if (mpirank != mpiroot){
        // send row by row
        for (int iRow = 1; iRow <= nRowsLocal; ++iRow){
            MPI_Recv(&currGrid[iRow][1], nCols, MPI_INT, mpiroot, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    // rand values
//    for (int iRow = 1; iRow <= nRowsLocal; ++iRow){
//        for (int iCol = 1; iCol <= nCols; ++iCol){
//            currGrid[iRow][iCol] = rand() % 2;
//        }
//    }

    int upperNeighbor = (mpirank == 0) ? mpisize - 1 : mpirank - 1;
    int lowerNeighbor = (mpirank == mpisize - 1) ? 0 : mpirank + 1;

    const int ALIVE = 1;
    const int DEAD = 0;

    double t_start = MPI_Wtime();
    // time loop
    for (int iTime = 0; iTime < nTime; ++iTime){
        MPI_Send(&currGrid[1][0], nColsWithGhost, MPI_INT, upperNeighbor, 0,
                MPI_COMM_WORLD);
        MPI_Send(&currGrid[nRowsLocal][0], nColsWithGhost, MPI_INT,
                lowerNeighbor, 0, MPI_COMM_WORLD);

        MPI_Recv(&currGrid[nRowsLocal + 1][0], nColsWithGhost, MPI_INT,
                lowerNeighbor, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&currGrid[0][0], nColsWithGhost, MPI_INT, upperNeighbor, 0,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if (mpirank == mpiroot){
            for (int iCol = 0; iCol < nColsWithGhost; ++iCol){
                currGrid[0][iCol] = 0;
            }
        }
        if (mpirank == mpisize - 1){
            for (int iCol = 0; iCol < nColsWithGhost; ++iCol){
                currGrid[nRowsLocalWithGhost-1][iCol] = 0;
            }
        }

        // ghost cols
        for (int iRow = 0; iRow < nRowsLocalWithGhost; ++iRow){
            currGrid[iRow][0] = 0;//currGrid[iRow][nCols];
            currGrid[iRow][nCols + 1] = 0;//currGrid[iRow][1];
        }



        // update the grid
        for (int iRow = 1; iRow <= nRowsLocal; ++iRow) {
            for (int iCol = 1; iCol <= nCols; ++iCol){
                int nAliveNeighbors = 0;

                for (int jRow = iRow - 1; jRow <= iRow + 1; ++jRow){
                    for (int jCol = iCol - 1; jCol <= iCol + 1; ++jCol){
                        if ( (jRow != iRow || jCol != iCol) &&
                        currGrid[jRow][jCol] == ALIVE){
                            ++nAliveNeighbors;
                        }
                    }
                }

                if (nAliveNeighbors < 2) {
                    nextGrid[iRow][iCol] = DEAD;
                }

                if (currGrid[iRow][iCol] == ALIVE && (nAliveNeighbors == 2 ||
                nAliveNeighbors == 3) ) {
                    nextGrid[iRow][iCol] = ALIVE;
                }

                if (nAliveNeighbors > 3){
                    nextGrid[iRow][iCol] = DEAD;
                }

                if (currGrid[iRow][iCol] == DEAD && nAliveNeighbors == 3){
                    nextGrid[iRow][iCol] = ALIVE;
                }
            }
        }

        // swap the grids
//        currGrid.swap(nextGrid);
        for (int iRow = 1; iRow <= nRowsLocal; ++iRow){
            for (int iCol = 1; iCol <= nCols; ++ iCol){
                currGrid[iRow][iCol] = nextGrid[iRow][iCol];
            }
        }

    }
    double secs = MPI_Wtime() - t_start;

    // Each MPI process sends its rank to reduction, root MPI process collects the result
    double time_sum = 0;
    double time_max = 0;
    double time_min = 0;
    MPI_Reduce(&secs, &time_sum, 1, MPI_DOUBLE, MPI_SUM, mpiroot,
               MPI_COMM_WORLD);
    MPI_Reduce(&secs, &time_max, 1, MPI_DOUBLE, MPI_MAX, mpiroot,
               MPI_COMM_WORLD);
    MPI_Reduce(&secs, &time_min, 1, MPI_DOUBLE, MPI_MIN, mpiroot,
               MPI_COMM_WORLD);
    if(mpirank == mpiroot)
    {
        double avg_time = time_sum / (double)mpisize;
        cout << "TIME: Min: " << time_min << " s " << "Avg: " << avg_time <<
        " s " << "Max: " << time_max << endl;
    }

    // display
    if (mpirank != mpiroot) {
        // send row by row
        for (int iRow = 1; iRow <= nRowsLocal; ++iRow) {
            MPI_Send(&currGrid[iRow][1], nCols, MPI_INT, mpiroot, 0,
                     MPI_COMM_WORLD);
        }
    }

    if (mpirank == mpiroot) {
        string out_name = string(argv[1]) + "." + to_string(nTime) + ".csv";
        FILE *fout = fopen(out_name.c_str(), "w"); // printing the result to a file
        // with

        // print its own data
        for (int iRow = 1; iRow <= nRowsLocal; ++iRow) {
            for (int iCol = 1; iCol <= nCols; ++iCol) {
                if (currGrid[iRow][iCol] == 1) {
                    fprintf(fout, "%d,%d\n", iRow-1, iCol-1);
//                    cout << iRow-1 << "," << iCol-1 << endl;
                }
            }
        }

        // receive from other ranks
        for (int sourceRank = 1; sourceRank < mpisize; ++sourceRank) {
            int nRecv = nRows / mpisize;
            int nRowsPerRank = nRecv;
            // last rank
            if (sourceRank == mpisize - 1) {
                nRecv += nRows % mpisize;
            }

            vector<int> buff(nCols, 0);
            for (int iRecv = 0; iRecv < nRecv; ++iRecv) {
                MPI_Recv(&buff[0], nCols, MPI_INT, sourceRank, 0,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                int col_cnt = 0;
                for (int i : buff) {
                    if (i == 1){
                        fprintf(fout, "%d,%d\n", sourceRank * nRowsPerRank +
                        iRecv, col_cnt);
//                        cout << sourceRank * nRowsPerRank + iRecv <<
//                             "," << col_cnt << endl;
                    }
                    col_cnt++;
                }
            }
        }

        fflush(fout);
        fclose(fout);
    }

    MPI_Finalize();
    return 0;
}