#include "life2d.decl.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
using namespace std;

/*readonly*/ CProxy_Main mainProxy;
/*readonly*/ int arrayDimX;
/*readonly*/ int arrayDimY;
/*readonly*/ int blockDimX;
/*readonly*/ int blockDimY;

// specify the number of worker chares in each dimension
/*readonly*/ int num_chare_x;
/*readonly*/ int num_chare_y;

/*readonly*/ int maxiterations;

char* file_name;

#define TOP			1
#define BOTTOM			2
#define LEFT			3
#define RIGHT			4

#define TOP_LEFT 5
#define TOP_RIGHT 6
#define BOTTOM_LEFT 7
#define BOTTOM_RIGHT 8

#define BROAD_CAST 9
#define GET_RESULT 10

#define ALIVE			1
#define DEAD			0

/** \class Main
 *
 */
class Main : public CBase_Main {

  double startTime;
  double endTime;

public:
  CProxy_GLife array;
  double max_error;
  Main(CkArgMsg* m) {
    if ( m->argc != 7) {
      CkPrintf("%s arg error\n", m->argv[0]);
      CkAbort("Abort");
    }

    // store the main proxy
    mainProxy = thisProxy;

    file_name = m->argv[1];
    maxiterations = atoi(m->argv[2]);
    arrayDimX = atoi(m->argv[3]);
    arrayDimY = atoi(m->argv[4]);
    num_chare_x = atoi(m->argv[5]);
    num_chare_y = atoi(m->argv[6]);


    blockDimX = ceil((double)arrayDimX / num_chare_x);
    blockDimY = ceil((double)arrayDimY / num_chare_y);
    arrayDimX = blockDimX * num_chare_x;
    arrayDimY = blockDimY * num_chare_y;

    array = CProxy_GLife::ckNew(num_chare_x, num_chare_y);
    // start computation
    array.run();
  }

  void done(int totalIter) {
    CkPrintf("Total time %.3f seconds. \n", CkWallTimer()-startTime);
    CkExit();
  }
};

/** \class GLife
 *
 */

class GLife: public CBase_GLife {
  GLife_SDAG_CODE

public:
    double **temperature;
    double **new_temperature;
    double **in_grid;
    int imsg;
    int iterations;
    int neighbors;
    int istart,ifinish,jstart,jfinish;
    double max_error;
    bool topBound, bottomBound, leftBound, rightBound;
    double ellapse_time;
    double start_time;
    // Constructor, initialize values
    GLife() {
      int i, j;
      temperature = new double*[blockDimX+2];
      new_temperature = new double*[blockDimX+2];

      for (i=0; i<blockDimX+2; i++) {
        temperature[i] = new double[blockDimY+2];
        new_temperature[i] = new double[blockDimY+2];
      }

      for(i=0; i<blockDimX+2; ++i) {
        for(j=0; j<blockDimY+2; ++j) {
          temperature[i][j] = 0.;
        } 
      }

      ellapse_time = 6.6;
      imsg = 0;
      iterations = 0;
      neighbors = 0;
      max_error = 0.;
      // determine border conditions
      topBound = bottomBound = leftBound = rightBound = false;
      istart = jstart = 1;
      ifinish = blockDimX+1;
      jfinish = blockDimY+1;

      if(thisIndex.x==0)
      {
        topBound = true;
        istart++;
      }
//      else
//        neighbors++;

      if(thisIndex.x==num_chare_x-1)
      {
        bottomBound = true;
        ifinish--;
      }
//      else
//        neighbors++;

      if(thisIndex.y==0)
      {
        leftBound = true;
        jstart++;
      }
//      else
//        neighbors++;

      if(thisIndex.y==num_chare_y-1)
      {
        rightBound = true;
        jfinish--;
      }
//      else
//        neighbors++;

      int num_bounds = topBound + bottomBound + leftBound + rightBound;
      if (num_bounds == 0)
        neighbors = 8;
      if (num_bounds == 1)
        neighbors = 5;
      if (num_bounds == 2)
        neighbors = 3;
      if (num_bounds == 3)
        neighbors = 1;
      if (num_bounds == 4)
        neighbors = 0;
//      CkPrintf("\nnum_bounds: %d, neighbors: %d\n", num_bounds, neighbors);

      // read data
      if (topBound && leftBound){
        in_grid = new double*[blockDimX * num_chare_x];
        for (i=0; i<blockDimX * num_chare_x; i++) {
          in_grid[i] = new double[blockDimY * num_chare_y];
        }

        for(i=0; i<blockDimX * num_chare_x; ++i) {
          for(j=0; j<blockDimY * num_chare_y; ++j) {
            in_grid[i][j] = 0.;
          }
        }
//        CkPrintf("read data\n");
        ifstream infile( file_name );
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
          in_grid[stoi(record[0])][stoi(record[1])] = 1;
        }
//        printGrid();
      }

//      MyconstrainBC();
//      CkPrintf("bound info: left: %d, right: %d, top: %d, bottom: %d\n",
//              leftBound, rightBound, topBound, bottomBound);
//      dumpMatrix(temperature);
    } // constructor

    void printGrid (void) {
      string out_name = string(file_name) + "." + to_string(maxiterations) + ".csv";
      FILE *fout = fopen(out_name.c_str(), "w"); // printing the result to a file
//      CkPrintf("\n");
      for(int i=0; i<blockDimX * num_chare_x; ++i) {
        for(int j=0; j<blockDimY * num_chare_y; ++j) {
//          CkPrintf("%d ", (int)in_grid[i][j]);
          if (in_grid[i][j] == 1){
            fprintf(fout, "%d,%d\n", i, j);
          }
        }
//        CkPrintf("\n");
      }
      fflush(fout);
      fclose(fout);
    }

    void updateOwnGrid (void) {
      for(int i=0; i<blockDimX; ++i) {
        for(int j=0; j<blockDimY; ++j) {
          in_grid[i][j] = temperature[i+1][j+1];
        }
      }
    }

    void pup(PUP::er &p)
    {
      int i,j;
      p|imsg;
      p|iterations;
      p|neighbors;
      p|istart; p|ifinish; p|jstart; p|jfinish;
      p|leftBound; p|rightBound; p|topBound; p|bottomBound;
      p|ellapse_time;
      p|max_error;

      if (p.isUnpacking()) {
        temperature = new double*[blockDimX+2];
        new_temperature = new double*[blockDimX+2];
        for (i=0; i<blockDimX+2; i++) {
          temperature[i] = new double[blockDimY + 2];
          new_temperature[i] = new double[blockDimY + 2];
        }
      }
      for(i=0;i<blockDimX+2; i++) {
        for(j=0;j<blockDimY+2; j++) {
          p|temperature[i][j];
          p|new_temperature[i][j];
        }
      }
    }

    GLife(CkMigrateMessage* m) { }

    ~GLife() {
      for (int i=0; i<blockDimX+2; i++) {
        delete [] temperature[i];
        delete [] new_temperature[i];
      }
      delete [] temperature; 
      delete [] new_temperature; 
    }

    void send_arrs(void) {
      if (topBound && leftBound) {
//        CkPrintf("send arrays\n");
        for (int x_idx = 0; x_idx < num_chare_x; ++x_idx) {
          for (int y_idx = 0; y_idx < num_chare_y; ++y_idx) {
            int nRows = num_chare_x * blockDimX;
            int nCols = num_chare_y * blockDimY;
            double *out_arr =  new double[blockDimX*blockDimY];
            int out_arr_idx = 0;
            for(int i=x_idx*blockDimX; i<(x_idx+1)*blockDimX; ++i) {
              for(int j=y_idx*blockDimY; j<(y_idx+1)*blockDimY; ++j){
                out_arr[out_arr_idx] = in_grid[i][j];
                out_arr_idx++;
              }
            }
            thisProxy(x_idx, y_idx).receiveGhosts(iterations, BROAD_CAST,
                    blockDimX*blockDimY, out_arr);
          }
        }
      }
    }

    // Debug
    void receive_arrs(void) {
      if (!(topBound && leftBound)) {
        CkPrintf("receive arrays\n");
        dumpMatrix(temperature);
      }
    }

    // Send ghost faces to the six neighbors
    void begin_iteration(void) {
      iterations++;

      if(!topBound)
      {
        double *topGhost =  new double[blockDimY];
        for(int j=0; j<blockDimY; ++j) 
          topGhost[j] = temperature[1][j+1];
        thisProxy(thisIndex.x-1, thisIndex.y).receiveGhosts(iterations, BOTTOM,
                blockDimY, topGhost);
        delete [] topGhost;
      }
      if(!bottomBound)
      {
        double *bottomGhost =  new double[blockDimY];
        for(int j=0; j<blockDimY; ++j) 
          bottomGhost[j] = temperature[blockDimX][j+1];
        thisProxy(thisIndex.x+1, thisIndex.y).receiveGhosts(iterations, TOP,
                blockDimY, bottomGhost);
        delete [] bottomGhost;
      }
      if(!leftBound)
      {
        double *leftGhost =  new double[blockDimX];
        for(int i=0; i<blockDimX; ++i) 
          leftGhost[i] = temperature[i+1][1];
        thisProxy(thisIndex.x, thisIndex.y-1).receiveGhosts(iterations, RIGHT,
                blockDimX, leftGhost);
        delete [] leftGhost;
      }
      if(!rightBound)
      {
        double *rightGhost =  new double[blockDimX];
        for(int i=0; i<blockDimX; ++i) 
          rightGhost[i] = temperature[i+1][blockDimY];
        thisProxy(thisIndex.x, thisIndex.y+1).receiveGhosts(iterations, LEFT,
                blockDimX, rightGhost);
        delete [] rightGhost;
      }
      if( (!topBound) && (!leftBound) ){
        double top_leftGhost_value = temperature[1][1];
        double *top_leftGhost = &top_leftGhost_value;
        thisProxy(thisIndex.x-1, thisIndex.y-1).receiveGhosts(iterations,
                BOTTOM_RIGHT, 1, top_leftGhost);
      }
      if( (!topBound) && (!rightBound) ){
        double top_rightGhost_value = temperature[1][blockDimY];
        double *top_rightGhost = &top_rightGhost_value;
        thisProxy(thisIndex.x-1, thisIndex.y+1).receiveGhosts(iterations,
                BOTTOM_LEFT, 1, top_rightGhost);
      }
      if( (!bottomBound) && (!leftBound) ){
        double bottom_leftGhost_value = temperature[blockDimX][1];
        double *bottom_leftGhost = &bottom_leftGhost_value;
        thisProxy(thisIndex.x+1, thisIndex.y-1).receiveGhosts(iterations,
                TOP_RIGHT, 1, bottom_leftGhost);
      }
      if( (!bottomBound) && (!rightBound) ){
        double bottom_rightGhost_value = temperature[blockDimX][blockDimY];
        double *bottom_rightGhost = &bottom_rightGhost_value;
        thisProxy(thisIndex.x+1, thisIndex.y+1).receiveGhosts(iterations,
                TOP_LEFT, 1, bottom_rightGhost);
      }
    }

    void processGhosts(int dir, int size, double gh[]) {
      switch(dir) {
      case TOP:
        for(int j=0; j<size; ++j) {
          temperature[0][j+1] = gh[j];
        }
        break;
      case BOTTOM:
        for(int j=0; j<size; ++j) {
          temperature[blockDimX+1][j+1] = gh[j];
        }
        break;
      case LEFT:
        for(int i=0; i<size; ++i) {
          temperature[i+1][0] = gh[i];
        }
        break;
      case RIGHT:
        for(int i=0; i<size; ++i) {
          temperature[i+1][blockDimY+1] = gh[i];
        }
        break;
      case TOP_LEFT:
        temperature[0][0] = *gh;
        break;
      case TOP_RIGHT:
        temperature[0][blockDimY+1] = *gh;
        break;
      case BOTTOM_LEFT:
        temperature[blockDimX+1][0] = *gh;
        break;
      case BOTTOM_RIGHT:
        temperature[blockDimX+1][blockDimY+1] = *gh;
        break;
      case BROAD_CAST:
        for(int i=0; i<size; ++i) {
          int row_idx = (int) i / blockDimY + 1;
          int col_idx = (int) i % blockDimY + 1;
          temperature[row_idx][col_idx] = gh[i];
          new_temperature[row_idx][col_idx] = gh[i];
        }
//        CkPrintf("get BROAD_CAST value: %f\n", *gh);
        break;
      default:
        CkAbort("ERROR\n");
      }
    }

    void check_and_compute() {
//      CkPrintf("\ncheck_and_compute iter: %d\n", iterations);
//      dumpMatrix(temperature);
      double temperatureIth = 0.;
      double difference = 0.;
      double **tmp;

      max_error = 0.;
      // When all neighbor values have been received, we update our values and proceed
//      for(int i=istart; i<ifinish; ++i) {
//        for(int j=jstart; j<jfinish; ++j) {
//
//          temperatureIth=temperature[i][j];
//          new_temperature[i][j] = temperatureIth + 1;
//
//          temperatureIth=(temperature[i][j]
//            + temperature[i-1][j]
//            +  temperature[i+1][j]
//            +  temperature[i][j-1]
//            +  temperature[i][j+1]) * 0.2;
//
//          // update relative error
//          difference = temperatureIth-temperature[i][j];
//          // fix sign without fabs overhead
//          if(difference<0) difference *= -1.0;
//          max_error=(max_error>difference) ? max_error : difference;
//          new_temperature[i][j] = temperatureIth;
//        }
//      }


      // update the grid
      for (int iRow = 1; iRow <= blockDimX; ++iRow) {
        for (int iCol = 1; iCol <= blockDimY; ++iCol){
          int nAliveNeighbors = 0;

          for (int jRow = iRow - 1; jRow <= iRow + 1; ++jRow){
            for (int jCol = iCol - 1; jCol <= iCol + 1; ++jCol){
              if ( (jRow != iRow || jCol != iCol) &&
                   temperature[jRow][jCol] == ALIVE){
                ++nAliveNeighbors;
              }
            }
          }

          if (nAliveNeighbors < 2) {
            new_temperature[iRow][iCol] = DEAD;
          }

          if (temperature[iRow][iCol] == ALIVE && (nAliveNeighbors == 2 ||
                                                nAliveNeighbors == 3) ) {
            new_temperature[iRow][iCol] = ALIVE;
          }

          if (nAliveNeighbors > 3){
            new_temperature[iRow][iCol] = DEAD;
          }

          if (temperature[iRow][iCol] == DEAD && nAliveNeighbors == 3){
            new_temperature[iRow][iCol] = ALIVE;
          }
        }
      }

      // new_temperature must be synchronize with
      // temperature when initializing. Because only some values of
      // new_temperature are updated in each iteration and these two matrices are
      // then swapped.
      for (int iRow = 1; iRow <= blockDimX; ++iRow){
        for (int iCol = 1; iCol <= blockDimY; ++ iCol){
          temperature[iRow][iCol] = new_temperature[iRow][iCol];
        }
      }
//      tmp = temperature;
//      temperature = new_temperature;
//      new_temperature = tmp;
    }

    void send_results(void){
      double *out_arr =  new double[blockDimX*blockDimY];
      int out_arr_idx = 0;
      for(int i=0; i<blockDimX; ++i) {
        for(int j=0; j<blockDimY; ++j){
          out_arr[out_arr_idx] = temperature[i+1][j+1];
          out_arr_idx++;
        }
      }
      thisProxy(0,0).receiveResults(iterations, thisIndex.x, thisIndex.y,
              blockDimX*blockDimY, out_arr);
//      CkPrintf("send results\n");
    }

    void processResults(int x_idx, int y_idx, int size, double gh1[]){
//      CkPrintf("receive results from [%d,%d]\n", x_idx, y_idx);
      int out_arr_idx = 0;
      for(int i=x_idx*blockDimX; i<(x_idx+1)*blockDimX; ++i) {
        for(int j=y_idx*blockDimY; j<(y_idx+1)*blockDimY; ++j){
          in_grid[i][j] = gh1[out_arr_idx];
          out_arr_idx++;
//          CkPrintf("%f ", gh1[out_arr_idx]);
        }
//        CkPrintf("\n");
      }
    }

    // for debugging
    void dumpMatrix(double **matrix)
    {
      CkPrintf("[%d,%d]\n",thisIndex.x, thisIndex.y);
      for(int i=0; i<blockDimX+2;++i)
      {
        for(int j=0; j<blockDimY+2;++j)
        {
          CkPrintf("%d ",(int)matrix[i][j]);
        }
        CkPrintf("\n");
      }
    }
};


#include "life2d.def.h"
