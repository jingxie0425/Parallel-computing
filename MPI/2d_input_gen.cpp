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
    int nRows = atoi(argv[2]); int nCols = atoi(argv[3]);
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

    string out_name = string(argv[1]) + ".txt";
    FILE *fout = fopen(out_name.c_str(), "w"); // printing the result to a file

    for (int iRow = 0; iRow < nRows; ++iRow) {
        for (int iCol = 0; iCol < nCols; ++iCol) {
            fprintf(fout, "%d", srcGrid[iRow][iCol]);
        }
        fprintf(fout, "\n");
    }

    fflush(fout);
    fclose(fout);

    return 0;
}