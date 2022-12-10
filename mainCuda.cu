#include "main.cu"
#include <stdio.h>
// #include <sys/time.h>
#include <stdlib.h>
#include <time.h>
#include <chrono>
// #include "hrtime.h"
using namespace std;

int step;
double steps, endtime;


int main(int argc, char** argv){
    char oFlag = 0;
    if (argc > 2)
        oFlag = 1;

    SetStartTime(0);
    timestep = 0.25;
    LoadFile(argv[1]);

    SetOutToFile("test.txt", 10);

    initialize();

    endtime = 30000;
    steps = endtime/timestep;

    double elapsed = 0;

    auto t_start = chrono::high_resolution_clock::now();
    // the work...

    for(step = 0; step < steps; step++){
        if (oFlag)
            Output();
        Step();
    }
    cudaDeviceSynchronize();
    auto t_end = chrono::high_resolution_clock::now();
    elapsed = (t_end - t_start).count();

    fprintf(stdout, "\n==============================\nCode took: %f seconds\n", elapsed / 1000000000.);

    fclose(stdout);
}
