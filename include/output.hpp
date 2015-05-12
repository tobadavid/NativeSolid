#pragma once

#ifdef HAVE_MPI
    #include "mpi.h"
#endif

#include <iostream>
#include <cstdlib>
#include <stdlib.h>

#include "MatVec.hpp"
#include "config.hpp"
#include "structure.hpp"
#include "integration.hpp"

using namespace std;

class Output{

public:
    Output(void);
    ~Output();
    void WriteHistory(Integration* solver, Structure* structure, ofstream* outputfile, const double & time);
    //void WriteRestart();
};
