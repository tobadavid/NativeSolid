#pragma once

#ifdef HAVE_MPI
    #include "mpi.h"
#endif

#include <iostream>
#include <cstdlib>
#include <stdlib.h>
#include <fstream>
#include <sstream>

#include "MatVec.hpp"
#include "config.hpp"


using namespace std;

class Structure{

protected:
    unsigned int nDof;
    double m;
    double Kh;
    double Ka;
    double Ch;
    double Ca;
    double c;
    double xf;
    double S;
    double Ia;
    CMatrix* M;
    CMatrix* C;
    CMatrix* K;
    //CVector* InitialSolution;


public:
    Structure(Config* config);
    virtual ~Structure();

    void SetStructuralMatrices(Config* config);
    //void SetSolutionInTime();
    CMatrix* GetM();
    CMatrix* GetC();
    CMatrix* GetK();
    unsigned int GetnDof();

};
