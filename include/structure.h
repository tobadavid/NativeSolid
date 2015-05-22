#pragma once

#include "MatVec.h"
#include "config.h"


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
