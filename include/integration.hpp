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

#define PI 3.14159265

using namespace std;

class Integration{

protected:
    double totTime;
    double deltaT;
    double beta;
    double gamma;
    double alpha_m;
    double alpha_f;
    double rho;
    double gammaPrime;
    double betaPrime;
    CVector* q;
    CVector* qdot;
    CVector* qddot;
    CVector* a;
    CVector* q_n;
    CVector* qdot_n;
    CVector* qddot_n;
    CVector* a_n;
    CVector* Loads;
    string algo;

public:
    Integration(Structure *structure);
    ~Integration();
    double GettotTime();
    double GetdeltaT();
    CVector* GetDisp() const;
    CVector* GetVel() const;
    CVector* GetAcc() const;
    void SetIntegrationParam(Config* config);
    void SetLoadsAtTime(Config* config, Structure* structure, const double & time);
    void SetInitialConditions(Config* config, Structure* structure);
    void ComputeResidual(CMatrix* M, CMatrix* C, CMatrix* K, CVector* res);
    void ComputeTangentOperator(Structure* structure, CMatrix* St);
    void TemporalIteration(Config *config, Structure *structure);
    void UpdateSolution();
};
