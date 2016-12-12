#pragma once

#include "MatVec.h"
#include "config.h"
#include "structure.h"

#define PI 3.14159265


class Integration{

protected:
    double totTime;
    double deltaT;
    unsigned long ExtIter;
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
    std::string algo;

public:
    Integration(Structure *structure);
    ~Integration();
    double GettotTime();
    double GetdeltaT();
    CVector* GetDisp() const;
    void SetDisp(double & val_disp);
    void SetExtIter(unsigned long val_ExtIter);
    unsigned long GetExtIter();
    //void SetLoads()
    CVector* GetVel() const;
    CVector* GetAcc() const;
    CVector* GetAccVar() const;
    CVector* GetDisp_n() const;
    CVector* GetVel_n() const;
    CVector* GetAcc_n() const;
    CVector* GetAccVar_n() const;
    CVector* GetLoads() const;
    void SetIntegrationParam(Config* config);
    void SetLoadsAtTime(Config* config, Structure* structure, const double & time, double FSI_Load);
    void SetStaticLoads(Config* config, Structure* structure);
    void SetInitialConditions(Config* config, Structure* structure);
    void ComputeResidual(CMatrix* M, CMatrix* C, CMatrix* K, CVector* res);
    void ComputeTangentOperator(Structure* structure, CMatrix* St);
    void TemporalIteration(Config *config, Structure *structure);
    void StaticIteration(Config *config, Structure *structure);
    void UpdateSolution();
};
