#pragma once

#include "MatVec.h"
#include "structure.h"
#include "config.h"

class Solver
{
protected:
    CVector q;
    CVector qdot;
    CVector qddot;
    CVector q_n;
    CVector qdot_n;
    CVector qddot_n;
    CVector Loads;
    CVector Loads_n;
    CVector a;
    CVector a_n;
    bool linear;

public:
    Solver(unsigned int nDof, bool bool_linear);
    virtual ~Solver();
    virtual void Iterate(double &t0, double &tf, Structure *structure);
    virtual CVector &GetDisp();
    virtual CVector &GetVel();
    virtual CVector &GetAcc();
    virtual CVector &GetDisp_n();
    virtual CVector &GetVel_n();
    virtual CVector &GetAcc_n();
    virtual CVector &GetLoads();
    virtual CVector &GetAccVar();
    virtual CVector &GetAccVar_n();
    virtual void ResetSolution();
    virtual void SaveToThePast();
    virtual void SetInitialState(Config *config, Structure *structure);
};

class AlphaGenSolver : public Solver
{
protected:
    double beta;
    double gamma;
    double alpha_m;
    double alpha_f;
    double rho;
    double gammaPrime;
    double betaPrime;

public:
    AlphaGenSolver(unsigned int nDof, double val_rho, bool bool_linear);
    ~AlphaGenSolver();
    CVector &GetAccVar();
    CVector &GetAccVar_n();
    virtual void Iterate(double &t0, double &tf, Structure *structure);
    void ComputeRHS(Structure *structure, CVector &RHS);
    void ComputeResidual(Structure *structure, CVector &res);
    void ComputeTangentOperator(Structure *structure, CMatrix &St);
    void ResetSolution();
    void SaveToThePast();
    virtual void SetInitialState(Config *config, Structure *structure);
};

class RK4Solver : public Solver
{
protected:
    unsigned int size;
    double lastTime;
    double currentTime;

public:
    RK4Solver(unsigned nDof, bool bool_linear);
    ~RK4Solver();
    virtual void Iterate(double &t0, double &tf, Structure *structure);
    void EvaluateStateDerivative(double tCurrent, CVector &state,
                                 CVector &stateDerivative,
                                 Structure *structure);
    void interpLoads(double &tCurrent, CVector &val_loads);
    virtual void SetInitialState(Config *config, Structure *structure);
    CVector SetState();
    CVector SetState_n();
};

class StaticSolver : public Solver
{
protected:
    unsigned int _nDof;
    CMatrix KK;

public:
    StaticSolver(unsigned nDof, bool bool_linear);
    ~StaticSolver();

    virtual void Iterate(double &t0, double &tf, Structure *structure);
    virtual void SetInitialState(Config *config, Structure *structure);
};
