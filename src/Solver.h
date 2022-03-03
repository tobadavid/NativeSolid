#pragma once

#include "MatVec.h"
#include "Structure.h"
#include "Config.h"

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
    double beta;
    double gamma;
    double alpha_m;
    double alpha_f;
    double rho;
    double gammaPrime;
    double betaPrime;

public:
    AlphaGenSolver(unsigned int nDof, double val_rho, bool bool_linear);
    virtual ~AlphaGenSolver() override;
    virtual CVector &GetAccVar() override;
    virtual CVector &GetAccVar_n() override;
    virtual void Iterate(double &t0, double &tf, Structure *structure) override;
    virtual void ResetSolution() override;
    virtual void SaveToThePast() override;
    virtual void SetInitialState(Config *config, Structure *structure) override;

private:
    void ComputeRHS(Structure *structure, CVector &RHS);
    void ComputeResidual(Structure *structure, CVector &res);
    void ComputeTangentOperator(Structure *structure, CMatrix &St);
};

class RK4Solver : public Solver
{
    unsigned int size;
    double lastTime;
    double currentTime;

public:
    RK4Solver(unsigned nDof, bool bool_linear);
    virtual ~RK4Solver() override;
    virtual void Iterate(double &t0, double &tf, Structure *structure) override;
    virtual void SetInitialState(Config *config, Structure *structure) override;

private:
    void EvaluateStateDerivative(double tCurrent, CVector &state,
                                 CVector &stateDerivative,
                                 Structure *structure);
    void interpLoads(double &tCurrent, CVector &val_loads);
    CVector SetState();
    CVector SetState_n();
};

class StaticSolver : public Solver
{
    unsigned int _nDof;
    CMatrix KK;

public:
    StaticSolver(unsigned nDof, bool bool_linear);
    virtual ~StaticSolver() override;

    virtual void Iterate(double &t0, double &tf, Structure *structure) override;
    virtual void SetInitialState(Config *config, Structure *structure) override;
};
