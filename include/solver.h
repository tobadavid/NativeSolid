#pragma once

#include "MatVec.h"
#include "structure.h"
#include "config.h"

class Solver{

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

public:
    Solver(unsigned int nDof);
    virtual ~Solver();
    virtual void Iterate(double& t0, double& tf, Structure* structure);
    virtual CVector & GetDisp();
    virtual CVector & GetVel();
    virtual CVector & GetAcc();
    virtual CVector & GetDisp_n();
    virtual CVector & GetVel_n();
    virtual CVector & GetAcc_n();
    virtual CVector & GetLoads();
    virtual CVector & GetAccVar();
    virtual CVector & GetAccVar_n();
    virtual void ResetSolution();
    virtual void SaveToThePast();
    virtual void SetInitialState(Config *config, Structure* structure);

};

class AlphaGenSolver : public Solver {

protected:
    double beta;
    double gamma;
    double alpha_m;
    double alpha_f;
    double rho;
    double gammaPrime;
    double betaPrime;

public:
  AlphaGenSolver(unsigned int nDof, double val_rho);
  ~AlphaGenSolver();
  CVector & GetAccVar();
  CVector & GetAccVar_n();
  virtual void Iterate(double& t0, double& tf, Structure* structure);
  void ComputeResidual(const CMatrix &M, const CMatrix &C, const CMatrix &K, CVector &res);
  void ComputeTangentOperator(Structure* structure, CMatrix & St);
  void ResetSolution();
  void SaveToThePast();
  virtual void SetInitialState(Config *config, Structure* structure);

};


class RK4Solver : public Solver {

protected:
  unsigned int size;
  //CVector* Q;
  //CVector* Q_dot;
  double lastTime;
  double currentTime;

public:
    RK4Solver(unsigned nDof);
    ~RK4Solver();
    virtual void Iterate(double &t0, double &tf, Structure* structure);
    void EvaluateStateDerivative(double tCurrent, CVector& state, CVector& stateDerivative, Structure* structure);
    void interpLoads(double& tCurrent, CVector& val_loads);
    virtual void SetInitialState(Config* config, Structure* structure);
    CVector SetState();
    CVector SetState_n();

};
