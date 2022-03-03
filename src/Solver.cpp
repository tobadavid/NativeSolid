#include "Solver.h"
#include "Structure.h"
#include "MatVec.h"
#include <iostream>
#include <fstream>
#include <cmath>

Solver::Solver(unsigned int nDof, bool bool_linear)
{
    q.Initialize(nDof, 0.0);
    qdot.Initialize(nDof, 0.0);
    qddot.Initialize(nDof, 0.0);
    q_n.Initialize(nDof, 0.0);
    qdot_n.Initialize(nDof, 0.0);
    qddot_n.Initialize(nDof, 0.0);
    Loads.Initialize(nDof, 0.0);
    Loads_n.Initialize(nDof, 0.0);

    ResetSolution();
    Loads.Reset();
    Loads_n.Reset();

    linear = bool_linear;
}

Solver::~Solver()
{
}

void Solver::Iterate(double &t0, double &tf, Structure *structure)
{
}

CVector &Solver::GetDisp()
{
    return q;
}

CVector &Solver::GetVel()
{
    return qdot;
}

CVector &Solver::GetAcc()
{
    return qddot;
}

CVector &Solver::GetDisp_n()
{
    return q_n;
}

CVector &Solver::GetVel_n()
{
    return qdot_n;
}

CVector &Solver::GetAcc_n()
{
    return qddot_n;
}

CVector &Solver::GetLoads()
{
    return Loads;
}

CVector &Solver::GetAccVar()
{
    return a;
}

CVector &Solver::GetAccVar_n()
{
    return a_n;
}

void Solver::ResetSolution()
{
    q.Reset();
    qdot.Reset();
    qddot.Reset();
    q_n.Reset();
    qdot_n.Reset();
    qddot_n.Reset();
}

void Solver::SaveToThePast()
{
    q_n = q;
    qdot_n = qdot;
    qddot_n = qddot;
    Loads_n = Loads;
}

void Solver::SetInitialState(Config *config, Structure *structure)
{
}

// CLASS ALPHAGENSOLVER
AlphaGenSolver::AlphaGenSolver(unsigned int nDof, double val_rho,
                               bool bool_linear) : Solver(nDof, bool_linear)
{
    a.Initialize(nDof, 0.0);
    a_n.Initialize(nDof, 0.0);

    rho = val_rho;
    alpha_m = (2 * rho - 1) / (rho + 1);
    alpha_f = rho / (rho + 1);
    gamma = 0.5 + alpha_f - alpha_m;
    beta = 0.25 * pow((gamma + 0.5), 2);
    std::cout << "Integration with the alpha-generalized algorithm :" << std::endl;
    std::cout << "rho : " << rho << std::endl;
    std::cout << "alpha_m : " << alpha_m << std::endl;
    std::cout << "alpha_f : " << alpha_f << std::endl;
    std::cout << "gamma : " << gamma << std::endl;
    std::cout << "beta : " << beta << std::endl;
}

AlphaGenSolver::~AlphaGenSolver() 
{
}

CVector &AlphaGenSolver::GetAccVar()
{
    return a;
}

CVector &AlphaGenSolver::GetAccVar_n()
{
    return a_n;
}

void AlphaGenSolver::Iterate(double &t0, double &tf, Structure *structure)
{

    double deltaT(tf - t0), epsilon(1e-6);
    int nMaxIter(1000), nIter(0);

    gammaPrime = gamma / (deltaT * beta);
    betaPrime = (1 - alpha_m) / (pow(deltaT, 2) * beta * (1 - alpha_f));

    //--- Prediction phase ---
    qddot.Reset();
    a.Reset();

    a += ScalVecProd(alpha_f / (1 - alpha_m), qddot_n);
    a -= ScalVecProd(alpha_m / (1 - alpha_m), a_n);

    q = q_n;
    q += ScalVecProd(deltaT, qdot_n);
    q += ScalVecProd((0.5 - beta) * deltaT * deltaT, a_n);
    q += ScalVecProd(deltaT * deltaT * beta, a);

    qdot = qdot_n;
    qdot += ScalVecProd((1 - gamma) * deltaT, a_n);
    qdot += ScalVecProd(deltaT * gamma, a);

    //--- Tangent operator and corrector computation ---
    CVector res(qddot.GetSize(), 0.0);
    CVector Deltaq(qddot.GetSize(), 0.0);
    CMatrix St(qddot.GetSize(), qddot.GetSize(), 0.0);
    ComputeResidual(structure, res);
    while (res.norm() >= epsilon && nIter < nMaxIter)
    {
        St.Reset();
        ComputeTangentOperator(structure, St);
        SolveSys(St, res);
        //*res -= ScalVecProd(double(2),res); //=deltaq
        Deltaq.Reset();
        Deltaq += ScalVecProd(-1, res);
        q += Deltaq;
        qdot += ScalVecProd(gammaPrime, Deltaq);
        qddot += ScalVecProd(betaPrime, Deltaq);
        res.Reset();
        ComputeResidual(structure, res);
        nIter++;
    }
    a += ScalVecProd((1 - alpha_f) / (1 - alpha_m), qddot);
}

void AlphaGenSolver::ComputeRHS(Structure *structure, CVector &RHS)
{

    unsigned long size = q.GetSize();
    CMatrix CC(size, size, 0.0);
    CMatrix KK(size, size, 0.0);
    CVector NonLinTerm(size, 0.0);

    double cos_a;

    if (structure->GetnDof() == 1)
    {
        KK.SetElm(1, 1, structure->Get_Kh());
        CC.SetElm(1, 1, structure->Get_Ch());
    }
    else if (structure->GetnDof() == 2)
    {
        if (linear)
        {
            cos_a = 1.0;
        }
        else
        {
            cos_a = cos(q[1]);
            NonLinTerm[0] = (structure->Get_S()) * sin(q[1]) * pow(qdot[1], 2);
        }
        CC.SetElm(1, 1, structure->Get_Ch());
        CC.SetElm(2, 2, structure->Get_Ca());
        KK.SetElm(1, 1, structure->Get_Kh());
        KK.SetElm(2, 2, structure->Get_Ka());
    }

    RHS += Loads;
    RHS -= MatVecProd(CC, qdot);
    RHS -= MatVecProd(KK, q);
    RHS += NonLinTerm;
}

void AlphaGenSolver::ComputeResidual(Structure *structure, CVector &res)
{

    res.Reset();

    unsigned long size = q.GetSize();
    CMatrix MM(size, size, 0.0);
    double cos_a;

    if (structure->GetnDof() == 1)
    {
        MM.SetElm(1, 1, structure->Get_m());
    }
    else if (structure->GetnDof() == 2)
    {
        if (linear)
        {
            cos_a = 1.0;
        }
        else
        {
            cos_a = cos(q[1]);
        }
        MM.SetElm(1, 1, structure->Get_m());
        MM.SetElm(1, 2, (structure->Get_S()) * cos_a);
        MM.SetElm(2, 1, (structure->Get_S()) * cos_a);
        MM.SetElm(2, 2, structure->Get_If());
    }

    CVector RHS(size, 0.0);
    ComputeRHS(structure, RHS);

    res = MatVecProd(MM, qddot) - RHS;
}

void AlphaGenSolver::ComputeTangentOperator(Structure *structure, CMatrix &St)
{

    St.Reset();

    unsigned long size = q.GetSize();
    CMatrix MM(size, size, 0.0);
    CMatrix Ct(size, size, 0.0);
    CMatrix Kt(size, size, 0.0);

    if (structure->GetnDof() == 1)
    {
        MM.SetElm(1, 1, structure->Get_m());
        Kt.SetElm(1, 1, (structure->Get_Kh()));
        Ct.SetElm(1, 1, (structure->Get_Ch()));
    }
    else if (structure->GetnDof() == 2)
    {
        if (linear)
        {
            MM.SetElm(1, 1, structure->Get_m());
            MM.SetElm(1, 2, (structure->Get_S()));
            MM.SetElm(2, 1, (structure->Get_S()));
            MM.SetElm(2, 2, structure->Get_If());
            Ct.SetElm(1, 1, (structure->Get_Ch()));
            Ct.SetElm(2, 2, (structure->Get_Ca()));
            Kt.SetElm(1, 1, (structure->Get_Kh()));
            Kt.SetElm(2, 2, (structure->Get_Ka()));
        }
        else
        {
            MM.SetElm(1, 1, structure->Get_m());
            MM.SetElm(1, 2, (structure->Get_S()) * cos(q[1]));
            MM.SetElm(2, 1, (structure->Get_S()) * cos(q[1]));
            MM.SetElm(2, 2, structure->Get_If());
            Ct.SetElm(1, 1, -(structure->Get_Ch()));
            Ct.SetElm(1, 2, (structure->Get_S()) * sin(q[1]) * 2 * qdot[1]);
            Ct.SetElm(2, 2, -(structure->Get_Ca()));
            Kt.SetElm(1, 1, -(structure->Get_Kh()));
            Kt.SetElm(1, 2, (structure->Get_S()) * cos(q[1]) * pow(qdot[1], 2));
            Kt.SetElm(2, 2, -(structure->Get_Ka()));
        }
    }

    St += ScalMatProd(betaPrime, MM);
    St += ScalMatProd(gammaPrime, Ct);
    St += Kt;
}

void AlphaGenSolver::ResetSolution()
{
    Solver::ResetSolution();
    a.Reset();
    a_n.Reset();
}

void AlphaGenSolver::SaveToThePast()
{

    Solver::SaveToThePast();
    a_n = a;
}

void AlphaGenSolver::SetInitialState(Config *config, Structure *structure)
{

    if (config->GetRestartSol() == "YES")
    {
        std::string InputFileName = config->GetRestartFile();
        std::string text_line;
        std::string token, tempString;
        size_t pos;
        std::string delimiter = "\t";
        std::ifstream InputFile;
        InputFile.open(InputFileName.c_str(), std::ios::in);
        double buffer[(4 * structure->GetnDof()) + 1];
        int kk = 0;
        int jj;
        while (getline(InputFile, text_line))
        {
            tempString = text_line;
            jj = 0;
            if (kk == 1)
            {
                while ((pos = tempString.find(delimiter)) != std::string::npos)
                {
                    token = tempString.substr(0, pos);
                    tempString.erase(0, pos + delimiter.length());
                    buffer[jj] = atof(token.c_str());
                    jj += 1;
                }
                buffer[jj] = atof(tempString.c_str());

                if (structure->GetnDof() == 1)
                {
                    q_n[0] = buffer[1];
                    qdot_n[0] = buffer[2];
                    qddot_n[0] = buffer[3];
                    a_n[0] = buffer[4];
                }
                else if (structure->GetnDof() == 2)
                {
                    q_n[0] = buffer[1];
                    q_n[1] = buffer[2];
                    qdot_n[0] = buffer[3];
                    qdot_n[1] = buffer[4];
                    qddot_n[0] = buffer[5];
                    qddot_n[1] = buffer[6];
                    a_n[0] = buffer[7];
                    a_n[1] = buffer[8];
                }
                q_n.print();
                qdot_n.print();
                qddot_n.print();
                a_n.print();
            }
            else if (kk == 2)
            {
                while ((pos = tempString.find(delimiter)) != std::string::npos)
                {
                    token = tempString.substr(0, pos);
                    tempString.erase(0, pos + delimiter.length());
                    buffer[jj] = atof(token.c_str());
                    jj += 1;
                }
                buffer[jj] = atof(tempString.c_str());

                if (structure->GetnDof() == 1)
                {
                    q[0] = buffer[1];
                    qdot[0] = buffer[2];
                    qddot[0] = buffer[3];
                    a[0] = buffer[4];
                }
                else if (structure->GetnDof() == 2)
                {
                    q[0] = buffer[1];
                    q[1] = buffer[2];
                    qdot[0] = buffer[3];
                    qdot[1] = buffer[4];
                    qddot[0] = buffer[5];
                    qddot[1] = buffer[6];
                    a[0] = buffer[7];
                    a[1] = buffer[8];
                }
                q.print();
                qdot.print();
                qddot.print();
                a.print();
            }
            kk += 1;
        }
        InputFile.close();
    }
    else
    {
        std::cout << "Setting basic initial conditions for alpha-Gen" << std::endl;
        q.Reset();
        q_n.Reset();
        std::cout << "Read initial configuration" << std::endl;
        q[0] = config->GetInitialDisp();
        if (structure->GetnDof() == 2)
            q[1] = config->GetInitialAngle();
        std::cout << "Initial plunging displacement : " << q[0] << std::endl;
        std::cout << "Initial pitching displacement : " << q[1] << std::endl;

        qdot.Reset();
        qddot.Reset();

        unsigned long size = q.GetSize();
        CVector RHS(size, 0.0);
        CMatrix MM(size, size, 0.0);

        double cos_a;

        if (structure->GetnDof() == 1)
        {
            MM.SetElm(1, 1, structure->Get_m());
        }
        else if (structure->GetnDof() == 2)
        {
            if (linear)
            {
                cos_a = 1.0;
            }
            else
            {
                cos_a = cos(q[1]);
            }
            MM.SetElm(1, 1, structure->Get_m());
            MM.SetElm(1, 2, (structure->Get_S()) * cos_a);
            MM.SetElm(2, 1, (structure->Get_S()) * cos_a);
            MM.SetElm(2, 2, structure->Get_If());
        }

        ComputeRHS(structure, RHS);
        SolveSys(MM, RHS);
        qddot = RHS;
        a = qddot;
    }
}

// CLASS RK4 SOLVER
RK4Solver::RK4Solver(unsigned nDof, bool bool_linear) : Solver(nDof, bool_linear)
{
    size = nDof;
    lastTime = 0.0;
    currentTime = 0.0;
}

RK4Solver::~RK4Solver() {}

void RK4Solver::Iterate(double &t0, double &tf, Structure *structure)
{

    double h = tf - t0;
    lastTime = t0;
    currentTime = tf;

    CVector k1(2 * size);
    CVector k2(2 * size);
    CVector k3(2 * size);
    CVector k4(2 * size);

    CVector state0(2 * size);
    CVector statef(2 * size);
    CVector statef_dot(2 * size);
    state0 = SetState_n();

    CVector TEMP(2 * size);

    EvaluateStateDerivative(lastTime, state0, k1, structure);
    TEMP = state0 + (k1 * (h / 2.0));
    EvaluateStateDerivative(lastTime + h / 2.0, TEMP, k2, structure);
    TEMP = state0 + (k2 * (h / 2.0));
    EvaluateStateDerivative(lastTime + h / 2.0, TEMP, k3, structure);
    TEMP = state0 + (k3 * h);
    EvaluateStateDerivative(lastTime + h, TEMP, k4, structure);

    statef = state0 + ((k1 + k2 * 2.0 + k3 * 2.0 + k4) * (h / 6.0));

    EvaluateStateDerivative(lastTime + h, statef, statef_dot, structure);

    if (structure->GetnDof() == 1)
    {
        q[0] = statef[0];
        qdot[0] = statef[1];
        qddot[0] = statef_dot[1];
    }
    else if (structure->GetnDof() == 2)
    {
        q[0] = statef[0];
        q[1] = statef[1];
        qdot[0] = statef[2];
        qdot[1] = statef[3];
        qddot[0] = statef_dot[2];
        qddot[1] = statef_dot[3];
    }
}

void RK4Solver::EvaluateStateDerivative(double tCurrent, CVector &state, CVector &stateDerivative, Structure *structure)
{

    CVector stateLoads(size, 0.0);
    interpLoads(tCurrent, stateLoads);

    CMatrix MM(size, size, 0.0);
    CMatrix CC(size, size, 0.0);
    CMatrix KK(size, size, 0.0);
    CVector NonLinTerm(size, 0.0);
    CVector RHS(size, 0.0);

    CVector q_current(size, 0.0);
    CVector qdot_current(size, 0.0);
    CVector qddot_current(size, 0.0);

    if (structure->GetnDof() == 1)
    {
        MM.SetElm(1, 1, structure->Get_m());
        KK.SetElm(1, 1, structure->Get_Kh());
        CC.SetElm(1, 1, structure->Get_Ch());
        q_current[0] = state[0];
        qdot_current[0] = state[1];
        stateDerivative[0] = state[1];
    }
    else if (structure->GetnDof() == 2)
    {
        double cos_a;
        if (linear)
        {
            cos_a = 1.0;
        }
        else
        {
            cos_a = cos(state[1]);
            NonLinTerm[0] = (structure->Get_S()) * sin(state[1]) * pow(state[3], 2);
        }
        MM.SetElm(1, 1, structure->Get_m());
        MM.SetElm(1, 2, (structure->Get_S()) * cos_a);
        MM.SetElm(2, 1, (structure->Get_S()) * cos_a);
        MM.SetElm(2, 2, structure->Get_If());
        CC.SetElm(1, 1, structure->Get_Ch());
        CC.SetElm(2, 2, structure->Get_Ca());
        KK.SetElm(1, 1, structure->Get_Kh());
        KK.SetElm(2, 2, structure->Get_Ka());
        q_current[0] = state[0];
        q_current[1] = state[1];
        qdot_current[0] = state[2];
        qdot_current[1] = state[3];
        stateDerivative[0] = state[2];
        stateDerivative[1] = state[3];
    }

    RHS += stateLoads;
    RHS -= MatVecProd(CC, qdot_current);
    RHS -= MatVecProd(KK, q_current);
    RHS += NonLinTerm;

    SolveSys(MM, RHS);
    qddot_current = RHS;

    if (structure->GetnDof() == 1)
    {
        stateDerivative[1] = qddot_current[0];
    }
    else if (structure->GetnDof() == 2)
    {
        stateDerivative[2] = qddot_current[0];
        stateDerivative[3] = qddot_current[1];
    }
}

void RK4Solver::interpLoads(double &tCurrent, CVector &val_loads)
{

    if (lastTime != currentTime)
    {
        val_loads[0] = (Loads[0] - Loads_n[0]) / (currentTime - lastTime) * (tCurrent - lastTime) + Loads_n[0];
        if (size == 2)
            val_loads[1] = (Loads[1] - Loads_n[1]) / (currentTime - lastTime) * (tCurrent - lastTime) + Loads_n[1];
    }
    else
    {
        val_loads[0] = Loads[0];
        if (size == 2)
            val_loads[1] = Loads[1];
    }
}

void RK4Solver::SetInitialState(Config *config, Structure *structure)
{
    if (config->GetRestartSol() == "YES")
    {
    }
    else
    {
        std::cout << "Setting basic initial conditions for RK4" << std::endl;
        q.Reset();
        q_n.Reset();
        std::cout << "Read initial configuration" << std::endl;
        q[0] = config->GetInitialDisp();
        if (structure->GetnDof() == 2)
            q[1] = config->GetInitialAngle();
        std::cout << "Initial plunging displacement : " << q[0] << std::endl;
        std::cout << "Initial pitching displacement : " << q[1] << std::endl;
        qdot.Reset();

        lastTime = 0.0;
        currentTime = 0.0;

        qddot.Reset();
        CVector state(2 * size);
        CVector state_dot(2 * size);

        state = SetState();
        EvaluateStateDerivative(0.0, state, state_dot, structure);

        if (structure->GetnDof() == 1)
        {
            qddot[0] = state_dot[1];
        }
        else if (structure->GetnDof() == 2)
        {
            qddot[0] = state_dot[2];
            qddot[1] = state_dot[3];
        }
    }
}

CVector RK4Solver::SetState()
{
    CVector state(2 * size);

    if (size == 1)
    {
        state[0] = q[0];
        state[1] = qdot[0];
    }
    else if (size == 2)
    {
        state[0] = q[0];
        state[1] = q[1];
        state[2] = qdot[0];
        state[3] = qdot[1];
    }
    else
    {
    }

    return state;
}

CVector RK4Solver::SetState_n()
{
    CVector state(2 * size);

    if (size == 1)
    {
        state[0] = q_n[0];
        state[1] = qdot_n[0];
    }
    else if (size == 2)
    {
        state[0] = q_n[0];
        state[1] = q_n[1];
        state[2] = qdot_n[0];
        state[3] = qdot_n[1];
    }
    else
    {
    }

    return state;
}

// CLASS STATIC SOLVER
StaticSolver::StaticSolver(unsigned nDof, bool bool_linear) : Solver(nDof, bool_linear)
{
    _nDof = nDof;
    KK.Initialize(_nDof, _nDof, 0.0);
}

StaticSolver::~StaticSolver()
{
    std::cout << "NativeSolid::~StaticSolver()" << std::endl;
}

void StaticSolver::SetInitialState(Config *config, Structure *structure)
{
    // Reset displacement, velocity and acceleration
    if (config->GetRestartSol() == "YES")
    {
    }
    else
    {
        std::cout << "Setting basic initial conditions for Static" << std::endl;
        q.Reset();
        q_n.Reset();
        std::cout << "Read initial configuration" << std::endl;
        q[0] = config->GetInitialDisp();
        if (_nDof == 2)
            q[1] = config->GetInitialAngle();
        std::cout << "Initial plunging displacement : " << q[0] << std::endl;
        std::cout << "Initial pitching displacement : " << q[1] << std::endl;
        qdot.Reset();
        qddot.Reset();
    }
    // Fill stiffness matrix
    if (_nDof == 1)
        KK.SetElm(1, 1, structure->Get_Kh());
    else if (_nDof == 2)
    {
        KK.SetElm(1, 1, structure->Get_Kh());
        KK.SetElm(2, 2, structure->Get_Ka());
    }
    else
    {
        std::cerr << "Error in NativeSolid::StaticSolver: Number of degrees of freedom is out of range. nDof = " << _nDof << std::endl;
        throw(-1);
    }
}

void StaticSolver::Iterate(double &t0, double &tf, Structure *structure)
{
    // Solve KK*q = Loads
    q_n = q; // save previous state (used to compute rotation in NativeSolidSolver::computeInterfacePosVel)
    CVector RHS(_nDof, 0.);
    RHS += Loads;
    SolveSys(KK, RHS);
    q = RHS;
}
