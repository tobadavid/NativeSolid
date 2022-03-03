#include "integration.h"
#include "solver.h"
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <stdlib.h>
#include <fstream>
#include <cmath>

Integration::Integration(Config *config, Structure *structure)
{
    solver = NULL;
    bool linear;

    linear = (config->GetLinearize()) == "YES";

    if (config->GetUnsteady() == "YES")
    {
        if (config->GetIntegrationAlgo() == "ALPHAGEN")
        {
            solver = new AlphaGenSolver(structure->GetnDof(), config->GetRho(), linear);
        }
        else if (config->GetIntegrationAlgo() == "RK4")
        {
            solver = new RK4Solver(structure->GetnDof(), linear);
        }
    }
    else
        solver = new StaticSolver(structure->GetnDof(), linear);
}

Integration::~Integration()
{
    if (solver != NULL)
        delete solver;
}

Solver *Integration::GetSolver()
{

    return solver;
}

double Integration::GettotTime()
{
    return totTime;
}

double Integration::GetdeltaT()
{
    return deltaT;
}

void Integration::SetExtIter(unsigned long val_ExtIter)
{

    ExtIter = val_ExtIter;
}

unsigned long Integration::GetExtIter()
{

    return ExtIter;
}

void Integration::SetInitialConditions(Config *config, Structure *structure)
{
    ExtIter = 0;
    solver->SetInitialState(config, structure);
    structure->SetCenterOfRotation_Y((solver->GetDisp())[0]);
}

void Integration::TemporalIteration(double &t0, double &tf, Structure *structure)
{
    int rank = 0;

#ifdef HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    double currentTime(tf);
    deltaT = tf - t0;
    solver->Iterate(t0, tf, structure);

    if (rank == 0)
    {
        if (structure->GetnDof() == 1)
        {
            std::cout << std::fixed
                 << std::setw(10) << "Time"
                 << std::setw(15) << "Displacement"
                 << std::setw(15) << "Velocity"
                 << std::setw(15) << "Acceleration" << std::endl;
            std::cout << std::fixed
                 << std::setw(10) << currentTime
                 << std::setw(15) << (solver->GetDisp())[0]
                 << std::setw(15) << (solver->GetVel())[0]
                 << std::setw(15) << (solver->GetAcc())[0] << std::endl;
        }
        else if (structure->GetnDof() == 2)
        {
            std::cout << std::fixed
                 << std::setw(10) << "Time"
                 << std::setw(15) << "Displacement1"
                 << std::setw(15) << "Displacement2"
                 << std::setw(15) << "Velocity1"
                 << std::setw(15) << "Velocity2"
                 << std::setw(15) << "Acceleration1"
                 << std::setw(15) << "Acceleration2" << std::endl;
            std::cout << std::fixed
                 << std::setw(10) << currentTime
                 << std::setw(15) << (solver->GetDisp())[0]
                 << std::setw(15) << (solver->GetDisp())[1]
                 << std::setw(15) << (solver->GetVel())[0]
                 << std::setw(15) << (solver->GetVel())[1]
                 << std::setw(15) << (solver->GetAcc())[0]
                 << std::setw(15) << (solver->GetAcc())[1] << std::endl;
        }
    }
}

void Integration::StaticIteration(Structure *structure)
{
    int rank = 0;

#ifdef HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    double t0 = 0.;
    double tf = 0.;
    deltaT = 0.;
    solver->Iterate(t0, tf, structure);

    if (rank == 0)
    {
        if (structure->GetnDof() == 1)
        {
            std::cout << std::fixed
                 << std::setw(10) << "Time"
                 << std::setw(15) << "Displacement" << std::endl;
            std::cout << std::fixed
                 << std::setw(10) << tf
                 << std::setw(15) << (solver->GetDisp())[0] << std::endl;
        }
        else if (structure->GetnDof() == 2)
        {
            std::cout << std::fixed
                 << std::setw(10) << "Time"
                 << std::setw(15) << "Displacement1"
                 << std::setw(15) << "Displacement2" << std::endl;
            std::cout << std::fixed
                 << std::setw(10) << tf
                 << std::setw(15) << (solver->GetDisp())[0]
                 << std::setw(15) << (solver->GetDisp())[1] << std::endl;
        }
    }
}

void Integration::UpdateSolution()
{
    solver->SaveToThePast();
}
