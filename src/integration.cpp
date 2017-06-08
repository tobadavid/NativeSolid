#include "../include/integration.h"
#include "../include/solver.h"
#include <iostream>
#include <cstdlib>
#include <stdlib.h>
#include <fstream>
#include <cmath>

using namespace std;

Integration::Integration(Config *config, Structure *structure){

  solver = NULL;
  bool linear;

  linear = (config->GetLinearize()) == "YES";

  if(config->GetIntegrationAlgo() == "ALPHAGEN"){
    solver = new AlphaGenSolver(structure->GetnDof(), config->GetRho(), linear);
  }
  else if(config->GetIntegrationAlgo() == "RK4"){
    solver = new RK4Solver(structure->GetnDof(), linear);
  }
  else{

  }
}

Integration::~Integration(){

  if (solver != NULL) delete solver;
}

Solver* Integration::GetSolver(){

  return solver;
}

double Integration::GettotTime(){
  return totTime;
}

double Integration::GetdeltaT(){
  return deltaT;
}

void Integration::SetExtIter(unsigned long val_ExtIter){

    ExtIter = val_ExtIter;
}

unsigned long Integration::GetExtIter(){

    return ExtIter;
}

void Integration::SetInitialConditions(Config* config, Structure* structure){

  ExtIter = 0;
  solver->SetInitialState(config, structure);
  structure->SetCenterOfRotation_Y((solver->GetDisp())[0]);

}

void Integration::TemporalIteration(double& t0, double& tf, Structure *structure){

  int rank = 0;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  double currentTime(tf);
  deltaT = tf - t0;
  solver->Iterate(t0, tf, structure);

  if(rank == 0){
  if(structure->GetnDof() == 1){
    cout << " Time" << "\t" << "Displacement" << "\t" << "Velocity" << "\t" << "Acceleration" << "\t" << endl;
    cout << " " << currentTime << "\t" << (solver->GetDisp())[0] << "\t" << (solver->GetVel())[0] << "\t" << (solver->GetAcc())[0] << endl;
  }
  else if(structure->GetnDof() == 2){
    cout << " Time" << "\t" << "Displacement 1" << "\t" << "Displacement 2" << "\t" << "Velocity 1"  << "\t" << "Velocity 2" << "\t" << "Acceleration 1" << "\t" << "Acceleration 2" << endl;
    cout << currentTime << "\t" << (solver->GetDisp())[0] << "\t" << (solver->GetDisp())[1] << "\t" << (solver->GetVel())[0] << "\t" << (solver->GetVel())[1] << "\t" << (solver->GetAcc())[0] << "\t" << (solver->GetAcc())[1] << endl;
  }
  }
}

/*
void Integration::StaticIteration(Config *config, Structure *structure){
  q->Reset();
  *q += *Loads;
  SolveSys(structure->GetK(),q); //q will be updated with the new solution...
}
*/

void Integration::UpdateSolution(){

  solver->SaveToThePast();

}
