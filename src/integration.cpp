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

  if(config->GetIntegrationAlgo() == "ALPHAGEN"){
    solver = new AlphaGenSolver(structure->GetnDof(), config->GetRho());
  }
  else if(config->GetIntegrationAlgo() == "RK4"){
    solver = new RK4Solver(structure->GetnDof());
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

/*
void Integration::SetIntegrationParam(Config* config){

  totTime = config->GetStopTime();
  deltaT = config->GetDeltaT();
  ExtIter = 0;
  algo = config->GetIntegrationAlgo();

  if(algo == "NEWMARK"){
    gamma = 0.5;
    beta = 0.25*pow((gamma+0.5),2);
    alpha_m = 0;
    alpha_f = 0;
    cout << "Integration with the Newmark algorithm :" << endl;
    cout << "gamma : " << gamma << endl;
    cout << "beta : " << beta << endl;
  }
  else if(algo == "HHT"){
    alpha_f = config->GetAlpha();
    alpha_m = 0;
    gamma = 0.5+alpha_f;
    beta = 0.25*pow((gamma+0.5),2);
    cout << "Integration with the Hilbert-Hugues-Taylor algorithm :" << endl;
    cout << "alpha : " << alpha_f << endl;
    cout << "gamma : " << gamma << endl;
    cout << "beta : " << beta << endl;
  }
  else if(algo == "ALPHAGEN"){
    rho = config->GetRho();
    alpha_m = (2*rho-1)/(rho+1);
    alpha_f = rho/(rho+1);
    gamma = 0.5+alpha_f-alpha_m;
    beta = 0.25*pow((gamma+0.5),2);
    cout << "Integration with the alpha-generalized algorithm :" << endl;
    cout << "rho : " << rho << endl;
    cout << "alpha_m : " << alpha_m << endl;
    cout << "alpha_f : " << alpha_f << endl;
    cout << "gamma : " << gamma << endl;
    cout << "beta : " << beta << endl;
  }
  else{
    cout << "Specified integration algorithm is not recognized !" << endl;
    //Faire sortir !!
  }

  gammaPrime = gamma/(deltaT*beta);
  betaPrime = (1-alpha_m)/(pow(deltaT,2)*beta*(1-alpha_f));
  cout << "gammaPrime : " << gammaPrime << endl;
  cout << "betaPrime : " << betaPrime << endl;
}
*/

/*
void Integration::SetLoadsAtTime(Config* config, Structure* structure, const double & time, double FSI_Load){
  if(config->GetForceInputType() == "FILE"){
    string textLine;
    string ForceFileName = config->GetForceInputFile();
    ifstream InputFile;
    InputFile.open(ForceFileName.c_str(), ios::in);
    int kk(0);
    while (getline(InputFile,textLine)){
      if(kk>=structure->GetnDof()) break;
      Loads[kk] = atof(textLine.c_str());
      kk++;
    }
    InputFile.close();
  }
  else if (config->GetForceInputType() == "ANALYTICAL"){
    if(config->GetUnsteady() == "YES" && config->GetAnalyticalFunction() == "SINE"){
      (*Loads)[0] = config->GetAmplitude()*sin(2*PI*config->GetFrequency()*time);
      if(structure->GetnDof() == 2) (*Loads)[1] = config->GetAmplitude()/2.0*sin(2*PI*config->GetFrequency()*time-0.5*PI);
      //if(structure->GetnDof() == 2) (*Loads)[1] = 0;
      //cout << (*Loads)[0] << endl;
    }
    else if (config->GetUnsteady() != "YES" && config->GetAnalyticalFunction() == "CONSTANT"){
    }
  }
  else if (config->GetForceInputType() == "FSI"){
    (*Loads)[0] = FSI_Load;
    if(structure->GetnDof() == 2) (*Loads)[1] = 0.0;
  }
  else cout << "Option for FORCE_INPUT_TYPE has to be specified (FILE or ANALYTICAL)" << endl;
}
*/

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
