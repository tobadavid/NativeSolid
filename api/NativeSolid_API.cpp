#include "NativeSolid_API.h"

using namespace std;

const int MASTER_NODE = 0;


NativeSolidSolver::NativeSolidSolver(string str):confFile(str){

}

NativeSolidSolver::~NativeSolidSolver(){}

void NativeSolidSolver::initialize(bool FSIComp){

  int rank = MASTER_NODE;
  int size = 1;
  //double currentTime(0),totTime(0),deltaT(0);

  outputFile.open("NativeHistory.dat", ios::out);

/*--- MPI initialization, and buffer setting ---*/
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  if(rank == MASTER_NODE){
    if(FSIComp)cout << endl << "***************************** Setting NativeSolid for FSI simulation *****************************" << endl;
    else cout << endl <<"***************************** Setting NativeSolid for CSD simulation *****************************" << endl;
  }

    /*--- Initialize the main containers ---*/
    config = new Config(confFile);
    output = new Output();

    /*--- Read CSD configuration file ---*/
    cout << endl << "\n----------------------- Reading Native configuration file ----------------------" << endl;
    config->ReadConfig();

    /*--- Initialize structural container and create the structural model ---*/
    cout << endl << "\n----------------------- Creating the structural model ----------------------" << endl;
    structure = new Structure(config);
    structure->SetStructuralMatrices(config);

    /*--- Initialize the temporal integrator ---*/
    cout << endl << "\n----------------------- Setting integration parameter ----------------------" << endl;
    solver = new Integration(structure);
    if(config->GetUnsteady() == "YES"){
      solver->SetIntegrationParam(config);
      solver->SetInitialConditions(config, structure);
    }

  if(rank == MASTER_NODE){
    if(FSIComp)cout << endl << "***************************** NativeSolid is set for FSI simulation *****************************" << endl;
    else cout << endl <<"***************************** NativeSolid is set for CSD simulation *****************************" << endl;
  }

}

void NativeSolidSolver::exit(){

  int rank = MASTER_NODE;
  int size = 1;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  outputFile.close();
  if(rank == MASTER_NODE){
    cout << "Solid history file is closed." << endl;
    cout << endl << "***************************** Exit NativeSolid *****************************" << endl;
  }

  delete config;
  delete structure;
  delete solver;
  delete output;

}

void NativeSolidSolver::inputFluidLoads(double currentTime, double FSI_Load){
  solver->SetLoadsAtTime(config, structure, currentTime, FSI_Load);
}

void NativeSolidSolver::timeIteration(double currentTime){

  int rank = MASTER_NODE;
  int size = 1;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  solver->TemporalIteration(config, structure);

  if(rank == MASTER_NODE){
  if(structure->GetnDof() == 1){
    cout << " Time" << "\t" << "Displacement" << "\t" << "Velocity" << "\t" << "Acceleration" << "\t" << endl;
    cout << " " << currentTime << "\t" << (*(solver->GetDisp()))[0] << "\t" << (*(solver->GetVel()))[0] << "\t" << (*(solver->GetAcc()))[0] << endl;
  }
  else if(structure->GetnDof() == 2){
    cout << " Time" << "\t" << "Displacement 1" << "\t" << "Displacement 2" << "\t" << "Velocity 1"  << "\t" << "Velocity 2" << "\t" << "Acceleration 1" << "\t" << "Acceleration 2" << endl;
    cout << currentTime << "\t" << (*(solver->GetDisp()))[0] << "\t" << (*(solver->GetDisp()))[1] << "\t" << (*(solver->GetVel()))[0] << "\t" << (*(solver->GetVel()))[1] << "\t" << (*(solver->GetAcc()))[0] << "\t" << (*(solver->GetAcc()))[1] << endl;
  }
  }

}

void NativeSolidSolver::writeSolution(double currentTime){

  int rank = MASTER_NODE;
  int size = 1;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  //output->WriteHistory(solver, structure, &outputFile, currentTime);
  if(rank == MASTER_NODE){
  if(structure->GetnDof() == 1){
    //cout << "\"Time\"" << "\t" << "\"Displacement\"" << "\t" << "\"Velocity\"" << "\t" << "\"Acceleration\"" << "\t" << endl;
    //cout << currentTime << "\t" << (*(solver->GetDisp()))[0] << "\t" << (*(solver->GetVel()))[0] << "\t" << (*(solver->GetAcc()))[0] << endl;
    if(currentTime == 0) outputFile << "\"Time\"" << "\t" << "\"Displacement\"" << "\t" << "\"Velocity\"" << "\t" << "\"Acceleration\"" << "\t" << endl;
    outputFile << currentTime << "\t" << (*(solver->GetDisp()))[0] << "\t" << (*(solver->GetVel()))[0] << "\t" << (*(solver->GetAcc()))[0] << endl;
  }
  else if(structure->GetnDof() == 2){
    //cout << "\"Time\"" << "\t" << "\"Displacement 1\"" << "\t" << "\"Displacement 2\"" << "\t" << "\"Velocity 1\""  << "\t" << "\"Velocity 2\"" << "\t" << "\"Acceleration 1\"" << "\t" << "\"Acceleration 2\"" << endl;
    //cout << time << "\t" << (*(solver->GetDisp()))[0] << "\t" << (*(solver->GetDisp()))[1] << "\t" << (*(solver->GetVel()))[0] << "\t" << (*(solver->GetVel()))[1] << "\t" << (*(solver->GetAcc()))[0] << "\t" << (*(solver->GetAcc()))[1] << endl;
    if(currentTime == 0) outputFile << "\"Time\"" << "\t" << "\"Displacement 1\"" << "\t" << "\"Displacement 2\"" << "\t" << "\"Velocity 1\""  << "\t" << "\"Velocity 2\"" << "\t" << "\"Acceleration 1\"" << "\t" << "\"Acceleration 2\"" << endl;
    outputFile << currentTime << "\t" << (*(solver->GetDisp()))[0] << "\t" << (*(solver->GetDisp()))[1] << "\t" << (*(solver->GetVel()))[0] << "\t" << (*(solver->GetVel()))[1] << "\t" << (*(solver->GetAcc()))[0] << "\t" << (*(solver->GetAcc()))[1] << endl;
  }
  }


}

void NativeSolidSolver::updateSolution(){
  solver->UpdateSolution();
}

double NativeSolidSolver::outputDisplacements(){
  return (*(solver->GetDisp()))[0];
}

double NativeSolidSolver::displacementPredictor(){

  double deltaT, q_n, qdot_n, qdot_nM1, q_nP1, alpha0, alpha1;
  deltaT = config->GetDeltaT();
  q_n = (*(solver->GetDisp()))[0];
  qdot_n = (*(solver->GetVel()))[0];
  qdot_nM1 = (*(solver->GetVel_n()))[0];
  alpha0 = 1.0;
  alpha1 = 0.5;

  q_nP1 = q_n + alpha0*deltaT*qdot_n + alpha1*deltaT*(qdot_n - qdot_nM1);
}
