#include "NativeSolid_API.h"
#include <cmath>

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

    cout << endl << "\n----------------------- Setting setting FSI features ----------------------" << endl;
    q_uM1 = new CVector(structure->GetnDof());
    q_uM1->Reset();
    globalFluidLoads = new double[3];
    globalFluidLoads[0] = 0.0;
    globalFluidLoads[1] = 0.0;
    globalFluidLoads[2] = 0.0;

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
  delete q_uM1;
  delete [] globalFluidLoads;

}

void NativeSolidSolver::inputFluidLoads(double currentTime, double FSI_Load){

  solver->SetLoadsAtTime(config, structure, currentTime, FSI_Load);

}

double* NativeSolidSolver::getGlobalFluidLoadsArray() const{

  return globalFluidLoads;

}

void NativeSolidSolver::applyGlobalFluidLoads(){

  if(config->GetStructType() == "SPRING_HOR"){
    ((solver->GetLoads())->GetVec())[0] = globalFluidLoads[1];
  }
  else if(config->GetStructType() == "SPRING_VER"){
    ((solver->GetLoads())->GetVec())[0] = globalFluidLoads[0];
  }
  else if(config->GetStructType() == "AIRFOIL"){
    ((solver->GetLoads())->GetVec())[0] = -globalFluidLoads[0];
    ((solver->GetLoads())->GetVec())[1] = -globalFluidLoads[2];
  }
  else{
    cerr << "Wrong structural type for applying global fluild loads !" << endl;
    throw(-1);
  }

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

void NativeSolidSolver::staticComputation(){

  int rank = MASTER_NODE;
  int size = 1;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  solver->StaticIteration(config,structure);

  if(rank == MASTER_NODE){
    if(structure->GetnDof() == 1){
      cout << "Static Displacement : " << (*(solver->GetDisp()))[0] << endl;
    }
    else if(structure->GetnDof() == 2){
      cout << "Static Displacement 1 : " << (*(solver->GetDisp()))[0] << endl;
      cout << "Static Displacement 2 : " << (*(solver->GetDisp()))[1] << endl;
    }
  }

}

void NativeSolidSolver::writeSolution(double currentTime, double currentFSIIter){

  int rank = MASTER_NODE;
  int size = 1;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  if(rank == MASTER_NODE){
  if(structure->GetnDof() == 1){
    if(config->GetUnsteady() == "YES"){
      if(currentTime == 0) outputFile << "\"Time\"" << "\t" << "\"Displacement\"" << "\t" << "\"Velocity\"" << "\t" << "\"Acceleration\"" << "\t" << endl;
      outputFile << currentTime << "\t" << (*(solver->GetDisp()))[0] << "\t" << (*(solver->GetVel()))[0] << "\t" << (*(solver->GetAcc()))[0] << endl;
    }
    else{
      if(currentFSIIter == 0) outputFile << "\"FSI Iteration\"" << "\t" << "\"Displacement\"" << endl;
      outputFile << currentFSIIter << "\t" << (*(solver->GetDisp()))[0] << endl;
    }
  }
  else if(structure->GetnDof() == 2){
    if(config->GetUnsteady() == "YES"){
      if(currentTime == 0) outputFile << "\"Time\"" << "\t" << "\"Displacement 1\"" << "\t" << "\"Displacement 2\"" << "\t" << "\"Velocity 1\""  << "\t" << "\"Velocity 2\"" << "\t" << "\"Acceleration 1\"" << "\t" << "\"Acceleration 2\"" << endl;
      outputFile << currentTime << "\t" << (*(solver->GetDisp()))[0] << "\t" << (*(solver->GetDisp()))[1] << "\t" << (*(solver->GetVel()))[0] << "\t" << (*(solver->GetVel()))[1] << "\t" << (*(solver->GetAcc()))[0] << "\t" << (*(solver->GetAcc()))[1] << endl;
    }
    else{
      if(currentFSIIter == 0) outputFile << "\"FSI Iteration\"" << "\t" << "\"Displacement 1\"" << "\t" << "\"Displacement 1\"" << endl;
      outputFile << currentFSIIter << "\t" << (*(solver->GetDisp()))[0] << "\t" << (*(solver->GetDisp()))[1] << endl;
    }
  }
  }

}

void NativeSolidSolver::updateSolution(){

  if(config->GetUnsteady() == "YES")
    solver->UpdateSolution();
  else
    *q_uM1 = (*(solver->GetDisp()));
}

void NativeSolidSolver::outputDisplacements(double* interfRigidDispArray){

  double disp(0.0), dAlpha(0.0);

  interfRigidDispArray[0] = 0.0;
  interfRigidDispArray[1] = 0.0;
  interfRigidDispArray[2] = 0.0;
  interfRigidDispArray[3] = 0.0;
  interfRigidDispArray[4] = 0.0;
  interfRigidDispArray[5] = 0.0;

  if(config->GetUnsteady() == "YES")
    disp =  ( (*(solver->GetDisp()))[0] - (*(solver->GetDisp_n()))[0]);
    if (structure->GetnDof() == 2)
      dAlpha =  ( (*(solver->GetDisp()))[1] - (*(solver->GetDisp_n()))[1]);
  else
    disp =  ( (*(solver->GetDisp()))[0] - (*q_uM1)[0] );

  if(config->GetStructType() == "SPRING_HOR")
    interfRigidDispArray[0] = disp;
  else if (config->GetStructType() == "SPRING_VER")
    interfRigidDispArray[1] = disp;
  else if (config->GetStructType() == "AIRFOIL"){
    interfRigidDispArray[1] = -disp;
    interfRigidDispArray[5] = -dAlpha;
    cout << "Displacement communicated :" << endl;
    cout << interfRigidDispArray[1] << endl;
    cout << interfRigidDispArray[5] << endl;
  }
}

void NativeSolidSolver::displacementPredictor(double* interfRigidDispArray){

  double deltaT, q_n, qdot_n, qdot_nM1, q_nP1;
  double alpha_n, alphadot_n, alphadot_nM1, alpha_nP1;
  double alpha0, alpha1;
  double disp, dAlpha;

  deltaT = config->GetDeltaT();
  alpha0 = 1.0;
  alpha1 = 0.5; //Second order prediction

  q_n = (*(solver->GetDisp()))[0];
  qdot_n = (*(solver->GetVel()))[0];
  qdot_nM1 = (*(solver->GetVel_n()))[0];
  q_nP1 = q_n + alpha0*deltaT*qdot_n + alpha1*deltaT*(qdot_n - qdot_nM1);
  disp = q_nP1-q_n;

  if(structure->GetnDof() == 2){
    alpha_n = (*(solver->GetDisp()))[1];
    alphadot_n = (*(solver->GetVel()))[1];
    alphadot_nM1 = (*(solver->GetVel_n()))[1];
    alpha_nP1 = alpha_n + alpha0*deltaT*alphadot_n + alpha1*deltaT*(alphadot_n - alphadot_nM1);
    dAlpha = alpha_nP1-alpha_n;
  }

  if(config->GetStructType() == "SPRING_HOR")
    interfRigidDispArray[0] = disp;
  else if (config->GetStructType() == "SPRING_VER")
    interfRigidDispArray[1] = disp;
  else if (config->GetStructType() == "AIRFOIL"){
    interfRigidDispArray[1] = -disp;
    interfRigidDispArray[5] = -dAlpha;
  }
}

void NativeSolidSolver::setAitkenCoefficient(unsigned long FSIIter){

  double u_uM1, u_pred, u_rel, q_pred, q_rel, q_nM1, lambda(0.0);

  q_pred = (*(solver->GetDisp()))[0];

  if(config->GetUnsteady() == "YES"){


  q_nM1 = (*(solver->GetDisp_n()))[0];
  u_pred = q_pred - q_nM1;

  if(FSIIter == 0){
    omega = 1.0;
    (*q_uM1)[0] = q_pred;
    q_rel = q_pred;

  }
  else{
    u_uM1 = (*q_uM1)[0] - q_nM1;

    lambda = - (u_uM1)*(u_pred)/(pow(abs(u_pred - u_uM1),2));
    omega *= lambda;
    cout << omega << endl;
    cout << u_uM1 << endl;

    u_rel = omega*u_pred + (1-omega)*u_uM1;
    u_uM1 = u_rel;

    (*q_uM1)[0] = u_uM1 + q_nM1;
    q_rel = u_rel + q_nM1;
  }

  }
  else{
    if(FSIIter == 0){
      omega = 1.0;
      q_rel = q_pred;
    }
    else{

      lambda = -((*q_uM1)[0])*(q_pred)/(pow(abs(q_pred-(*q_uM1)[0]),2));
      omega *= lambda;
      cout << omega << endl;
      cout << u_uM1 << endl;

      q_rel = omega*q_pred + (1-omega)*((*q_uM1)[0]);
    }
  }

  solver->SetDisp(q_rel);

}

double NativeSolidSolver::aitkenRelaxation(){

}
