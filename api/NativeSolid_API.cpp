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

  historyFile.open("NativeHistory.dat", ios::out);

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

    /*--- Read a SU2 native mesh file ---*/
    cout << endl << "\n----------------------- Reading SU2 based mesh file ----------------------" << endl;
    geometry = new Geometry(config);

    /*--- Initialize structural container and create the structural model ---*/
    cout << endl << "\n----------------------- Creating the structural model ----------------------" << endl;
    structure = new Structure(config);
    structure->SetStructuralMatrices(config);

    cout << endl << "\n----------------------- Creating the FSI interface ----------------------" << endl;
    double* Coord;
    unsigned long iMarker(0), iPoint;

    while(iMarker < geometry->GetnMarkers()){
      if(geometry->GetMarkersMoving(iMarker)){
        nSolidInterfaceVertex = geometry->nVertex[iMarker];
        solidInterface = new double*[nSolidInterfaceVertex];
        for(int iVertex=0; iVertex<nSolidInterfaceVertex; iVertex++){
          solidInterface[iVertex] = new double[4];
          solidInterface[iVertex][0] = geometry->vertex[iMarker][iVertex];
          iPoint = solidInterface[iVertex][0];
          Coord = geometry->node[iPoint]->GetCoord();
          solidInterface[iVertex][1] = Coord[0];
          solidInterface[iVertex][2] = Coord[1];
          solidInterface[iVertex][3] = Coord[2];
        }
      break;
      }
      iMarker++;
    }

    cout << nSolidInterfaceVertex << " nodes on the moving interface have to be tracked." << endl;

    /*for(int iVertex = 0; iVertex<nSolidInterfaceVertex; iVertex++){
      cout << solidInterface[iVertex][0] << ";" << solidInterface[iVertex][1] << ";" << solidInterface[iVertex][2] << ";" << solidInterface[iVertex][3] << endl;
    }*/

    /*--- Initialize the temporal integrator ---*/
    cout << endl << "\n----------------------- Setting integration parameter ----------------------" << endl;
    solver = new Integration(structure);
    if(config->GetUnsteady() == "YES"){
      solver->SetIntegrationParam(config);
      solver->SetInitialConditions(config, structure);
      mapRigidBodyMotion(false, true);          //The mesh is already set in the solid side
    }

    cout << endl << "\n----------------------- Setting FSI features ----------------------" << endl;
    q_uM1 = new CVector(structure->GetnDof());
    q_uM1->Reset();
    globalFluidLoads = new double[3];
    globalFluidLoads[0] = 0.0;
    globalFluidLoads[1] = 0.0;
    globalFluidLoads[2] = 0.0;

    if(rank == MASTER_NODE){
      if(structure->GetnDof() == 1){
        if(config->GetUnsteady() == "YES"){
          historyFile << "\"Time\"" << "\t" << "\"Displacement\"" << "\t" << "\"Velocity\"" << "\t" << "\"Acceleration\"" << "\t" << "\"Acceleration variable\"" << endl;
        }
        else{
          historyFile << "\"FSI Iteration\"" << "\t" << "\"Displacement\"" << endl;
        }
      }
      else if(structure->GetnDof() == 2){
        if(config->GetUnsteady() == "YES"){
          historyFile << "\"Time\"" << "\t" << "\"Displacement 1\"" << "\t" << "\"Displacement 2\"" << "\t" << "\"Velocity 1\""  << "\t" << "\"Velocity 2\"" << "\t" << "\"Acceleration 1\"" << "\t" << "\"Acceleration 2\"" << "\t" << "\"Acceleration variable 1\"" << "\t" << "\"Acceleration variable 2\"" << endl;
        }
        else{
          historyFile << "\"FSI Iteration\"" << "\t" << "\"Displacement 1\"" << "\t" << "\"Displacement 1\"" << endl;
        }
      }
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

  historyFile.close();
  if(rank == MASTER_NODE){
    cout << "Solid history file is closed." << endl;
    cout << endl << "***************************** Exit NativeSolid *****************************" << endl;
  }

  delete config;
  delete geometry;
  delete structure;
  delete solver;
  delete output;
  delete q_uM1;
  delete [] globalFluidLoads;

  if(solidInterface != NULL){
    for(int iVertex=0; iVertex<nSolidInterfaceVertex; iVertex++){
      delete [] solidInterface[iVertex];
    }

    delete [] solidInterface;
  }

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

  mapRigidBodyMotion(false, false);

}

void NativeSolidSolver::mapRigidBodyMotion(bool prediction, bool initialize){

  double *Coord, *Coord_n, newCoord[3], Center[3], Center_n[3], newCenter[3], rotCoord[3], r[3], varCoord[3];
  double rotMatrix[3][3] = {{0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
  double dTheta, dPhi, dPsi;
  double cosTheta, sinTheta, cosPhi, sinPhi, cosPsi, sinPsi;
  unsigned short iMarker, iVertex, nDim(3);
  unsigned long iPoint;

  double deltaT, q_n, qdot_n, qdot_nM1, q_nP1;
  double alpha_n, alphadot_n, alphadot_nM1, alpha_nP1;
  double alpha0, alpha1;
  double disp, dAlpha;

  /*--- Get the current center of rotation (can vary at each iteration) ---*/
  Center[0] = structure->GetCenterOfRotation_x();
  Center[1] = structure->GetCenterOfRotation_y();
  Center[2] = structure->GetCenterOfRotation_z();

  /*--- Get the center of rotation from previous time step ---*/
  Center_n[0] = structure->GetCenterOfRotation_n_x();
  Center_n[1] = structure->GetCenterOfRotation_n_y();
  Center_n[2] = structure->GetCenterOfRotation_n_z();

  if(!prediction){
    dTheta = 0.0;
    dPhi = 0.0;
    if (config->GetStructType() == "AIRFOIL"){
      dPsi = -( (*(solver->GetDisp()))[1] - (*(solver->GetDisp_n()))[1]);
      newCenter[0] = Center[0];
      newCenter[1] = -(*(solver->GetDisp()))[0];
      newCenter[2] = Center[2];
      disp = newCenter[1] - Center_n[1];
    }
    else if(config->GetStructType() == "SPRING_HOR"){
      dPsi = 0.0;
      newCenter[0] = (*(solver->GetDisp()))[0];
      newCenter[1] = Center[1];
      newCenter[2] = Center[2];
    }

    //cout << "Center of rotation at previous time step : " << Center_n[0] << " " << Center_n[1] << " " << Center_n[2] << endl;
    //cout << "Current center of rotation : " << Center[0] << " " << Center[1] << " " << Center[2] << endl;



    //cout << "New center of rotation : " << newCenter[0] << " " << newCenter[1] << " " << newCenter[2] << endl;
  }
  else{
    deltaT = config->GetDeltaT();
    alpha0 = 1.0;
    alpha1 = 0.5; //Second order prediction

    q_n = (*(solver->GetDisp()))[0];
    qdot_n = (*(solver->GetVel()))[0];
    qdot_nM1 = (*(solver->GetVel_n()))[0];
    q_nP1 = q_n + alpha0*deltaT*qdot_n + alpha1*deltaT*(qdot_n - qdot_nM1);
    //cout << "Predicted position q : " << q_nP1 << endl;
    disp = q_nP1-q_n;

    if(structure->GetnDof() == 2){
      alpha_n = (*(solver->GetDisp()))[1];
      alphadot_n = (*(solver->GetVel()))[1];
      alphadot_nM1 = (*(solver->GetVel_n()))[1];
      alpha_nP1 = alpha_n + alpha0*deltaT*alphadot_n + alpha1*deltaT*(alphadot_n - alphadot_nM1);
      //cout << "Predicted position alpha : " << alpha_nP1 << endl;
      dAlpha = alpha_nP1-alpha_n;
    }
    dTheta = 0.0;
    dPhi = 0.0;
    if (config->GetStructType() == "AIRFOIL"){
      dPsi = -dAlpha;
      newCenter[0] = Center[0];
      newCenter[1] = -q_nP1;
      newCenter[2] = Center[2];
    }
    else{
      dPsi = 0.0;
    }
  }

  /*if (config->GetStructType() == "AIRFOIL"){
    cout << "Displacement communicated :" << endl;
    cout << disp << endl;
    cout << dPsi << endl;
  }*/


  cosTheta = cos(dTheta);  cosPhi = cos(dPhi);  cosPsi = cos(dPsi);
  sinTheta = sin(dTheta);  sinPhi = sin(dPhi);  sinPsi = sin(dPsi);

  /*--- Compute the rotation matrix. The implicit
  ordering is rotation about the x-axis, y-axis, then z-axis. ---*/

  rotMatrix[0][0] = cosPhi*cosPsi;
  rotMatrix[1][0] = cosPhi*sinPsi;
  rotMatrix[2][0] = -sinPhi;

  rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
  rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
  rotMatrix[2][1] = sinTheta*cosPhi;

  rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
  rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
  rotMatrix[2][2] = cosTheta*cosPhi;

  for(iMarker = 0; iMarker < geometry->GetnMarkers(); iMarker++){
    if (geometry->markersMoving[iMarker] == true){
      for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++){
        iPoint = geometry->vertex[iMarker][iVertex];
        Coord = geometry->node[iPoint]->GetCoord();
        Coord_n = geometry->node[iPoint]->GetCoord_n();

        if(config->GetUnsteady() == "YES"){
          for (int iDim=0; iDim < nDim; iDim++){
            r[iDim] = Coord_n[iDim] - Center_n[iDim];
          }
        }
        else{
          for (int iDim=0; iDim < nDim; iDim++){
            r[iDim] = Coord[iDim] - Center[iDim];
          }
        }
        //cout << Center[2] << endl;
        //cout << Coord[2] << endl;
        //if((r[0]*r[0]+r[1]*r[1]+r[2]*r[2])==0) cout << "PROBLEEEEEEEME" << endl;

        rotCoord[0] = rotMatrix[0][0]*r[0]
              + rotMatrix[0][1]*r[1]
              + rotMatrix[0][2]*r[2];

        rotCoord[1] = rotMatrix[1][0]*r[0]
              + rotMatrix[1][1]*r[1]
              + rotMatrix[1][2]*r[2];

        rotCoord[2] = rotMatrix[2][0]*r[0]
              + rotMatrix[2][1]*r[1]
              + rotMatrix[2][2]*r[2];

              //cout << rotCoord[2] << endl;
        if(solidInterface[iVertex][0] == 0) cout << "VarCoord : " << endl;
        for(int iDim=0; iDim < nDim; iDim++){
                newCoord[iDim] = newCenter[iDim] + rotCoord[iDim];
                solidInterface[iVertex][iDim+1] = newCoord[iDim];
                varCoord[iDim] = newCoord[iDim] - Coord[iDim];
                if(solidInterface[iVertex][0] == 0)cout << varCoord[iDim] << endl;
        }

        /*--- Apply change of coordinates to the node on the moving interface ---*/
        //if(!prediction) geometry->node[iPoint]->SetCoord(newCoord);
        geometry->node[iPoint]->SetCoord(newCoord);
        /*--- At initialisation, propagate the initial position of the inteface in the past ---*/
        if(initialize) geometry->node[iPoint]->SetCoord_n(newCoord);
        //cout << "**********************************************" << endl;
        //geometry->node[iPoint]->PrintCoord();
        //cout << solidInterface[iVertex][0] << ";" << solidInterface[iVertex][1] << ";" << solidInterface[iVertex][2] << ";" << solidInterface[iVertex][3] << endl;
        //solidInterface[iVertex][2] = Coord[1];
        //solidInterface[iVertex][3] = Coord[2];

      }
    }
  }

  /*--- Update the position of the center of rotation ---*/
  //if(!prediction){
  structure->SetCenterOfRotation_X(newCenter[0]);
  structure->SetCenterOfRotation_Y(newCenter[1]);
  structure->SetCenterOfRotation_Z(newCenter[2]);
  //}

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

  mapRigidBodyMotion(false,false);

}

void NativeSolidSolver::writeSolution(double currentTime, double currentFSIIter, unsigned long ExtIter, unsigned long NbExtIter){

  int rank = MASTER_NODE;
  int size = 1;
  double DeltaT;
  string restartFileName;
  unsigned long DeltaIter = config->GetDeltaIterWrite();

  DeltaT = config->GetDeltaT();

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  if(rank == MASTER_NODE){
  if(structure->GetnDof() == 1){
    if(config->GetUnsteady() == "YES"){
      //if(currentTime == config->GetStartTime()) historyFile << "\"Time\"" << "\t" << "\"Displacement\"" << "\t" << "\"Velocity\"" << "\t" << "\"Acceleration\"" << "\t" << "\"Acceleration variable\"" << endl;
      historyFile << currentTime << "\t" << (*(solver->GetDisp()))[0] << "\t" << (*(solver->GetVel()))[0] << "\t" << (*(solver->GetAcc()))[0] << "\t" << (*(solver->GetAccVar()))[0] << endl;
    }
    else{
      //if(currentFSIIter == config->GetStartTime()) historyFile << "\"FSI Iteration\"" << "\t" << "\"Displacement\"" << endl;
      historyFile << currentFSIIter << "\t" << (*(solver->GetDisp()))[0] << endl;
    }
  }
  else if(structure->GetnDof() == 2){
    if(config->GetUnsteady() == "YES"){
      //if(currentTime == config->GetStartTime()) historyFile << "\"Time\"" << "\t" << "\"Displacement 1\"" << "\t" << "\"Displacement 2\"" << "\t" << "\"Velocity 1\""  << "\t" << "\"Velocity 2\"" << "\t" << "\"Acceleration 1\"" << "\t" << "\"Acceleration 2\"" << "\t" << "\"Acceleration variable 1\"" << "\t" << "\"Acceleration variable 2\"" << endl;
      historyFile << currentTime << "\t" << (*(solver->GetDisp()))[0] << "\t" << (*(solver->GetDisp()))[1] << "\t" << (*(solver->GetVel()))[0] << "\t" << (*(solver->GetVel()))[1] << "\t" << (*(solver->GetAcc()))[0] << "\t" << (*(solver->GetAcc()))[1] << "\t" << (*(solver->GetAccVar()))[0] << "\t" << (*(solver->GetAccVar()))[1] << endl;
    }
    else{
      //if(currentFSIIter == config->GetStartTime()) historyFile << "\"FSI Iteration\"" << "\t" << "\"Displacement 1\"" << "\t" << "\"Displacement 1\"" << endl;
      historyFile << currentFSIIter << "\t" << (*(solver->GetDisp()))[0] << "\t" << (*(solver->GetDisp()))[1] << endl;
    }
  }
  }


  if(config->GetUnsteady() == "YES"){
  if ( ExtIter%DeltaIter == 0 || ExtIter == NbExtIter ){
    restartFileName = config->GetRestartFile();
    restartFile.open(restartFileName.c_str(), ios::out);
    if (structure->GetnDof() == 1){
      restartFile << "\"Time\"" << "\t" << "\"Displacement\"" << "\t" << "\"Velocity\"" << "\t" << "\"Acceleration\"" << "\t" << "\"Acceleration variable\"" << endl;
      restartFile << currentTime-DeltaT << "\t" << (*(solver->GetDisp_n()))[0] << "\t" << (*(solver->GetVel_n()))[0] << "\t" << (*(solver->GetAcc_n()))[0] << "\t" << (*(solver->GetAccVar_n()))[0] << endl;
      restartFile << currentTime << "\t" << (*(solver->GetDisp()))[0] << "\t" << (*(solver->GetVel()))[0] << "\t" << (*(solver->GetAcc()))[0] << "\t" << (*(solver->GetAccVar()))[0] << endl;
    }
    else if (structure->GetnDof() == 2){
      restartFile << "\"Time\"" << "\t" << "\"Displacement 1\"" << "\t" << "\"Displacement 2\"" << "\t" << "\"Velocity 1\""  << "\t" << "\"Velocity 2\"" << "\t" << "\"Acceleration 1\"" << "\t" << "\"Acceleration 2\"" << "\t" << "\"Acceleration variable 1\"" << "\t" << "\"Acceleration variable 2\"" << endl;
      restartFile << currentTime-DeltaT << "\t" << (*(solver->GetDisp_n()))[0] << "\t" << (*(solver->GetDisp_n()))[1] << "\t" << (*(solver->GetVel_n()))[0] << "\t" << (*(solver->GetVel_n()))[1] << "\t" << (*(solver->GetAcc_n()))[0] << "\t" << (*(solver->GetAcc_n()))[1] << "\t" << (*(solver->GetAccVar_n()))[0] << "\t" << (*(solver->GetAccVar_n()))[1] << endl;
      restartFile << currentTime << "\t" << (*(solver->GetDisp()))[0] << "\t" << (*(solver->GetDisp()))[1] << "\t" << (*(solver->GetVel()))[0] << "\t" << (*(solver->GetVel()))[1] << "\t" << (*(solver->GetAcc()))[0] << "\t" << (*(solver->GetAcc()))[1] << "\t" << (*(solver->GetAccVar()))[0] << "\t" << (*(solver->GetAccVar()))[1] << endl;
    }
    restartFile.close();
  }
  }

}

void NativeSolidSolver::updateSolution(){

  if(config->GetUnsteady() == "YES"){
    solver->UpdateSolution();
  }
  else
    *q_uM1 = (*(solver->GetDisp()));
}

void NativeSolidSolver::updateGeometry(){

  if(config->GetUnsteady() == "YES"){
    geometry->UpdateGeometry();
    structure->SetCenterOfRotation_n_X(structure->GetCenterOfRotation_x());
    structure->SetCenterOfRotation_n_Y(structure->GetCenterOfRotation_y());
    structure->SetCenterOfRotation_n_Z(structure->GetCenterOfRotation_z());
  }
}

void NativeSolidSolver::outputDisplacements(double* interfRigidDispArray, bool initialize){

  double disp(0.0), dAlpha(0.0);

  interfRigidDispArray[0] = 0.0;
  interfRigidDispArray[1] = 0.0;
  interfRigidDispArray[2] = 0.0;
  interfRigidDispArray[3] = 0.0;
  interfRigidDispArray[4] = 0.0;
  interfRigidDispArray[5] = 0.0;

  if(initialize){
    if(config->GetUnsteady() == "YES"){
      disp =  (*(solver->GetDisp()))[0];
      if (structure->GetnDof() == 2)
        dAlpha = (*(solver->GetDisp()))[1];
    }
    else{
      disp = (*(solver->GetDisp()))[0];
    }
  }
  else{
    if(config->GetUnsteady() == "YES"){
      disp =  ( (*(solver->GetDisp()))[0] - (*(solver->GetDisp_n()))[0]);
      if (structure->GetnDof() == 2)
        dAlpha =  ( (*(solver->GetDisp()))[1] - (*(solver->GetDisp_n()))[1]);
    }
    else{
      disp =  ( (*(solver->GetDisp()))[0] - (*q_uM1)[0] );
    }
  }

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

void NativeSolidSolver::displacementPredictor_Old(double* interfRigidDispArray){

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
  cout << "Predicted position q : " << q_nP1 << endl;

  if(structure->GetnDof() == 2){
    alpha_n = (*(solver->GetDisp()))[1];
    alphadot_n = (*(solver->GetVel()))[1];
    alphadot_nM1 = (*(solver->GetVel_n()))[1];
    alpha_nP1 = alpha_n + alpha0*deltaT*alphadot_n + alpha1*deltaT*(alphadot_n - alphadot_nM1);
    cout << "Predicted position alpha : " << alpha_nP1 << endl;
    dAlpha = alpha_nP1-alpha_n;
  }

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

void NativeSolidSolver::displacementPredictor(){

  mapRigidBodyMotion(true, false);

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

double** NativeSolidSolver::getSolidInterface() const{
  return solidInterface;
}

const double* NativeSolidSolver::getCenterCoordinate() const{
  return structure->GetCenterOfRotation();
}

unsigned long NativeSolidSolver::getnSolidInterfaceVertex() const{
  return nSolidInterfaceVertex;
}
