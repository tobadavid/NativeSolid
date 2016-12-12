#include "NativeSolid_API.h"
#include <cmath>

using namespace std;

const int MASTER_NODE = 0;


NativeSolidSolver::NativeSolidSolver(string str, bool FSIComp):confFile(str){

int rank = MASTER_NODE;
  int size = 1;

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
      break;
      }
      iMarker++;
    }

    varCoordNorm = 0.0;

    cout << nSolidInterfaceVertex << " nodes on the moving interface have to be tracked." << endl;

    /*--- Initialize the temporal integrator ---*/
    cout << endl << "\n----------------------- Setting integration parameter ----------------------" << endl;
    solver = new Integration(structure);
    if(config->GetUnsteady() == "YES"){
      solver->SetIntegrationParam(config);
      solver->SetInitialConditions(config, structure);
    }

    cout << endl << "\n----------------------- Setting FSI features ----------------------" << endl;
    q_uM1 = new CVector(structure->GetnDof());
    q_uM1->Reset();

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

NativeSolidSolver::~NativeSolidSolver(){}

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

}

void NativeSolidSolver::preprocessIteration(unsigned long ExtIter){

    solver->SetExtIter(ExtIter);
}

/*double NativeSolidSolver::getVarCoordNorm() const{

  return varCoordNorm;

}*/

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

  //mapRigidBodyMotion(false, false);
  computeInterfacePosVel(false);

}

void NativeSolidSolver::mapRigidBodyMotion(bool prediction, bool initialize){

  double *Coord, *Coord_n, newCoord[3], Center[3], Center_n[3], newCenter[3], rotCoord[3], r[3];
  double varCoord[3] = {0.0, 0.0, 0.0};
  double rotMatrix[3][3] = {{0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
  double dTheta, dPhi, dPsi;
  double cosTheta, sinTheta, cosPhi, sinPhi, cosPsi, sinPsi;
  unsigned short iMarker, iVertex, nDim(3);
  unsigned long iPoint;

  double deltaT, q_n, qdot_n, qdot_nM1, q_nP1;
  double alpha_n, alphadot_n, alphadot_nM1, alpha_nP1;
  double alpha0, alpha1;
  double disp, dAlpha;
  double varCoordNorm2(0.0);

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
  }
  else{
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

        rotCoord[0] = rotMatrix[0][0]*r[0]
              + rotMatrix[0][1]*r[1]
              + rotMatrix[0][2]*r[2];

        rotCoord[1] = rotMatrix[1][0]*r[0]
              + rotMatrix[1][1]*r[1]
              + rotMatrix[1][2]*r[2];

        rotCoord[2] = rotMatrix[2][0]*r[0]
              + rotMatrix[2][1]*r[1]
              + rotMatrix[2][2]*r[2];

        for(int iDim=0; iDim < nDim; iDim++){
                newCoord[iDim] = newCenter[iDim] + rotCoord[iDim];
                varCoord[iDim] = newCoord[iDim] - Coord[iDim];
        }

        varCoordNorm2 += varCoord[0]*varCoord[0] + varCoord[1]*varCoord[1] + varCoord[2]*varCoord[2];

        /*--- Apply change of coordinates to the node on the moving interface ---*/
        geometry->node[iPoint]->SetCoord(newCoord);

        /*--- At initialisation, propagate the initial position of the inteface in the past ---*/
        if(initialize) geometry->node[iPoint]->SetCoord_n(newCoord);
      }
    }
  }

  varCoordNorm = sqrt(varCoordNorm2);

  /*--- Update the position of the center of rotation ---*/
  structure->SetCenterOfRotation_X(newCenter[0]);
  structure->SetCenterOfRotation_Y(newCenter[1]);
  structure->SetCenterOfRotation_Z(newCenter[2]);
}

void NativeSolidSolver::computeInterfacePosVel(bool initialize){

    double *Coord, *Coord_n, newCoord[3], newVel[3], Center[3], Center_n[3], newCenter[3], rotCoord[3], r[3];
    double varCoord[3] = {0.0, 0.0, 0.0};
    double rotMatrix[3][3] = {{0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
    double dTheta, dPhi, dPsi;
    double psidot;
    double cosTheta, sinTheta, cosPhi, sinPhi, cosPsi, sinPsi;
    unsigned short iMarker, iVertex, nDim(3);
    unsigned long iPoint;
    double varCoordNorm2(0.0);

    /*--- Get the current center of rotation (can vary at each iteration) ---*/
    Center[0] = structure->GetCenterOfRotation_x();
    Center[1] = structure->GetCenterOfRotation_y();
    Center[2] = structure->GetCenterOfRotation_z();

    /*--- Get the center of rotation from previous time step ---*/
    Center_n[0] = structure->GetCenterOfRotation_n_x();
    Center_n[1] = structure->GetCenterOfRotation_n_y();
    Center_n[2] = structure->GetCenterOfRotation_n_z();

    dTheta = 0.0;
    dPhi = 0.0;
    if (config->GetStructType() == "AIRFOIL"){
      dPsi = -( (*(solver->GetDisp()))[1] - (*(solver->GetDisp_n()))[1]);
      psidot = (*(solver->GetVel()))[1];
      newCenter[0] = Center[0];
      newCenter[1] = -(*(solver->GetDisp()))[0];
      newCenter[2] = Center[2];
    }
    else if(config->GetStructType() == "SPRING_HOR"){
      dPsi = 0.0;
      newCenter[0] = (*(solver->GetDisp()))[0];
      newCenter[1] = Center[1];
      newCenter[2] = Center[2];
    }

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

          rotCoord[0] = rotMatrix[0][0]*r[0]
                + rotMatrix[0][1]*r[1]
                + rotMatrix[0][2]*r[2];

          rotCoord[1] = rotMatrix[1][0]*r[0]
                + rotMatrix[1][1]*r[1]
                + rotMatrix[1][2]*r[2];

          rotCoord[2] = rotMatrix[2][0]*r[0]
                + rotMatrix[2][1]*r[1]
                + rotMatrix[2][2]*r[2];

          for(int iDim=0; iDim < nDim; iDim++){
                  newCoord[iDim] = newCenter[iDim] + rotCoord[iDim];
                  varCoord[iDim] = newCoord[iDim] - Coord[iDim];
          }
          newVel[0] = psidot*(newCoord[1]-newCenter[1]);
          newVel[1] = -psidot*(newCoord[0]-newCenter[0]);
          newVel[2] = 0.0;

          varCoordNorm2 += varCoord[0]*varCoord[0] + varCoord[1]*varCoord[1] + varCoord[2]*varCoord[2];

          /*--- Apply change of coordinates to the node on the moving interface ---*/
          geometry->node[iPoint]->SetCoord(newCoord);
          geometry->node[iPoint]->SetVel(newVel);

          /*--- At initialisation, propagate the initial position of the inteface in the past ---*/
          if(initialize){
              geometry->node[iPoint]->SetCoord_n(newCoord);
              geometry->node[iPoint]->SetVel_n(newVel);
          }
        }
      }
    }

    varCoordNorm = sqrt(varCoordNorm2);

    /*--- Update the position of the center of rotation ---*/
    structure->SetCenterOfRotation_X(newCenter[0]);
    structure->SetCenterOfRotation_Y(newCenter[1]);
    structure->SetCenterOfRotation_Z(newCenter[2]);

}

void NativeSolidSolver::setInitialDisplacements(){
  //mapRigidBodyMotion(false, true);
  computeInterfacePosVel(true);
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

  //mapRigidBodyMotion(false,false);
  computeInterfacePosVel(false);

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

void NativeSolidSolver::saveSolution(){

    int rank = MASTER_NODE;
    int size = 1;
    double DeltaT;
    //string restartFileName;
    //unsigned long DeltaIter = config->GetDeltaIterWrite();
    unsigned long ExtIter = solver->GetExtIter();

    DeltaT = config->GetDeltaT();

  #ifdef HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
  #endif

    if(rank == MASTER_NODE){
    if(structure->GetnDof() == 1){
      if(config->GetUnsteady() == "YES"){
        historyFile << ExtIter << "\t" << (*(solver->GetDisp()))[0] << "\t" << (*(solver->GetVel()))[0] << "\t" << (*(solver->GetAcc()))[0] << "\t" << (*(solver->GetAccVar()))[0] << endl;
      }
      else{
        historyFile << ExtIter << "\t" << (*(solver->GetDisp()))[0] << endl;
      }
    }
    else if(structure->GetnDof() == 2){
      if(config->GetUnsteady() == "YES"){
        historyFile << ExtIter << "\t" << (*(solver->GetDisp()))[0] << "\t" << (*(solver->GetDisp()))[1] << "\t" << (*(solver->GetVel()))[0] << "\t" << (*(solver->GetVel()))[1] << "\t" << (*(solver->GetAcc()))[0] << "\t" << (*(solver->GetAcc()))[1] << "\t" << (*(solver->GetAccVar()))[0] << "\t" << (*(solver->GetAccVar()))[1] << endl;
      }
      else{
        historyFile << ExtIter << "\t" << (*(solver->GetDisp()))[0] << "\t" << (*(solver->GetDisp()))[1] << endl;
      }
    }
    }
}

void NativeSolidSolver::updateSolution(){

  if(config->GetUnsteady() == "YES"){
    solver->UpdateSolution();
    geometry->UpdateGeometry();
    structure->SetCenterOfRotation_n_X(structure->GetCenterOfRotation_x());
    structure->SetCenterOfRotation_n_Y(structure->GetCenterOfRotation_y());
    structure->SetCenterOfRotation_n_Z(structure->GetCenterOfRotation_z());
  }
  else
    *q_uM1 = (*(solver->GetDisp()));
}

/*void NativeSolidSolver::updateGeometry(){

  if(config->GetUnsteady() == "YES"){
    geometry->UpdateGeometry();
    structure->SetCenterOfRotation_n_X(structure->GetCenterOfRotation_x());
    structure->SetCenterOfRotation_n_Y(structure->GetCenterOfRotation_y());
    structure->SetCenterOfRotation_n_Z(structure->GetCenterOfRotation_z());
  }
}*/

/*void NativeSolidSolver::displacementPredictor(){

  mapRigidBodyMotion(true, false);

}*/

unsigned short NativeSolidSolver::getFSIMarkerID(){

  unsigned short iMarker(0), IDtoSend;


  while(iMarker < geometry->GetnMarkers()){
      if(geometry->GetMarkersMoving(iMarker)){
          IDtoSend =  iMarker;
          break;
      }
      iMarker++;
  }
  return IDtoSend;
}

unsigned long NativeSolidSolver::getNumberOfSolidInterfaceNodes(unsigned short iMarker){
  return nSolidInterfaceVertex;
}

unsigned int NativeSolidSolver::getInterfaceNodeGlobalIndex(unsigned short iMarker, unsigned short iVertex){

  unsigned long iPoint;

  iPoint = geometry->vertex[iMarker][iVertex];

  return iPoint;
}

double NativeSolidSolver::getInterfaceNodePosX(unsigned short iMarker, unsigned short iVertex){

    unsigned long iPoint;
    double *Coord;

    iPoint = geometry->vertex[iMarker][iVertex];
    Coord = geometry->node[iPoint]->GetCoord();

    return Coord[0];
}

double NativeSolidSolver::getInterfaceNodePosY(unsigned short iMarker, unsigned short iVertex){

    unsigned long iPoint;
    double *Coord;

    iPoint = geometry->vertex[iMarker][iVertex];
    Coord = geometry->node[iPoint]->GetCoord();

    return Coord[1];
}

double NativeSolidSolver::getInterfaceNodePosZ(unsigned short iMarker, unsigned short iVertex){

    unsigned long iPoint;
    double *Coord;

    iPoint = geometry->vertex[iMarker][iVertex];
    Coord = geometry->node[iPoint]->GetCoord();

    return 0.0; //3D is not really implemented in this solver...
}

double NativeSolidSolver::getInterfaceNodePosX0(unsigned short iMarker, unsigned short iVertex){

    unsigned long iPoint;
    double *Coord;

    iPoint = geometry->vertex[iMarker][iVertex];
    Coord = geometry->node[iPoint]->GetCoord0();

    return Coord[0];
}

double NativeSolidSolver::getInterfaceNodePosY0(unsigned short iMarker, unsigned short iVertex){

    unsigned long iPoint;
    double *Coord;

    iPoint = geometry->vertex[iMarker][iVertex];
    Coord = geometry->node[iPoint]->GetCoord0();

    return Coord[1];
}

double NativeSolidSolver::getInterfaceNodePosZ0(unsigned short iMarker, unsigned short iVertex){

    unsigned long iPoint;
    double *Coord;

    iPoint = geometry->vertex[iMarker][iVertex];
    Coord = geometry->node[iPoint]->GetCoord0();

    return 0.0;
}

double NativeSolidSolver::getInterfaceNodeDispX(unsigned short iMarker, unsigned short iVertex){

    unsigned long iPoint;
    double *Coord, *Coord0;

    iPoint = geometry->vertex[iMarker][iVertex];
    Coord = geometry->node[iPoint]->GetCoord();
    Coord0 = geometry->node[iPoint]->GetCoord0();

    return Coord[0] - Coord0[0];
}

double NativeSolidSolver::getInterfaceNodeDispY(unsigned short iMarker, unsigned short iVertex){

    unsigned long iPoint;
    double *Coord, *Coord0;

    iPoint = geometry->vertex[iMarker][iVertex];
    Coord = geometry->node[iPoint]->GetCoord();
    Coord0 = geometry->node[iPoint]->GetCoord0();
    return Coord[1] - Coord0[1];
}

double NativeSolidSolver::getInterfaceNodeDispZ(unsigned short iMarker, unsigned short iVertex){

    unsigned long iPoint;
    double *Coord, *Coord0;

    iPoint = geometry->vertex[iMarker][iVertex];
    Coord = geometry->node[iPoint]->GetCoord();
    Coord0 = geometry->node[iPoint]->GetCoord0();

    return Coord[2] - Coord0[2];
}

double NativeSolidSolver::getInterfaceNodeVelX(unsigned short iMarker, unsigned short iVertex){
    unsigned long iPoint;
    double *Vel;

    iPoint = geometry->vertex[iMarker][iVertex];
    Vel = geometry->node[iPoint]->GetVel();

    return Vel[0];
}

double NativeSolidSolver::getInterfaceNodeVelY(unsigned short iMarker, unsigned short iVertex){
    unsigned long iPoint;
    double *Vel;

    iPoint = geometry->vertex[iMarker][iVertex];
    Vel = geometry->node[iPoint]->GetVel();

    return Vel[1];
}

double NativeSolidSolver::getInterfaceNodeVelZ(unsigned short iMarker, unsigned short iVertex){
    unsigned long iPoint;
    double *Vel;

    iPoint = geometry->vertex[iMarker][iVertex];
    Vel = geometry->node[iPoint]->GetVel();

    return Vel[2];
}

double NativeSolidSolver::getInterfaceNodeVelXNm1(unsigned short iMarker, unsigned short iVertex){
    unsigned long iPoint;
    double *Vel_n;

    iPoint = geometry->vertex[iMarker][iVertex];
    Vel_n = geometry->node[iPoint]->GetVel_n();

    return Vel_n[0];
}

double NativeSolidSolver::getInterfaceNodeVelYNm1(unsigned short iMarker, unsigned short iVertex){
    unsigned long iPoint;
    double *Vel_n;

    iPoint = geometry->vertex[iMarker][iVertex];
    Vel_n = geometry->node[iPoint]->GetVel_n();

    return Vel_n[1];
}

double NativeSolidSolver::getInterfaceNodeVelZNm1(unsigned short iMarker, unsigned short iVertex){
    unsigned long iPoint;
    double *Vel_n;

    iPoint = geometry->vertex[iMarker][iVertex];
    Vel_n = geometry->node[iPoint]->GetVel_n();

    return Vel_n[2];
}

double NativeSolidSolver::getRotationCenterPosX(){

    return structure->GetCenterOfRotation_x();

}

double NativeSolidSolver::getRotationCenterPosY(){

    return structure->GetCenterOfRotation_y();
}

double NativeSolidSolver::getRotationCenterPosZ(){

    return structure->GetCenterOfRotation_z();
}

void NativeSolidSolver::setGeneralisedForce(){

  unsigned short iVertex, iMarker;
  unsigned long iPoint;
  double ForceX(0.0), ForceY(0.0), ForceZ(0.0);
  double* force;

  iMarker = getFSIMarkerID();

  for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++){
      iPoint = geometry->vertex[iMarker][iVertex];
      force = geometry->node[iPoint]->GetForce();
      ForceX += force[0];
      ForceY += force[1];
      ForceZ += force[2];
  }

  if(config->GetStructType() == "SPRING_HOR"){
    ((solver->GetLoads())->GetVec())[0] = ForceX;
  }
  else if(config->GetStructType() == "SPRING_VER"){
    ((solver->GetLoads())->GetVec())[0] = ForceY;
  }
  else if(config->GetStructType() == "AIRFOIL"){
    ((solver->GetLoads())->GetVec())[0] = -ForceY;
  }
  else{
    cerr << "Wrong structural type for applying global fluild loads !" << endl;
    throw(-1);
  }

}

void NativeSolidSolver::setGeneralisedMoment(){

  unsigned short iVertex, iMarker;
  unsigned long iPoint;
  double Moment(0.0), CenterX, CenterY, CenterZ;
  double* Force;
  double* Coord;

  iMarker = getFSIMarkerID();

  for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++){
      iPoint = geometry->vertex[iMarker][iVertex];
      Force = geometry->node[iPoint]->GetForce();
      Coord = geometry->node[iPoint]->GetCoord();
      CenterX = getRotationCenterPosX();
      CenterY = getRotationCenterPosY();
      Moment += (Force[1]*(Coord[0]-CenterX) - Force[0]*(Coord[1]-CenterY));
  }

  if(config->GetStructType() == "AIRFOIL"){
    ((solver->GetLoads())->GetVec())[1] = -Moment;
  }
  else if(config->GetStructType() == "SPRING_VER"){}
  else if(config->GetStructType() == "SPRING_HOR"){}
  else{
    cerr << "Wrong structural type for applying global fluild loads !" << endl;
    throw(-1);
  }

}

void NativeSolidSolver::applyload(unsigned short iVertex, double Fx, double Fy, double Fz){

    unsigned short iMarker;
    unsigned long iPoint;
    double Force[3];

    Force[0] = Fx;
    Force[1] = Fy;
    Force[2] = Fz;

    iMarker = getFSIMarkerID();
    iPoint = geometry->vertex[iMarker][iVertex];
    geometry->node[iPoint]->SetForce(Force);

}
