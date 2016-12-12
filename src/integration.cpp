#include "../include/integration.h"
#include <iostream>
#include <cstdlib>
#include <stdlib.h>
#include <fstream>
#include <cmath>

using namespace std;

Integration::Integration(Structure *structure){

  unsigned int nDof = structure->GetnDof();
  q = new CVector(nDof);
  qdot = new CVector(nDof);
  qddot = new CVector(nDof);
  a = new CVector(nDof);
  q_n = new CVector(nDof);
  qdot_n = new CVector(nDof);
  qddot_n = new CVector(nDof);
  a_n = new CVector(nDof);

  Loads = new CVector(nDof);
}

Integration::~Integration(){
  delete q;
  delete qdot;
  delete qddot;
  delete a;
  delete q_n;
  delete qdot_n;
  delete qddot_n;
  delete a_n;

  delete Loads;
}

double Integration::GettotTime(){
  return totTime;
}

double Integration::GetdeltaT(){
  return deltaT;
}

CVector* Integration::GetDisp() const{
  return q;
}

void Integration::SetDisp(double & val_disp){
  (*q)[0] = val_disp;
}

void Integration::SetExtIter(unsigned long val_ExtIter){

    ExtIter = val_ExtIter;
}

unsigned long Integration::GetExtIter(){

    return ExtIter;
}

CVector* Integration::GetVel() const{
  return qdot;
}
CVector* Integration::GetAcc() const{
  return qddot;
}

CVector* Integration::GetAccVar() const{
  return a;
}

CVector* Integration::GetDisp_n() const{
  return q_n;
}
CVector* Integration::GetVel_n() const{
  return qdot_n;
}

CVector* Integration::GetAcc_n() const{
  return qddot_n;
}

CVector* Integration::GetAccVar_n() const{
  return a_n;
}

CVector* Integration::GetLoads() const{

  return Loads;

}

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

void Integration::SetStaticLoads(Config* config, Structure* structure){
  if(config->GetForceInputType() == "FILE"){
    cout << "Setting applied force from a file : " << config->GetForceInputFile() << endl;
    string delimiter = "\t";
    size_t pos;
    string ForceFileName = config->GetForceInputFile();
    string text_line, tempString, token;
    ifstream InputFile;
    InputFile.open(ForceFileName.c_str(), ios::in);
        while (getline(InputFile,text_line)){
            tempString = text_line;
        }
        InputFile.close();
    pos = tempString.find(delimiter);
    token = tempString.substr(0,pos);
    tempString.erase(0,pos+delimiter.length());
    if (config->GetStructType() == "SPRING_HOR") (*Loads)[0] = atof(tempString.c_str());
    else if (config->GetStructType() == "SPRING_VER") (*Loads)[0] = atof(token.c_str());
    else{
      cerr << "Not ready for AIRFOIL implementation" << endl;
      throw(-1);
    }
  }
  else if(config->GetForceInputType() == "ANALYTICAL" && config->GetAnalyticalFunction() == "CONSTANT"){
    (*Loads)[0] = config->GetConstantForceValue();
    if(structure->GetnDof() == 2) (*Loads)[1] = 0.0;
  }
  else{
    cerr << "For static computation, the analytical function for the applied load must be CONSTANT" << endl;
    throw(-1);
  }
}

void Integration::SetInitialConditions(Config* config, Structure* structure){
  string delimiter = "\t";
  size_t pos;
  string textLine;
  ifstream Inputfile;
  string RestartFile;
  //double* qddot[nDof];
  CVector* RHS = new CVector(structure->GetnDof(),double(0));

  if(config->GetRestartSol() == "YES"){
    string InputFileName = config->GetRestartFile();
    string text_line;
    string token, tempString;
    size_t pos;
    string delimiter = "\t";
    ifstream InputFile;
    InputFile.open(InputFileName.c_str(), ios::in);
    double buffer[(4*structure->GetnDof())+1];
    int kk = 0;
    int jj;
	while (getline(InputFile,text_line)){
      tempString = text_line;
      jj = 0;
      if (kk == 1){
        while ((pos = tempString.find(delimiter)) != string::npos){
          token = tempString.substr(0,pos);
          tempString.erase(0,pos+delimiter.length());
          buffer[jj] = atof(token.c_str());
          jj += 1;
        }
	    buffer[jj] = atof(tempString.c_str());

	    if(structure->GetnDof() == 1){
	      (*q_n)[0] = buffer[1];
	      (*qdot_n)[0] = buffer[2];
          (*qddot_n)[0] = buffer[3];
	      (*a_n)[0] = buffer[4];
	    }
	    else if (structure->GetnDof() == 2){
	      (*q_n)[0] = buffer[1];
	      (*q_n)[1] = buffer[2];
	      (*qdot_n)[0] = buffer[3];
	      (*qdot_n)[1] = buffer[4];
	      (*qddot_n)[0] = buffer[5];
	      (*qddot_n)[1] = buffer[6];
	      (*a_n)[0] = buffer[7];
	      (*a_n)[1] = buffer[8];
	    }
	    q_n->print();
	    qdot_n->print();
	    qddot_n->print();
	    a_n->print();
      }
      else if (kk == 2){
        while ((pos = tempString.find(delimiter)) != string::npos){
          token = tempString.substr(0,pos);
          tempString.erase(0,pos+delimiter.length());
          buffer[jj] = atof(token.c_str());
          jj += 1;
        }
	    buffer[jj] = atof(tempString.c_str());

	    if(structure->GetnDof() == 1){
	      (*q)[0] = buffer[1];
	      (*qdot)[0] = buffer[2];
          (*qddot)[0] = buffer[3];
	      (*a)[0] = buffer[4];
	    }
	    else if (structure->GetnDof() == 2){
	      (*q)[0] = buffer[1];
	      (*q)[1] = buffer[2];
	      (*qdot)[0] = buffer[3];
	      (*qdot)[1] = buffer[4];
	      (*qddot)[0] = buffer[5];
	      (*qddot)[1] = buffer[6];
	      (*a)[0] = buffer[7];
	      (*a)[1] = buffer[8];
	    }
	    q->print();
	    qdot->print();
	    qddot->print();
	    a->print();
      }
      kk += 1;
    }
    InputFile.close();
  }
  else{
    cout << "Setting basic initial conditions" << endl;
    SetLoadsAtTime(config, structure, 0.0, 0.0);
    ExtIter = 0;
    q->Reset();
    q_n->Reset();
    cout << "Read initial configuration" << endl;
    (*q)[0] = config->GetInitialDisp();
    if(structure->GetnDof() == 2) (*q)[1] = config->GetInitialAngle();
    cout << "Initial plunging displacement : " << (*q)[0] << endl;
    cout << "Initial pitching displacement : " << (*q)[1] << endl;

    qdot->Reset();
    qddot->Reset();
    *RHS += *Loads;
    *RHS -= MatVecProd(structure->GetC(),qdot);
    *RHS -= MatVecProd(structure->GetK(),q);
    SolveSys(structure->GetM(),RHS);
    *qddot = *RHS;
    *a = *qddot;
    //cout << (*q)[0] << endl;
    //cout << (*q)[1] << endl;
  }

  structure->SetCenterOfRotation_Y((*q)[0]);

  delete RHS;
  RHS = NULL;
}

void Integration::ComputeResidual(CMatrix* M, CMatrix* C, CMatrix* K, CVector* res){
  *res += MatVecProd(M,qddot);
  *res += MatVecProd(C,qdot);
  *res += MatVecProd(K,q);
  *res -= *Loads;
}

void Integration::ComputeTangentOperator(Structure* structure, CMatrix* St){
  *St += ScalMatProd(betaPrime,structure->GetM());
  *St += ScalMatProd(gammaPrime,structure->GetC());
  *St += *(structure->GetK());
}

void Integration::TemporalIteration(Config *config, Structure *structure){

  double epsilon = 1e-6;

  /*--- Prediction phase ---*/
  qddot->Reset();

  a->Reset();
  *a += ScalVecProd(alpha_f/(1-alpha_m),qddot_n);
  *a -= ScalVecProd(alpha_m/(1-alpha_m),a_n);

  *q = *q_n;
  *q += ScalVecProd(deltaT,qdot_n);
  *q += ScalVecProd((0.5-beta)*deltaT*deltaT,a_n);
  *q += ScalVecProd(deltaT*deltaT*beta,a);

  *qdot = *qdot_n;
  *qdot += ScalVecProd((1-gamma)*deltaT,a_n);
  *qdot += ScalVecProd(deltaT*gamma,a);

  /*--- Tangent operator and corrector computation ---*/
  CVector* res;
  CVector* Deltaq;
  CMatrix* St;
  res = new CVector(qddot->GetSize(),double(0));
  Deltaq = new CVector(qddot->GetSize(),double(0));
  St = new CMatrix(qddot->GetSize(),qddot->GetSize(),double(0));
  ComputeResidual(structure->GetM(),structure->GetC(),structure->GetK(),res);
  while (res->norm() >= epsilon){
    St->Reset();
    ComputeTangentOperator(structure,St);
    SolveSys(St,res);
    //*res -= ScalVecProd(double(2),res); //=deltaq
    *Deltaq += ScalVecProd(-1,res);
    *q += *Deltaq;
    *qdot += ScalVecProd(gammaPrime,Deltaq);
    *qddot += ScalVecProd(betaPrime,Deltaq);
    res->Reset();
    Deltaq->Reset();
    ComputeResidual(structure->GetM(),structure->GetC(),structure->GetK(),res);
  }
  *a += ScalVecProd((1-alpha_f)/(1-alpha_m),qddot);

  delete res;
  delete Deltaq;
  delete St;
  Deltaq = NULL;
  res = NULL;
  St = NULL;
}

void Integration::StaticIteration(Config *config, Structure *structure){
  q->Reset();
  *q += *Loads;
  SolveSys(structure->GetK(),q); //q will be updated with the new solution...
}

void Integration::UpdateSolution(){

  *q_n = *q;
  q->Reset();
  *qdot_n = *qdot;
  qdot->Reset();
  *qddot_n = *qddot;
  qddot->Reset();
  *a_n = *a;
  a->Reset();
}
