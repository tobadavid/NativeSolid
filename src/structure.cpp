#include "../include/structure.h"
#include <iostream>
#include <cstdlib>
#include <stdlib.h>
#include <fstream>
#include <sstream>


using namespace std;

Structure::Structure(Config* config){
  if(config->GetStructType() == "SPRING_HOR" || config->GetStructType() == "SPRING_VER" ){
    nDof = 1;
    m = config->GetSpringMass();
    Kh = config->GetSpringStiffness();
    Ch = config->GetSpringDamping();
    Ka = 0;
    Ca = 0;
    xf = 0;
    xCG = 0;
    c = 0;
    b = 0;
    S = 0;
    ICG = 0;
    If = 0;
    cout << "Setting mass-spring-damper system" << endl;
    cout << "Number of DOF : " << nDof << endl;
  }
  else if (config->GetStructType() == "AIRFOIL"){
    nDof = 2;
    m = config->GetSpringMass();
    Kh = config->GetSpringStiffness();
    Ch = config->GetSpringDamping();
    Ka = config->GetTorsionalStiffness();
    Ca = config->GetTorsionalDamping();
    xf = config->GetFlexuralAxis();
    xCG = config->GetGravityCenter();
    c = config->GetCord();
    b = c/2.0;
    S = m*(xCG-xf);
    //If = ICG + m*pow((xCG-xf),2);
    If = config->GetInertiaFlexural();
    //S = m*(c/2.0-xf);
    //Ia = 1.0/3.0*m*(c*c-3*c*xf+3*xf*xf);
    cout << "Setting pitching-plunging airfoil system" << endl;
    cout << "Number of DOF : " << nDof << endl;
  }
  else nDof = 0;

  M = new CMatrix(nDof,nDof);
  C = new CMatrix(nDof, nDof);
  K = new CMatrix(nDof, nDof);
}

Structure::~Structure(){
  delete M;
  delete C;
  delete K;
}

void Structure::SetStructuralMatrices(Config* config){
  if(config->GetStructType() == "SPRING_HOR" || config->GetStructType() == "SPRING_VER" ){
    K->SetElm(1,1,Kh);
    C->SetElm(1,1,Ch);
    M->SetElm(1,1,m);

    cout << "Plunging mass : " << m << " [kg]" << endl;
    cout << "Plunging damping : " << Ch << " [Ns/m]" << endl;
    cout << "Plunging stiffness : " << Kh << " [N/m]" << endl;
  }
  else if (config->GetStructType() == "AIRFOIL"){
    M->SetElm(1,1,m);
    M->SetElm(1,2,S);
    M->SetElm(2,1,S);
    M->SetElm(2,2,If);
    K->SetElm(1,1,Kh);
    K->SetElm(2,2,Ka);
    C->SetElm(1,1,Ch);
    C->SetElm(2,2,Ca);

    cout << "Airfoil mass : " << m << " [kg]" << endl;
    cout << "Airfoil cord : " << c << " [m]" << endl;
    cout << "Position of the flexural axis : " << xf << " [m]" << endl;
    cout << "Inertia around the flexural axis : " << If << " [kg m²]" << endl;
    cout << "Plunging damping : " << Ch << " [Ns/m]" << endl;
    cout << "Plunging stiffness : " << Kh << " [N/m]" << endl;
    cout << "Pitching damping : " << Ca << " [Ns]" << endl;
    cout << "Pitching stiffness : " << Ka << " [N]" << endl;
    cout << "Position of the center of gravity : " << xCG << " [m]" << endl;
    cout << "Static unbalance : " << S << " [kg m]" << endl;
  }
  else{
    cerr << "Invalid structural type. Available choices are : SPRIN_HOR, SPRING_VER and AIRFOIL." << endl;
    throw(-1);
  }
}

CMatrix* Structure::GetM(){
  return M;
}

CMatrix* Structure::GetC(){
  return C;
}

CMatrix* Structure::GetK(){
  return K;
}

unsigned int Structure::GetnDof(){
  return nDof;
}
