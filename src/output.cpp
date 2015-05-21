#include "../include/output.hpp"

Output::Output(){}

Output::~Output(){}

void Output::WriteHistory(Integration* solver, Structure* structure, ofstream* outputfile, const double & time){

  if(structure->GetnDof() == 1){
    if(time == 0){
      cout << "\"Time\"" << "\t" << "\"Displacement\"" << "\t" << "\"Velocity\"" << "\t" << "\"Acceleration\"" << "\t" << endl;
      outputfile[0] << "\"Time\"" << "\t" << "\"Displacement\"" << "\t" << "\"Velocity\"" << "\t" << "\"Acceleration\"" << "\t" << endl;
    }
    cout << time << "\t" << (*(solver->GetDisp()))[0] << "\t" << (*(solver->GetVel()))[0] << "\t" << (*(solver->GetAcc()))[0] << endl;
    outputfile[0] << time << "\t" << (*(solver->GetDisp()))[0] << "\t" << (*(solver->GetVel()))[0] << "\t" << (*(solver->GetAcc()))[0] << endl;
  }
  else if(structure->GetnDof() == 2){
    if(time == 0){
      cout << "\"Time\"" << "\t" << "\"Displacement 1\"" << "\t" << "\"Displacement 2\"" << "\t" << "\"Velocity 1\""  << "\t" << "\"Velocity 2\"" << "\t" << "\"Acceleration 1\"" << "\t" << "\"Acceleration 2\"" << endl;
      outputfile[0] << "\"Time\"" << "\t" << "\"Displacement 1\"" << "\t" << "\"Displacement 2\"" << "\t" << "\"Velocity 1\""  << "\t" << "\"Velocity 2\"" << "\t" << "\"Acceleration 1\"" << "\t" << "\"Acceleration 2\"" << endl;
    }
    cout << time << "\t" << (*(solver->GetDisp()))[0] << "\t" << (*(solver->GetDisp()))[1] << "\t" << (*(solver->GetVel()))[0] << "\t" << (*(solver->GetVel()))[1] << "\t" << (*(solver->GetAcc()))[0] << "\t" << (*(solver->GetAcc()))[1] << endl;
    outputfile[0] << time << "\t" << (*(solver->GetDisp()))[0] << "\t" << (*(solver->GetDisp()))[1] << "\t" << (*(solver->GetVel()))[0] << "\t" << (*(solver->GetVel()))[1] << "\t" << (*(solver->GetAcc()))[0] << "\t" << (*(solver->GetAcc()))[1] << endl;
  }
}

void Output::WriteRestart(Integration* solver, Structure* structure){
  ofstream RestartFile;
  RestartFile.open("Nat_solution_restart.out", ios::out);

  if(structure->GetnDof() == 1){
    RestartFile << "\"Displacement\"" << "\t" << "\"Velocity\"" << "\t" << "\"Acceleration\"" << "\t" << "\"Acceleration variable\"" << endl;
    RestartFile << (*(solver->GetDisp()))[0] << "\t" << (*(solver->GetVel()))[0] << "\t" << (*(solver->GetAcc()))[0] << "\t" << (*(solver->GetAccVar()))[0] << endl;
  }
  else if(structure->GetnDof() == 2){
    RestartFile << "\"Displacement 1\"" << "\t" << "\"Displacement 2\"" << "\t" << "\"Velocity 1\""  << "\t" << "\"Velocity 2\"" << "\t" << "\"Acceleration 1\"" << "\t" << "\"Acceleration 2\"" << "\t" << "\"Acceleration variable 1\"" << "\t" << "\"Acceleration variable 2\"" << endl;
    RestartFile << time << "\t" << (*(solver->GetDisp()))[0] << "\t" << (*(solver->GetDisp()))[1] << "\t" << (*(solver->GetVel()))[0] << "\t" << (*(solver->GetVel()))[1] << "\t" << (*(solver->GetAcc()))[0] << "\t" << (*(solver->GetAcc()))[1]<< "\t" << (*(solver->GetAccVar()))[0] << "\t" << (*(solver->GetAccVar()))[1] << endl;
  }
}

void Output::WriteStaticSolution(Config* config, Integration* solver, Structure* structure, ofstream* outputfile){
  if(structure->GetnDof() == 1){
    cout << "Static displacement is : " << (*(solver->GetDisp()))[0] << " [m]" << endl;
    cout << "Writing displacement into a solution file" << endl;
    if(config->GetStructType() == "SPRING_HOR"){
      outputfile[0] << -1.000 << "\t" << (*(solver->GetDisp()))[0] << "\t" << 0.000;
    }
    else if(config->GetStructType() == "SPRING_VER"){
      outputfile[0] << -1.000 << "\t" << 0.000 << "\t" << (*(solver->GetDisp()))[0];
    }
  }
  else if(structure->GetnDof() == 2){
    cout << "Plunging displacement is :" << (*(solver->GetDisp()))[0] << " [m]" << endl;
    cout << "Pitching rotation is :" << (*(solver->GetDisp()))[1] << " [rad]" << endl;
    cout << "Writing displacement into a solution file" << endl;
    outputfile[0] << -1.000 << "\t" << (*(solver->GetDisp()))[0] << "\t" << (*(solver->GetDisp()))[1];
  }
}
