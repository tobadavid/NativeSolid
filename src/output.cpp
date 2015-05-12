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
