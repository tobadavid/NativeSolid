/*!
 * \file Native.cpp
 * \brief Main file for the Native structural solver
 * \author THOMAS David, University of Liège, Belgium. Aerospace and Mechanical Engineering Department
 * \version BETA
*/

#include "../include/Native.h"
#include "config.h"
#include "structure.h"
#include "integration.h"
#include "output.h"
#include "MatVec.h"
#include <iostream>

#ifdef HAVE_MPI
  #include "mpi.h"
#endif // HAVE_MPI

using namespace std;

//int main(int argc, char* argv[])
int LaunchNative(string FileName)
{
    /*if (argc <= 1){
      cerr << "No configuration file name !" << endl;
      throw(-1);
    }*/

    int rank = MASTER_NODE;
    int size = SINGLE_NODE;
    double currentTime(0),totTime(0),deltaT(0);

    Config* config;
    Structure* structure;
    Integration* solver;
    Output* output;

    ofstream outputfile;
    outputfile.open("NativeOutput.out", ios::out);

    /*--- MPI initialization, and buffer setting ---*/
#ifdef HAVE_MPI
  int *bptr, bl;
  MPI_Init(&argc,&argv);
  MPI_Buffer_attach( malloc(BUFSIZE), BUFSIZE );
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

    /*--- Initialize the main containers ---*/
    //config = new Config(argv[1]);
    config = new Config(FileName);
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
      totTime = solver->GettotTime();
      deltaT = solver->GetdeltaT();

      cout << endl << "\n----------------------- Begin temporal integration ----------------------" << endl;
      solver->SetLoadsAtTime(config, structure, currentTime);
      solver->SetInitialConditions(config, structure);
      output->WriteHistory(solver, structure, &outputfile, currentTime);
      solver->UpdateSolution();
      currentTime += deltaT;

      while(currentTime <= totTime+deltaT){
        solver->SetLoadsAtTime(config, structure, currentTime);
        solver->TemporalIteration(config, structure);
        output->WriteHistory(solver, structure, &outputfile, currentTime);
        solver->UpdateSolution();
        currentTime += deltaT;
      }
      cout << endl << "\n----------------------- Temporal integration successfully ended ----------------------" << endl;
    }
    else{
      cout << endl << "\n----------------------- Compute static displacement ----------------------" << endl;
      solver->SetStaticLoads(config,structure);
      solver->StaticIteration(config,structure);
      output->WriteStaticSolution(config, solver, structure, &outputfile);
      cout << endl << "\n----------------------- Successfull computation ----------------------" << endl;
    }
    outputfile.close();

    delete config;
    delete structure;
    delete solver;
    delete output;

#ifdef HAVE_MPI
  /*--- Finalize MPI parallelization ---*/
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Buffer_detach(&bptr,&bl);
  MPI_Finalize();
#endif

    return 0;
}
