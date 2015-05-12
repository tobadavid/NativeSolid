#ifdef HAVE_MPI
    #include "mpi.h"
#endif

#include <cstring>
#include <cstdlib>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cmath>
#include <algorithm>

#include "../include/Native.hpp"

using namespace std;

int main(int argc, char* argv[]){

  if (argc <= 1){
    cerr << "No configuration file name !" << endl;
    throw(-1);
  }

  LaunchNative(argv[1]);


  return EXIT_SUCCESS;
}
