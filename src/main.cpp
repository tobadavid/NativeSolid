#include <iostream>
#include <stdlib.h>

#include "../include/Native.h"

using namespace std;

int main(int argc, char* argv[]){

  if (argc <= 1){
    cerr << "No configuration file name !" << endl;
    throw(-1);
  }

  LaunchNative(argv[1]);


  return EXIT_SUCCESS;
}
