#pragma once

#include <fstream>

#include "MatVec.h"
#include "config.h"
#include "structure.h"
#include "integration.h"

class Output
{
public:
    Output();
    ~Output();
    // void WriteHistory(Integration* solver, Structure* structure, std::ofstream* outputfile, const double & time);
    // void WriteRestart(Integration* solver, Structure* structure);
    // void WriteStaticSolution(Config* config, Integration* solver, Structure* structure, std::ofstream* outputfile);
    // void WriteRestart();
};
