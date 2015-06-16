#pragma once

#include <string>
#include "../include/config.h"
#include "../include/structure.h"
#include "../include/integration.h"
#include "../include/output.h"
#include <fstream>
#include <iostream>
#include <sstream>


class NativeSolidSolver{

protected:
    std::string confFile;
    Config* config;
    Structure* structure;
    Integration* solver;
    Output* output;
    std::ofstream outputFile;

public:
    NativeSolidSolver(std::string str);
    ~NativeSolidSolver();
    void initialize(bool FSIComp);
    void exit();
    void inputFluidLoads(double currentTime, double FSI_Load);
    void timeIteration(double currentTime);
    void writeSolution(double currentTime);
    void updateSolution();
    double outputDisplacements();
    double displacementPredictor();
};
