#pragma once

#include <string>
#include "../include/config.h"
#include "../include/structure.h"
#include "../include/integration.h"
#include "../include/output.h"
#include "../include/MatVec.h"
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
    CVector* q_uM1; //The displacement at the previous FSI iteration
    double omega;
    double* globalFluidLoads;

public:
    NativeSolidSolver(std::string str);
    ~NativeSolidSolver();
    void initialize(bool FSIComp);
    void exit();
    void inputFluidLoads(double currentTime, double FSI_Load);
    double* getGlobalFluidLoadsArray() const;
    void applyGlobalFluidLoads();
    void timeIteration(double currentTime);
    void staticComputation();
    void writeSolution(double currentTime, double currentFSIIter);
    void updateSolution();
    void outputDisplacements(double* interfRigidDispArray);
    void displacementPredictor(double* interfRigidDispArray);
    void setAitkenCoefficient(unsigned long FSIIter);
    double aitkenRelaxation();

};
