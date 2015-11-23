#pragma once

#include <string>
#include "../include/config.h"
#include "../include/structure.h"
#include "../include/integration.h"
#include "../include/output.h"
#include "../include/MatVec.h"
#include "../include/geometry.h"
#include <fstream>
#include <iostream>
#include <sstream>


class NativeSolidSolver{

protected:
    std::string confFile;
    Config* config;
    Geometry* geometry;
    Structure* structure;
    Integration* solver;
    Output* output;
    std::ofstream historyFile;
    std::ofstream restartFile;
    CVector* q_uM1; //The displacement at the previous FSI iteration
    double omega;
    double* globalFluidLoads; //Used for communications with FluidSolver !!
    double** solidInterface;
    double* solidInterfaceBuffer;
    unsigned long nSolidInterfaceVertex;

public:
    NativeSolidSolver(std::string str);
    ~NativeSolidSolver();
    void initialize(bool FSIComp);
    void exit();
    void inputFluidLoads(double currentTime, double FSI_Load);
    double* getGlobalFluidLoadsArray() const;
    double** getSolidInterface() const;
    const double* getCenterCoordinate() const;
    unsigned long getnSolidInterfaceVertex() const;
    void applyGlobalFluidLoads();
    void timeIteration(double currentTime);
    void mapRigidBodyMotion(bool predicition, bool initialize);
    void staticComputation();
    void writeSolution(double currentTime, double currentFSIIter, unsigned long ExtIter, unsigned long NbExtIter);
    void updateSolution();
    void updateGeometry();
    void outputDisplacements(double* interfRigidDispArray, bool initialize);
    //void outputSolidInterface(double* buffer);
    void displacementPredictor_Old(double* interfRigidDispArray);
    void displacementPredictor();
    void setAitkenCoefficient(unsigned long FSIIter);
    double aitkenRelaxation();

};
