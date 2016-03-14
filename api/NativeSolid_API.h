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
    unsigned long nSolidInterfaceVertex;
    double varCoordNorm;


public:
    /* NEW generation */
    NativeSolidSolver(std::string str, bool FSIComp);
    ~NativeSolidSolver();
    void exit();
    double getVarCoordNorm() const;
    void timeIteration(double currentTime);
    void mapRigidBodyMotion(bool predicition, bool initialize);
    void setInitialDisplacements();
    void staticComputation();
    void writeSolution(double currentTime, double currentFSIIter, unsigned long ExtIter, unsigned long NbExtIter);
    void updateSolution();
    void updateGeometry();
    void displacementPredictor();
    unsigned short getFSIMarkerID();
    unsigned long getNumberOfSolidInterfaceNodes(unsigned short iMarker);
    unsigned int getInterfaceNodeGlobalIndex(unsigned short iMarker, unsigned short iVertex);
    double getInterfaceNodePosX(unsigned short iMarker, unsigned short iVertex);
    double getInterfaceNodePosY(unsigned short iMarker, unsigned short iVertex);
    double getInterfaceNodePosZ(unsigned short iMarker, unsigned short iVertex);
    double getRotationCenterPosX();
    double getRotationCenterPosY();
    double getRotationCenterPosZ();
    void setGeneralisedForce(double ForceX , double ForceY);
    void setGeneralisedMoment(double Moment);
};
