#pragma once

#include "Config.h"
#include "Structure.h"
#include "Integration.h"
#include "MatVec.h"
#include "Geometry.h"
#include "Output.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

class NativeSolidSolver
{
    std::string confFile;
    Config *config;
    Geometry *geometry;
    Structure *structure;
    Integration *integrator;
    Output *output;
    std::ofstream historyFile;
    std::ofstream historyFile2;
    std::ofstream restartFile;
    CVector q_uM1; // The displacement at the previous FSI iteration
    double omega;
    unsigned long nSolidInterfaceVertex;
    double varCoordNorm;

public:
    NativeSolidSolver(std::string str, bool FSIComp);
    ~NativeSolidSolver();
    void exit();
    // double getVarCoordNorm() const;
    void preprocessIteration(unsigned long ExtIter);
    void timeIteration(double t0, double tf);
    // void mapRigidBodyMotion(bool predicition, bool initialize);
    void computeInterfacePosVel(bool initialize);
    void setInitialDisplacements();
    void staticComputation();
    void writeSolution(double currentTime, double lastTime, double currentFSIIter, unsigned long ExtIter, unsigned long NbExtIter);
    void writeSolution(double time, int FSIter);
    void saveSolution();
    void updateSolution();
    // void updateGeometry();
    // void displacementPredictor();
    unsigned short getFSIMarkerID();
    unsigned long getNumberOfSolidInterfaceNodes(unsigned short iMarker);
    unsigned int getInterfaceNodeGlobalIndex(unsigned short iMarker, unsigned short iVertex);
    double getInterfaceNodePosX(unsigned short iMarker, unsigned short iVertex);
    double getInterfaceNodePosY(unsigned short iMarker, unsigned short iVertex);
    double getInterfaceNodePosZ(unsigned short iMarker, unsigned short iVertex);
    double getInterfaceNodePosX0(unsigned short iMarker, unsigned short iVertex);
    double getInterfaceNodePosY0(unsigned short iMarker, unsigned short iVertex);
    double getInterfaceNodePosZ0(unsigned short iMarker, unsigned short iVertex);
    double getInterfaceNodeDispX(unsigned short iMarker, unsigned short iVertex);
    double getInterfaceNodeDispY(unsigned short iMarker, unsigned short iVertex);
    double getInterfaceNodeDispZ(unsigned short iMarker, unsigned short iVertex);
    double getInterfaceNodeVelX(unsigned short iMarker, unsigned short iVertex);
    double getInterfaceNodeVelY(unsigned short iMarker, unsigned short iVertex);
    double getInterfaceNodeVelZ(unsigned short iMarker, unsigned short iVertex);
    double getInterfaceNodeVelXNm1(unsigned short iMarker, unsigned short iVertex);
    double getInterfaceNodeVelYNm1(unsigned short iMarker, unsigned short iVertex);
    double getInterfaceNodeVelZNm1(unsigned short iMarker, unsigned short iVertex);
    double getRotationCenterPosX();
    double getRotationCenterPosY();
    double getRotationCenterPosZ();
    void setGeneralisedForce();
    void setGeneralisedForce(double Fx, double Fy);
    void setGeneralisedMoment();
    void setGeneralisedMoment(double M);
    void applyload(unsigned short iVertex, double Fx, double Fy, double Fz);
};
