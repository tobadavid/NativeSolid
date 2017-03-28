#pragma once

#include "MatVec.h"
#include "config.h"


class Structure{

protected:
    unsigned int nDof;
    double m;
    double Kh;
    double Ka;
    double Ch;
    double Ca;
    double c;
    double b;
    double xf;
    double xCG;
    double S;
    double ICG;
    double If;
    CMatrix M;
    CMatrix C;
    CMatrix K;
    double centerOfRotation[3];
    double centerOfRotation_n[3];
    //double* interfaceCoord;
    //CVector* InitialSolution;


public:
    Structure(Config* config);
    virtual ~Structure();
    //void readMeshSU2;
    void SetStructuralMatrices(Config* config);
    //void SetSolutionInTime();
    CMatrix & GetM();
    CMatrix & GetC();
    CMatrix & GetK();
    void SetCenterOfRotation_X(double coord_x);
    void SetCenterOfRotation_Y(double coord_y);
    void SetCenterOfRotation_Z(double coord_z);
    double GetCenterOfRotation_x() const;
    double GetCenterOfRotation_y() const;
    double GetCenterOfRotation_z() const;
    const double* GetCenterOfRotation() const;
    void SetCenterOfRotation_n_X(double coord_x);
    void SetCenterOfRotation_n_Y(double coord_y);
    void SetCenterOfRotation_n_Z(double coord_z);
    double GetCenterOfRotation_n_x() const;
    double GetCenterOfRotation_n_y() const;
    double GetCenterOfRotation_n_z() const;
    unsigned int GetnDof();

};
