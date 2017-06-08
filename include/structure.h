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
    double centerOfRotation[3];
    double centerOfRotation_n[3];


public:
    Structure(Config* config);
    virtual ~Structure();
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
    unsigned int GetnDof() const;
    double Get_m() const;
    double Get_Kh() const;
    double Get_Ka() const;
    double Get_Ch() const;
    double Get_Ca() const;
    double Get_S() const;
    double Get_If() const;

};
