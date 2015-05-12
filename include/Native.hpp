#pragma once

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

#include "config.hpp"
#include "structure.hpp"
#include "integration.hpp"
#include "output.hpp"
#include "MatVec.hpp"

//#define PI 3.14159265

const unsigned int BUFSIZE = 3000000;
const int MASTER_NODE = 0;
const int SINGLE_NODE = 1;


int LaunchNative(string FileName);

/*Function y = f(x) class prototype*/
/*class Function
{
public:
    Function();
    virtual ~Function();
    virtual Function* GetAddress();
    virtual void addX(double element);
    virtual void addY(double element);
    virtual void show();
    virtual double getX(int index);
    virtual double getY(int index);
    virtual int GetSizeX();
    virtual int GetSizeY();

protected:
    std::vector<double> x;
    std::vector<double> y;
};

/*Spring class prototype*/
/*class Spring
{
public:
    Spring(CCSDSolver *Structure);
    virtual ~Spring();
    virtual Spring* GetAddress();
    virtual void SetInitialCond(CCSDSolver *Structure);
    virtual void GetForce(CCSDSolver *Structure);
    virtual void ShowForce();
    virtual void Integrate(CCSDSolver *Structure);
    virtual void GetSpringParam(CCSDSolver *Structure);

protected:
    double M;
    double K;
    double C;
    double Q[4][2];
    Function TimeForce;
};
*/
