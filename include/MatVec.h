#pragma once

#include <cblas.h>
#include <lapacke.h>

#define ROW_MAJ 0
#define COL_MAJ 1

class CVector {

protected:
    unsigned long nElm;
    double* vec_val;

public:
    CVector(void);
    CVector(const unsigned long & size, const double & val = 0.0);
    CVector(const CVector & u);
    explicit CVector(const unsigned long & size, const double* u_array);
    virtual ~CVector();
    unsigned long GetSize() const;
    double* GetVec() const;
    void print() const;
    void Initialize(const double & val);
    CVector & operator=(const CVector & u);
    CVector & operator+(const CVector & u) const;
    CVector & operator+=(const CVector & u);
    CVector & operator-(const CVector & u) const;
    CVector & operator-=(const CVector & u);
    CVector & operator*(const double & val) const;
    //CVector & operator*(const double & val, const CVector & u);
    CVector & operator*=(const double & val);
    CVector & operator/(const double & val) const;
    CVector & operator/=(const double & val);
    double & operator[](const unsigned long & i) const;
    //const double & operator[](const unsigned long & i) const;
    double norm() const;
    double dotProd(const CVector & u, const CVector & v) const;
    void Reset();
};

class CMatrix {

protected:
    unsigned long nEq;
    unsigned long nVar;
    CVector* mat_val;

public:
    CMatrix(void);
    CMatrix(const unsigned long & val_nEq, const unsigned long & val_nVar, const double & val = 0.0);
    explicit CMatrix(const unsigned long & val_nEq, const unsigned long & val_nVar, const double* a_array);
    virtual ~CMatrix();
    void print() const;
    //void printt() const;
    CMatrix & operator=(const CMatrix & a);
    CMatrix & operator+(const CMatrix & a) const;
    CMatrix & operator+=(const CMatrix & a);
    unsigned long GetnEq() const;
    unsigned long GetnVar() const;
    double* GetMat() const;
    CVector* GetCVec() const;
    void SetElm(int i, int j, double val);
    double GetElm(int i, int j) const;
    double DiagProduct() const;
    double ComputeDet() const;
    void Reset();
};

CVector MatVecProd(CMatrix* A, CVector* b);
CVector ScalVecProd(const double & scal, CVector* b);
CMatrix ScalMatProd(const double & scal, CMatrix* A);
int SolveSys(CMatrix *A, CVector *b);

void MatrixToVec(int order, double** matrix, double* vecteur, int Nrow, int Ncol, int sizeVec);
void VecToMatrix(int order, double** matrix, double* vecteur, int Nrow, int Ncol, int sizeVec);
