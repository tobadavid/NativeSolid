#include "../include/MatVec.h"
#include <iostream>
#include <cmath>

using namespace std;

/*--- Definition of the CVector class ---*/

CVector::CVector(void){

  vec_val = NULL;

}

CVector::CVector(const unsigned long & size, const double & val){

  nElm = size;

  if(nElm <= 0){
    cerr << "CVector:CVector(const unsigned long &, const double): " << "invalid input : size = " << size << endl;
    throw(-1);
  }

  vec_val = new double[nElm];
  for(unsigned int i=0; i < nElm; i++)
    vec_val[i] = val;

}

CVector::CVector(const CVector & u){

  nElm = u.nElm;

  vec_val = new double[nElm];
  for(unsigned int i = 0; i < nElm; i++)
    vec_val[i] = u.vec_val[i];
}

CVector::CVector(const unsigned long & size, const double* u_array){

  nElm = size;

  if(nElm <= 0){
    cerr << "CVector::CVector(const unsigned long &, const double): " << "invalid input: size = " << size << endl;
    throw(-1);
  }

  vec_val = new double[nElm];
  for(unsigned int i=0; i<nElm; i++)
    vec_val[i] = u_array[i];
}

CVector::~CVector(){
  delete [] vec_val;
  vec_val = NULL;
  nElm = -1;
}


void CVector::Initialize(const double & val){
  if(nElm <= 0){
    cerr << "CVector::Initialize(const double &): " << "invalid number of element" << nElm << endl;
    throw(-1);
  }

  for(unsigned int i= 0; i < nElm; i++)
    vec_val[i] = val;
}

unsigned long CVector::GetSize() const{
  if(nElm <= 0){
    cerr << "CVector::GetSize() const: " << "invalid number of element" << nElm << endl;
    throw(-1);
  }

  return nElm;
}

double* CVector::GetVec() const{
  return vec_val;
}

void CVector::print() const{
  for(unsigned int i=0; i<nElm; i++)
    cout << vec_val[i] << endl;
}

CVector & CVector::operator=(const CVector & u){
  if(this == &u) return *this;

  if(nElm != u.nElm){
    delete [] vec_val;
    nElm = u.nElm;
    vec_val = new double[nElm];
  }

  for(unsigned int i=0; i<nElm; i++)
    vec_val[i] = u.vec_val[i];

  return *this;
}

CVector & CVector::operator+(const CVector & u) const{
   if(nElm != u.nElm){
    cerr << "CVector::operator+(const CVector &) const: " << "sizes do not match" << endl;
    throw(-1);
   }

   CVector sum(*this);
   sum += u;
   return sum;
 }

CVector & CVector::operator+=(const CVector & u){
   if(nElm != u.nElm){
    cerr << "CVector::operator+=(const CVector &) const: " << "sizes do not match" << endl;
    throw(-1);
   }

    for(unsigned int i=0; i < nElm; i++)
        vec_val[i] += u.vec_val[i];

    return *this;
 }

CVector & CVector::operator-(const CVector & u) const{
   if(nElm != u.nElm){
    cerr << "CVector::operator-(const CVector &) const: " << "sizes do not match" << endl;
    throw(-1);
   }

   CVector diff(*this);
   diff -= u;

   return diff;
 }

 CVector & CVector::operator-=(const CVector & u){
   if(nElm != u.nElm){
    cerr << "CVector::operator-=const CVector &) const: " << "sizes do not match" << endl;
    throw(-1);
   }

   for(unsigned int i=0; i<nElm; i++)
    vec_val[i] -= u.vec_val[i];

   return *this;
 }

CVector & CVector::operator*(const double & val) const{
  CVector prod(*this);
  prod *= val;
  return prod;
}

/*CVector & CVector::operator*(const double & val, const CVector & u){
  CVector prod(u);
  prod *= val;
  return prod;
}*/

CVector & CVector::operator*=(const double & val){
  for(unsigned int i=0; i<nElm; i++)
    vec_val[i] *= val;

  return *this;
}

CVector & CVector::operator/(const double & val) const{
  CVector quotient(*this);
  quotient /= val;
  return quotient;
}

CVector & CVector::operator/=(const double & val){
  for(unsigned int i=0; i<nElm; i++)
    vec_val[i] /= val;

  return *this;
}

double & CVector::operator[](const unsigned long & i) const{
  return vec_val[i];
}

double CVector::dotProd(const CVector & u, const CVector & v) const{
  if(u.nElm != v.nElm){
    cerr << "CVector::dotProd(const CVector &, const CVector &): " << "sizes do not math" << endl;
    throw(-1);
  }

  double dotProd(0);
  for(unsigned int i=0; i<u.nElm; i++)
    dotProd += u.vec_val[i]*v.vec_val[i];

  return dotProd;
}

double CVector::norm() const{
  double norm(0), tempVal(0);
  tempVal = dotProd(*this,*this);

  if(tempVal < 0){
    cerr << "CVector::norm(): " << "computed dot product is negative: " << tempVal << endl;
    throw(-1);
  }

  norm = sqrt(tempVal);

  return norm;
}

void CVector::Reset(){
  for(unsigned int i=0; i<nElm; i++)
    vec_val[i] = 0;
}

/*--- Definition of the CMatrix class ---*/

CMatrix::CMatrix(void){
  mat_val = NULL;
}

CMatrix::CMatrix(const unsigned long & val_nEq, const unsigned long & val_nVar, const double & val){
  nEq = val_nEq;
  nVar = val_nVar;
  unsigned long nElm = nEq*nVar;

  mat_val = new CVector(nElm,val);
}

CMatrix::CMatrix(const unsigned long & val_nEq, const unsigned long & val_nVar, const double* a_array){
  nEq = val_nEq;
  nVar = val_nVar;
  const unsigned long nElm = nEq*nVar;
  mat_val = new CVector(nElm, a_array);
}

CMatrix::~CMatrix(){
  delete mat_val;
  mat_val = NULL;
}

void CMatrix::print() const{
  mat_val->print();
}

CMatrix & CMatrix::operator=(const CMatrix & a){
  if(this == &a) return *this;

  if(nEq != a.nEq || nVar != a.nVar){
    delete mat_val;
    nEq = a.nEq;
    nVar = a.nVar;
    mat_val = new CVector(nEq*nVar);
  }

  *mat_val = *(a.mat_val);


  return *this;
}

CMatrix & CMatrix::operator+(const CMatrix & a) const{
  if(nEq != a.nEq || nVar != a.nVar){
    cerr << "CMatrix::operator+(const CMatrix &) const: " << "sizes do not match" << endl;
    throw(-1);
  }

  CMatrix sum(*this);
  sum += a;
  return sum;
}

CMatrix & CMatrix::operator+=(const CMatrix & a){
  if(nEq != a.nEq || nVar != a.nVar){
    cerr << "CMatrix::operator+(const CMatrix &) const: " << "sizes do not match" << endl;
    throw(-1);
  }

  *mat_val += *(a.mat_val);
}

unsigned long CMatrix::GetnEq() const{
  return nEq;
}

unsigned long CMatrix::GetnVar() const{
  return nVar;
}

double* CMatrix::GetMat() const{
  return mat_val->GetVec();
}

CVector* CMatrix::GetCVec() const{
  return mat_val;
}

void CMatrix::SetElm(int i, int j, double val){
  (*mat_val)[(j-1)*nEq+(i-1)] = val;
}

double CMatrix::GetElm(int i, int j) const{
  return (*mat_val)[(j-1)*nEq+(i-1)];
}

double CMatrix::DiagProduct() const{
  double prod(1);
  if (nEq == nVar){
    for(int i=1; i<=nEq; i++){
      prod *= GetElm(i,i);
    }
  }

  return prod;
}

double CMatrix::ComputeDet() const{
  double Det(1);
  int NbChange(0);
  double PLU[nEq*nVar];
  int IPIV[nEq];
  int INFO(-1);

  *PLU = *(mat_val->GetVec());
  INFO = LAPACKE_dgetrf(LAPACK_COL_MAJOR,nEq,nVar,PLU,nEq,IPIV);

  for(int i=0; i<nEq-1; i++){
    if(IPIV[i] != i+1) NbChange++;
  }
}

void CMatrix::Reset(){
  mat_val->Reset();
}

CVector MatVecProd(CMatrix* A, CVector* b){
  CVector x(A->GetnVar());
  cblas_dgemv(CblasColMajor,CblasNoTrans,A->GetnEq(),A->GetnVar(),1,A->GetMat(),A->GetnEq(),b->GetVec(),1,0,x.GetVec(),1);

   return x;
}

CVector ScalVecProd(const double & scal, CVector* b){
  CVector x(b->GetSize());
  cblas_daxpy(b->GetSize(),scal,b->GetVec(),1,x.GetVec(),1);
  return x;
}

CMatrix ScalMatProd(const double & scal, CMatrix* A){
  CMatrix U(A->GetnEq(),A->GetnVar(), A->GetMat());
  CVector* temp = U.GetCVec();
  *temp *= scal;

  return U;
}

int SolveSys(CMatrix *A, CVector *b){
  const int N = A->GetnEq();
  const int LDA = A->GetnVar();
  CMatrix* Mat;
  Mat = new CMatrix(N,LDA,A->GetMat());
  int NRHS(1);
  int *IPIV;
  int LDB = b->GetSize();
  int INFO(-1);

  IPIV = new int[N];

  INFO = LAPACKE_dgesv(LAPACK_COL_MAJOR,N,NRHS,Mat->GetMat(),LDA,IPIV,b->GetVec(),LDB);

  delete Mat;
  Mat = NULL;
}

void MatrixToVec(int order, double** matrix, double* vecteur, int Nrow, int Ncol, int sizeVec){
    if(order == COL_MAJ){
      for(int j=0; j<sizeVec; j++)
        vecteur[j] = matrix[j%Nrow][j/Nrow];
    }
    else if(order == ROW_MAJ){
      for(int j=0; j<sizeVec; j++)
        vecteur[j] = matrix[j/Ncol][j%Ncol];
    }
    else{
      cerr << "Wrong storage order" << endl;
      throw(-1);
    }
}

void VecToMatrix(int order, double** matrix, double* vecteur, int Nrow, int Ncol, int sizeVec){
    if(order == COL_MAJ){
      for(int i=0; i<Nrow;i++)
        for(int j=0; j<Ncol; j++)
          matrix[i][j] = vecteur[j*Nrow+i];
    }
    else if(order == ROW_MAJ){
      for(int i=0; i<Nrow;i++)
        for(int j=0; j<Ncol; j++)
          matrix[i][j] = vecteur[i*Ncol+j];
    }
    else{
      cerr << "Wrong storage order" << endl;
      throw(-1);
    }
}
