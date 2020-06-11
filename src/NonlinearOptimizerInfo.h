#ifndef NONLINEAROPTIMIZATIONINFO_H
#define NONLINEAROPTIMIZATIONINFO_H
#include "RobotInfo.h"
#include "snoptProblem.hpp"

using namespace std;

struct NonlinearOptimizerInfo
{
  // This struct serves as the high-level interface to call SNOPT or IPOPT
  NonlinearOptimizerInfo()
  {
    nS = 0;
    n = 0;          neF = 0;
    ObjRow = 0;     ObjAdd = 0;
    DerivativeFlag = 0;
  }
  void InnerVariableInitialize(const int& _n, const int& _neF)
  { // This function can only be called after n and neF have been given values
    n = _n;
    neF = _neF;
    x      = new double[_n];
    xlow   = new double[_n];
    xupp   = new double[_n];
    xmul   = new double[_n];
    xstate = new    int[_n];

    F      = new double[_neF];
    Flow   = new double[_neF];
    Fupp   = new double[_neF];
    Fmul   = new double[_neF];
    Fstate = new int[_neF];
  }
  void InnerGraidentInitialize(const int &_lenA, const int & _lenG)
  { // This function can only be called after lenA and lenG have been given values
    if (_lenA != 0)
    {
      // There are cases where all the constraints are nonlinear.
      lenA   = _lenA;
      iAfun = new int[lenA];
      jAvar = new int[lenA];
      A  = new double[lenA];
    }
    lenG   = _lenG;
    iGfun = new int[lenG];
    jGvar = new int[lenG];
    DerivativeFlag = 1;
  }
  void JacobianPositionUpdate(const int & _neG, const std::vector<int> &_iGfun, const std::vector<int> &_jGvar, const std::vector<int> &_iAfun, const std::vector<int> &_jAvar, const int & _neA, const std::vector<double> &_A)
  {
    // This function is used to save the position of the position of the Jacobian matrix for constraint
    for (int i = 0; i < _neG; i++) {
      iGfun[i] = _iGfun[i];
    }
    for (int i = 0; i < _neG; i++) {
      jGvar[i] = _jGvar[i];
    }
    if(_neA !=0)
    {
      // There are cases where all the constraints are nonlinear.
      for (int i = 0; i < lenA; i++) {
        iAfun[i] = _iAfun[i];
      }
      for (int i = 0; i < lenA; i++) {
        jAvar[i] = _jAvar[i];
      }
      for (int i = 0; i < lenA; i++) {
        A[i] = _A[i];
      }
    }
    neG = _neG;
    neA = _neA;
  }
  void JacobianPositionUpdate(const int& _neG, const int *_iGfun, const int *_jGvar, const int *_iAfun, const int * _jAvar, const int& _neA, const double * _A)
  {
    // This function is used to save the position of the position of the Jacobian matrix for constraint
    for (int i = 0; i < _neG; i++) {
      iGfun[i] = _iGfun[i];
    }
    for (int i = 0; i < _neG; i++) {
      jGvar[i] = _jGvar[i];
    }
    if(_neA !=0)
    {
      // There are cases where the constraints are all nonlinear.
      for (int i = 0; i < _neA; i++) {
        iAfun[i] = _iAfun[i];
      }
      for (int i = 0; i < _neA; i++) {
        jAvar[i] = _jAvar[i];
      }
      for (int i = 0; i < _neA; i++) {
        A[i] = _A[i];
      }
    }
    neG = _neG;
    neA = _neA;
  }
  void VariableBoundsUpdate(const std::vector<double> &_xlow, const std::vector<double> &_xupp)
  {
    // This function can only be called after the InnerVariableInitialize function
    // Set the upper and lower bounds.
    assert(_xlow.size() == _xupp.size());
    for (int i = 0; i < _xlow.size(); i++)
    {
      xlow[i] = _xlow[i];
      xupp[i] = _xupp[i];
      xstate[i] = 0;
    }
  }
  void ConstraintBoundsUpdate(const std::vector<double> &_Flow, const std::vector<double> &_Fupp)
  {
    // This function can only be called after the InnerGraidentInitialize function
    assert(_Flow.size() == _Fupp.size());
    for (int i = 0; i < _Flow.size(); i++)
    {
      Flow[i] = _Flow[i];
      Fupp[i] = _Fupp[i];
      Fmul[i] = 0;
    }
  }
  void SeedGuessUpdate(const std::vector<double> & _x)
  {
    // This function can only be called after n and neF have alrady been set.
    assert(n != 0);           // Make sure that n has already been evaluated.
    assert(n == _x.size());
    for (int i = 0; i < _x.size(); i++)
    {
      x[i] = _x[i];
    }
  }
  void ProblemNameUpdate(const std::string & ProName, const int& PrintFlag)
  {
    // Load the data for NonlinearProb ...
    NonlinearProb.initialize    ("", 1);      // no print file; summary on
    std::string PrintFileName = ProName + ".out";
    const char *PrintFileChar = PrintFileName.c_str();
    switch (PrintFlag) {
      case 1:
      {
        NonlinearProb.setPrintFile  (PrintFileChar); // oh wait, i want a print file
        NonlinearProb.setProbName   (ProName.c_str());
      }
      break;
      default:
      break;
    }
  }
  void ProblemOptionsUpdate(const int& DerivativeOrNot, const int& VerifyLevel)
  {
    NonlinearProb.setIntParameter("Derivative option", DerivativeOrNot);
    NonlinearProb.setIntParameter("Verify level ", VerifyLevel);
  }

  snoptProblemA NonlinearProb;

  int n, neF;

  int    ObjRow;  // = 0;
  double ObjAdd;  // = 0;

  int nS, nInf;
  double sumInf;

  // State
  double *x, *xlow, *xupp, *xmul;
  int    *xstate;

  int DerivativeFlag;

  // Constraint
  double *F, *Flow, *Fupp, *Fmul;
  int    *Fstate;

  // Gradient
  int lenA, *iAfun, *jAvar;
  double *A;
  int lenG, *iGfun, *jGvar;

  int neA, neG; // neA and neG must be defined when providing dervatives

  int StartType; //     (Cold = 0, Basis = 1, Warm = 2;)

  // The following variables are used as internal variables.
  static std::vector<LinkInfo> RobotLinkInfo;
  static SignedDistanceFieldInfo SDFInfo;
  static AnyCollisionGeometry3D TerrColGeom;
};

#endif
