#include "CommonHeader.h"
#include "RobotInfo.h"
#include "NonlinearOptimizerInfo.h"

static double EdgeProjMagnitude(const double & cur_s,  const Vector3 & init_dir, const Vector3 & goal_dir){
  // This function calculates the projection based on current s'length.
  double xProj = init_dir.dot(goal_dir);
  return (1.0 - cur_s) * xProj;
}

static Robot SimRobotObj;
static int SwingLinkInfoIndex;
static std::vector<int> SwingLinkChain;

static Matrix Jac;
static Vector3 GoalDir;

static std::vector<double> CurrentConfig;
static std::vector<double> CurrentVelocity;
static std::vector<double> CurrentBaseDelta(6);

static Vector3 EndEffectorInitxDir;
static Vector3 EndEffectorInityDir;

static double AlignmentTol = 1e-4;
static double EndEffectorTol = 1e-4;
static double K_delta_s = 10.0;
static double K_goal = 10.0;
static double K_fb = 0.001;
static double K_E = 100.0;
static double K_A = 100.0;

static int DOF;

static Vector3 d_x_s;
static double sCur;
static double delta_t;

struct TaskSpaceOpt: public NonlinearOptimizerInfo
{
  TaskSpaceOpt():NonlinearOptimizerInfo(){};

  // This struct inherits the NonlinearOptimizerInfo struct and we just need to defined the Constraint function
  static void ObjNConstraint(int    *Status, int *n,    double x[],
    int    *needF,  int *neF,  double F[],
    int    *needG,  int *neG,  double G[],
    char      *cu,  int *lencu,
    int    iu[],    int *leniu,
    double ru[],    int *lenru)
    {
      std::vector<double> x_vec(*n);
      for (int i = 0; i < *n; i++)
        x_vec[i] = x[i];

      std::vector<double> F_val = TaskSpaceOptNCons(*n, *neF, x_vec);
      for (int i = 0; i < *neF; i++)
        F[i] = F_val[i];
    }
  void Solve(std::vector<double> &RobotConfig)
  {
    int StartType = 0;
      NonlinearProb.solve(StartType, neF, n, ObjAdd, ObjRow, ObjNConstraint,
      xlow, xupp, Flow, Fupp,
      x, xstate, xmul, F, Fstate, Fmul,
      nS, nInf, sumInf);
      for (int i = 0; i < n; i++)
        RobotConfig[i] = x[i];

      delete []x;      delete []xlow;   delete []xupp;
      delete []xmul;   delete []xstate;

      delete []F;      delete []Flow;   delete []Fupp;
      delete []Fmul;   delete []Fstate;
  }
  static std::vector<double> TaskSpaceOptNCons(const int & nVar, const int & nObjNCons, const std::vector<double> & AccelerationPluss)
  {
    // This funciton provides the constraint for the configuration variable
    std::vector<double> F(nObjNCons);

    double delta_s = AccelerationPluss[SwingLinkChain.size()];
    double Ex = AccelerationPluss[SwingLinkChain.size() + 1];
    double Ey = AccelerationPluss[SwingLinkChain.size() + 2];
    double Ez = AccelerationPluss[SwingLinkChain.size() + 3];
    double Ax = AccelerationPluss[SwingLinkChain.size() + 4];
    double Ay = AccelerationPluss[SwingLinkChain.size() + 5];

    F[0] = -K_delta_s * delta_s + K_E * (Ex * Ex + Ey * Ey + Ez * Ez) + K_A * (Ax * Ax + Ay * Ay);

    // Velocity Constraint
    int ConstraintIndex = 1;
    std::vector<double> NextVelocity(SwingLinkChain.size());
    for (int i = 0; i < SwingLinkChain.size(); i++) {
      NextVelocity[i] = CurrentVelocity[SwingLinkChain[i]] + AccelerationPluss[i] * delta_t;
      F[ConstraintIndex] = NextVelocity[i];
      ConstraintIndex++;
    }

    std::vector<double> NextConfig = CurrentConfig;
    for (int i = 0; i < 6; i++)
      NextConfig[i]+=CurrentBaseDelta[i];
    for (int i = 0; i < SwingLinkChain.size(); i++) {
      double CurrentSwingLinkDelta = 0.5 * (CurrentVelocity[SwingLinkChain[i]] + NextVelocity[i]) * delta_t;
      NextConfig[SwingLinkChain[i]]+=CurrentSwingLinkDelta;
    }
    for (int i = 0; i < SwingLinkChain.size(); i++) {
      F[ConstraintIndex] = NextConfig[SwingLinkChain[i]];
      ConstraintIndex++;
    }

    std::vector<double> DeltaQ(DOF);
    for (int i = 0; i < DOF; i++)
      DeltaQ[i] = NextConfig[i] - CurrentConfig[i];

    // Vector3 Jac_delta_q;
    for (int i = 0; i < 3; i++) {
      double Jac_delta_q_i = 0.0;
      for (int j = 0; j < DOF; j++) {
        Jac_delta_q_i+=Jac(i,j) * DeltaQ[j];
      }
      double Jac_delta_q_delta_s = Jac_delta_q_i - delta_s * d_x_s[i] - AccelerationPluss[SwingLinkChain.size() + i + 1];
      F[ConstraintIndex] = Jac_delta_q_delta_s * Jac_delta_q_delta_s;
      ConstraintIndex++;
    }

    SimRobotObj.UpdateConfig(Config(NextConfig));

    double sNew = sCur + delta_s;

    double EdgexProjVal = EdgeProjMagnitude(sNew,  EndEffectorInitxDir,   GoalDir);
    double EdgeyProjVal = EdgeProjMagnitude(sNew,  EndEffectorInityDir,   GoalDir);

    Vector3 EndEffectorxDir, EndEffectoryDir;
    getEndEffectorXYAxes(SimRobotObj, SwingLinkInfoIndex, EndEffectorxDir, EndEffectoryDir);
    F[ConstraintIndex] = EndEffectorxDir.dot(GoalDir) - EdgexProjVal + Ax;
    ConstraintIndex+=1;
    F[ConstraintIndex] = EndEffectoryDir.dot(GoalDir) - EdgeyProjVal + Ay;
    ConstraintIndex+=1;

    return F;
  }
};

void TaskSpaceOptimazation( double & delta_s, std::vector<double> & Acceleration, std::vector<double> & SlackVec, const double & DampingRatio){
  // This function is used to optimize robot's velocity to reach a target as soon as possible.

  TaskSpaceOpt TaskSpaceOptProblem;

  // Static Variable Substitution
  int n = SwingLinkChain.size() + 1 + 3 + 2;      // acceleration + delta_s + Ex, Ey, Ez + Ax, Ay

  // Cost function on the norm difference between the reference avg position and the modified conta`ct position.
  int neF = 1;
  neF += 2 * SwingLinkChain.size();                                                 // Velocity
  neF += 3;                                                                         // Alignment
  neF += 2;                                                                         // End Effector Orientation Constraint
  TaskSpaceOptProblem.InnerVariableInitialize(n, neF);

  /*
    Initialize the bounds of variables
  */
  std::vector<double> xlow_vec(n), xupp_vec(n);
  for (int i = 0; i < SwingLinkChain.size(); i++){
    // Velocity
    xlow_vec[i] = -SimRobotObj.accMax(SwingLinkChain[i]);
    xupp_vec[i] = SimRobotObj.accMax(SwingLinkChain[i]);
  }
  xlow_vec[SwingLinkChain.size()] = 1e-3;
  xupp_vec[SwingLinkChain.size()] = 0.25;
  double SlackBound = 1e-3;
  for (int i = n - 5; i < n; i++) {
    xlow_vec[i] = -SlackBound;
    xupp_vec[i] = SlackBound;
  }
  TaskSpaceOptProblem.VariableBoundsUpdate(xlow_vec, xupp_vec);

  /*
    Initialize the bounds of variables
  */
  std::vector<double> Flow_vec(neF), Fupp_vec(neF);
  // Objective: delta_s
  Flow_vec[0] = -1e10;
  Fupp_vec[0] = 1e10;
  int neFCount = 1;
  // Constraint: velocity
  for (int i = 0; i <SwingLinkChain.size(); i++) {
    double velMax = SimRobotObj.velMax(SwingLinkChain[i]);
    Flow_vec[neFCount] = -velMax * DampingRatio;
    Fupp_vec[neFCount] = velMax * DampingRatio;
    neFCount++;
  }
  // Constraint: configuration
  for (int i = 0; i <SwingLinkChain.size(); i++) {
    Flow_vec[neFCount] = SimRobotObj.qMin(SwingLinkChain[i]);
    Fupp_vec[neFCount] = SimRobotObj.qMax(SwingLinkChain[i]);
    neFCount++;
  }
  // Constraint:: Alignment
  for (int i = 0; i < 3; i++) {
    Flow_vec[neFCount] = -AlignmentTol * AlignmentTol;
    Fupp_vec[neFCount] = AlignmentTol * AlignmentTol;
    neFCount++;
  }
  Flow_vec[neF-2] = -1.0 * EndEffectorTol;
  Fupp_vec[neF-2] = EndEffectorTol;
  Flow_vec[neF-1] = -1.0 * EndEffectorTol;
  Fupp_vec[neF-1] = EndEffectorTol;

  TaskSpaceOptProblem.ConstraintBoundsUpdate(Flow_vec, Fupp_vec);

  /*
    Initialize the seed guess
  */
  std::vector<double> SeedGuess = Acceleration;
  double SlackInit = 0.1;
  SeedGuess.push_back(0.1);
  for (int i = 0; i < 5; i++)
    SeedGuess.push_back(SlackInit);

  TaskSpaceOptProblem.SeedGuessUpdate(SeedGuess);

  /*
    Given a name of this problem for the output
  */
  TaskSpaceOptProblem.ProblemNameUpdate("TaskSpaceOptProblem", 0);

  // Here we would like allow much more time to be spent on IK
  TaskSpaceOptProblem.NonlinearProb.setIntParameter("Iterations limit", 1000);
  TaskSpaceOptProblem.NonlinearProb.setIntParameter("Major iterations limit", 100);
  TaskSpaceOptProblem.NonlinearProb.setIntParameter("Major print level", 0);
  TaskSpaceOptProblem.NonlinearProb.setIntParameter("Minor print level", 0);
  /*
    ProblemOptions seting
  */
  // Solve with Finite-Difference
  TaskSpaceOptProblem.ProblemOptionsUpdate(0, 3);
  TaskSpaceOptProblem.Solve(SeedGuess);

  delta_s = SeedGuess[SwingLinkChain.size()];
  for (int i = 0; i < SwingLinkChain.size(); i++)
    Acceleration[i] = SeedGuess[i];
  for (int i = 0; i < 5; i++)
    SlackVec.push_back(SeedGuess[SwingLinkChain.size() + i + 1]);
  return;
}

bool TaskTrajectoryPlanningInner( const double & _sVal, double & _sNew,
                                  Robot & _SimRobotInner,
                                  const Config & _PreVelocity,                      const std::vector<double> & _CurrentBaseDelta,
                                  const Config & _CurrentConfig,                    const Config & _CurrentVelocity,
                                  std::vector<double> & _NextConfig,                std::vector<double> & _NextVelocity,
                                  const EndEffectorPathInfo & EndEffectorPath,      const std::vector<int> & _SwingLinkChain,
                                  const Vector3 & _EndEffectorInitxDir,             const Vector3 & _EndEffectorInityDir,
                                  SimPara & _SimParaObj,
                                  const double & _StageTime,                        const double & _DampingRatio){

  SimRobotObj = _SimRobotInner;
  SwingLinkInfoIndex = _SimParaObj.getSwingLinkInfoIndex();
  SwingLinkChain = _SwingLinkChain;
  GoalDir = _SimParaObj.getGoalDirection();

  DOF = _SimRobotInner.q.size();

  CurrentConfig = _CurrentConfig;
  CurrentVelocity = _CurrentVelocity;
  CurrentBaseDelta = _CurrentBaseDelta;

  EndEffectorInitxDir = _EndEffectorInitxDir;
  EndEffectorInityDir = _EndEffectorInityDir;

  delta_t = _StageTime;

  Vector3 CurrentContactPos;
  _SimRobotInner.GetWorldPosition(    NonlinearOptimizerInfo::RobotLinkInfo[SwingLinkInfoIndex].AvgLocalContact,
                                      NonlinearOptimizerInfo::RobotLinkInfo[SwingLinkInfoIndex].LinkIndex, CurrentContactPos);
  _SimRobotInner.GetPositionJacobian( NonlinearOptimizerInfo::RobotLinkInfo[SwingLinkInfoIndex].AvgLocalContact,
                                      NonlinearOptimizerInfo::RobotLinkInfo[SwingLinkInfoIndex].LinkIndex, Jac);

  sCur = _sVal;
  Vector3 CurrentsPos, CurrentContactsSlope;
  EndEffectorPathInfo EndEffectorPathInner = EndEffectorPath;
  EndEffectorPathInner.PosNTang(sCur, CurrentsPos, CurrentContactsSlope);
  Vector3 FeedbackPos = CurrentsPos - CurrentContactPos;
  // d_x_s = CurrentContactsSlope + K_goal *  (_SimParaObj.getContactGoal() - CurrentContactPos) + K_fb * FeedbackPos;
  d_x_s = K_goal *  (_SimParaObj.getContactGoal() - CurrentContactPos);

  double delta_s = 0.25;      // Initial Guess
  std::vector<double> NextConfig    = CurrentConfig;
  std::vector<double> NextVelocity  = CurrentVelocity;
  std::vector<double> Acceleration(SwingLinkChain.size());   // Note that here the acceleration is only for the swing link!
  for (int i = 0; i < Acceleration.size(); i++)
    Acceleration[i] = (_CurrentVelocity[SwingLinkChain[i]] - _PreVelocity[SwingLinkChain[i]])/_StageTime;

  std::vector<double> SlackVec;
  TaskSpaceOptimazation(delta_s, Acceleration, SlackVec, _DampingRatio);
  _sNew = _sVal + delta_s;
  // Update current configuration
  for (int i = 0; i < 6; i++)
    NextConfig[i]+=_CurrentBaseDelta[i];
  for (int i = 0; i < SwingLinkChain.size(); i++) {
    NextConfig[SwingLinkChain[i]]+= CurrentVelocity[SwingLinkChain[i]] * delta_t + 0.5 * Acceleration[i] * delta_t * delta_t;
  }
  // Update current velocity
  for (int i = 0; i < 6; i++)
    NextVelocity[i] =_CurrentBaseDelta[i]/delta_t;
  for (int i = 0; i < SwingLinkChain.size(); i++)
    NextVelocity[SwingLinkChain[i]] = CurrentVelocity[SwingLinkChain[i]] + Acceleration[i] * delta_t;

  // double Ex = SlackVec[0];
  // double Ey = SlackVec[1];
  // double Ez = SlackVec[2];
  // double Ax = SlackVec[3];
  // double Ay = SlackVec[4];
  //
  // std::vector<double> DeltaQ(DOF);
  // for (int i = 0; i < DOF; i++)
  //   DeltaQ[i] = NextConfig[i] - CurrentConfig[i];
  //
  // // Vector3 Jac_delta_q;
  // for (int i = 0; i < 3; i++) {
  //   double Jac_delta_q_i = 0.0;
  //   for (int j = 0; j < DOF; j++) {
  //     Jac_delta_q_i+=Jac(i,j) * DeltaQ[j];
  //   }
  //   double Jac_delta_q = Jac_delta_q_i - delta_s * d_x_s[i] - SlackVec[i];
  //   std::cout<<"Jac_delta_q : " << Jac_delta_q << std::endl;
  // }
  //
  // double EdgexProjVal = EdgeProjMagnitude(_sNew,  EndEffectorInitxDir,   GoalDir);
  // double EdgeyProjVal = EdgeProjMagnitude(_sNew,  EndEffectorInityDir,   GoalDir);
  //
  // Vector3 EndEffectorxDir, EndEffectoryDir;
  // SimRobotObj.UpdateConfig(Config(NextConfig));
  //
  // getEndEffectorXYAxes(SimRobotObj, SwingLinkInfoIndex, EndEffectorxDir, EndEffectoryDir);
  // double AxDiff = EndEffectorxDir.dot(GoalDir) - EdgexProjVal + Ax;
  // double AyDiff  = EndEffectoryDir.dot(GoalDir) - EdgeyProjVal + Ay;
  //
  // std::cout<<"AxDiff : " <<AxDiff<<std::endl;
  // std::cout<<"AyDiff : " <<AyDiff<<std::endl;

  _NextConfig = NextConfig;
  _NextVelocity = NextVelocity;

  return true;
}
