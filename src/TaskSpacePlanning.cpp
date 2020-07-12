#include "CommonHeader.h"
#include "RobotInfo.h"
#include "NonlinearOptimizerInfo.h"

static Robot SimRobotObj;
static int SwingLinkInfoIndex;
static std::vector<int> SwingLinkChain;
static Matrix Jac;
static Vector3 GoalDir;
static std::vector<double> CurrentConfig;
static std::vector<double> CurrentVelocity;
static std::vector<double> NextVelocity;
static SelfLinkGeoInfo SelfLinkGeoObj;
static Vector3 EndEffectorInitxDir;
static Vector3 EndEffectorInityDir;

static EndEffectorPathInfo EndEffectorPathObj;

static double AlignmentTol = 1e-2;
static double EndEffectorTol = 1e-3;
static Vector3 d_q_s;
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
  static std::vector<double> TaskSpaceOptNCons(const int & nVar, const int & nObjNCons, const std::vector<double> & SwingLinkVelocityPluss)
  {
    // This funciton provides the constraint for the configuration variable
    std::vector<double> F(nObjNCons);
    for (int i = 0; i < SwingLinkChain.size(); i++)
        NextVelocity[SwingLinkChain[i]] = SwingLinkVelocityPluss[i];
    double delta_s = SwingLinkVelocityPluss.back();
    F[0] = -delta_s;

    // Acceleration Constraint
    int ConstraintIndex = 1;
    for (int i = 0; i < SwingLinkChain.size(); i++){
      F[ConstraintIndex] = (SwingLinkVelocityPluss[i] - CurrentVelocity[SwingLinkChain[i]])/delta_t;
      ConstraintIndex++;
    }

    std::vector<double> NextConfig = CurrentConfig;
    std::vector<double> delta_q(NextConfig.size());
    for (int i = 0; i < SwingLinkChain.size(); i++){
      double delta_q_i = 0.5 * delta_t * (CurrentVelocity[SwingLinkChain[i]] + NextVelocity[i]);
      delta_q[SwingLinkChain[i]] = delta_q_i;
      NextConfig[SwingLinkChain[i]] += delta_q_i;
    }
    // Vector3 Jac_delta_q;
    for (int i = 0; i < 3; i++) {
      double Jac_delta_q_i = 0.0;
      for (int j = 0; j < NextConfig.size(); j++) {
        Jac_delta_q_i+=Jac(i,j) * delta_q[j];
      }
      // Jac_delta_q[i] = Jac_delta_q_i;
      F[ConstraintIndex] = Jac_delta_q_i - delta_s * d_q_s[i];
      ConstraintIndex++;
    }
    SimRobotObj.UpdateConfig(Config(NextConfig));

    double sNew = sCur + delta_s;

    double EdgexProjVal = EdgeProjMagnitude(sNew,  EndEffectorInitxDir, GoalDir);
    double EdgeyProjVal = EdgeProjMagnitude(sNew,  EndEffectorInityDir, GoalDir);

    RobotLink3D EndEffectorLink = SimRobotObj.links[NonlinearOptimizerInfo::RobotLinkInfo[SwingLinkInfoIndex].LinkIndex];

    Vector3 EndEffectorxDir, EndEffectoryDir;   // Eventually these two directions should be orthgonal to goal direction.
    EndEffectorxDir.x = EndEffectorLink.T_World.R.data[0][0];
    EndEffectorxDir.y = EndEffectorLink.T_World.R.data[0][1];
    EndEffectorxDir.z = EndEffectorLink.T_World.R.data[0][2];
    F[ConstraintIndex] = EndEffectorxDir.dot(GoalDir) - EdgexProjVal;
    ConstraintIndex+=1;

    EndEffectoryDir.x = EndEffectorLink.T_World.R.data[1][0];
    EndEffectoryDir.y = EndEffectorLink.T_World.R.data[1][1];
    EndEffectoryDir.z = EndEffectorLink.T_World.R.data[1][2];
    F[ConstraintIndex] = EndEffectoryDir.dot(GoalDir) - EdgeyProjVal;
    ConstraintIndex+=1;

    return F;
  }
};

std::vector<double> TaskSpaceOptimazation(const Robot & SimRobot, double & delta_s){
  // This function is used to optimize robot's velocity to reach a target as soon as possible.

  TaskSpaceOpt TaskSpaceOptProblem;

  // Static Variable Substitution
  std::vector<double> SwingLinkChainVelocityGuess(SwingLinkChain.size());
  int n = SwingLinkChain.size() + 1;      // velocity + delta_s

  // Cost function on the norm difference between the reference avg position and the modified contact position.
  int neF = 1;
  neF += SwingLinkChain.size();                                                 // Acceleration
  neF += 3;                                                                     // Alignment
  neF += 2;                                                                     // End Effector Orientation Constraint
  TaskSpaceOptProblem.InnerVariableInitialize(n, neF);

  /*
    Initialize the bounds of variables
  */
  std::vector<double> xlow_vec(n), xupp_vec(n);
  for (int i = 0; i < n-1; i++)
  {
    // Velocity
    xlow_vec[i] = SimRobot.velMin(SwingLinkChain[i]);
    xupp_vec[i] = SimRobot.velMax(SwingLinkChain[i]);
    SwingLinkChainVelocityGuess[i] = CurrentVelocity[SwingLinkChain[i]];
  }
  xlow_vec[n-1] = 1e-3;
  xupp_vec[n-1] = 1.0;
  TaskSpaceOptProblem.VariableBoundsUpdate(xlow_vec, xupp_vec);

  /*
    Initialize the bounds of variables
  */
  std::vector<double> Flow_vec(neF), Fupp_vec(neF);
  // Objective: delta_s
  Flow_vec[0] = 1e-3;
  Fupp_vec[0] = 1.0;
  int neFCount = 1;
  // Constraint: acceleartions
  for (int i = 0; i <SwingLinkChain.size(); i++) {
    double accMax = SimRobot.accMax(SwingLinkChain[i]);
    Flow_vec[neFCount] = -accMax;
    Fupp_vec[neFCount] = accMax;
    neFCount++;
  }
  // Constraint:: Alignment
  for (int i = 0; i < 3; i++) {
    Flow_vec[neFCount] = -1.0 * AlignmentTol;
    Fupp_vec[neFCount] = AlignmentTol;
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
  std::vector<double> SeedGuess = SwingLinkChainVelocityGuess;
  SeedGuess.push_back(0.1);
  TaskSpaceOptProblem.SeedGuessUpdate(SeedGuess);

  /*
    Given a name of this problem for the output
  */
  TaskSpaceOptProblem.ProblemNameUpdate("TaskSpaceOptProblem", 0);

  // Here we would like allow much more time to be spent on IK
  TaskSpaceOptProblem.NonlinearProb.setIntParameter("Iterations limit", 250);
  TaskSpaceOptProblem.NonlinearProb.setIntParameter("Major iterations limit", 25);
  TaskSpaceOptProblem.NonlinearProb.setIntParameter("Major print level", 0);
  TaskSpaceOptProblem.NonlinearProb.setIntParameter("Minor print level", 0);
  /*
    ProblemOptions seting
  */
  // Solve with Finite-Difference
  TaskSpaceOptProblem.ProblemOptionsUpdate(0, 3);
  TaskSpaceOptProblem.Solve(SeedGuess);

  std::vector<double> NextConfig = CurrentConfig;
  std::vector<double> delta_q(CurrentConfig.size());
  for (int i = 0; i < SwingLinkChain.size(); i++){
    double delta_q_i = 0.5 * delta_t * (CurrentVelocity[SwingLinkChain[i]] + SeedGuess[i]);
    NextConfig[SwingLinkChain[i]] += delta_q_i;
    NextVelocity[SwingLinkChain[i]] = SeedGuess[i];
  }
  SimRobotObj.UpdateConfig(Config(NextConfig));

  NextConfig = YPRShifter(NextConfig);
  return NextConfig;
}

bool TaskSpacePlanning(const double & sVal, double & sNew, Robot & SimRobotInner,
                              const Config & _CurrentConfig, const Config & _CurrentVelocity,
                              std::vector<double> & _NextConfig, std::vector<double> & _NextVelocity,
                              const EndEffectorPathInfo & EndEffectorPathInner, const std::vector<int> & _SwingLinkChain,
                              const Vector3 & _EndEffectorInitxDir, const Vector3 & _EndEffectorInityDir,
                              SimPara & SimParaObj, const double & StageTime){

  SimRobotObj = SimRobotInner;
  SwingLinkInfoIndex = SimParaObj.getSwingLinkInfoIndex();
  SwingLinkChain = _SwingLinkChain;
  GoalDir = SimParaObj.getGoalDirection();

  CurrentConfig = _CurrentConfig;
  CurrentVelocity = _CurrentVelocity;
  NextVelocity = CurrentVelocity;

  EndEffectorInitxDir = _EndEffectorInitxDir;
  EndEffectorInityDir = _EndEffectorInityDir;

  delta_t = StageTime;

  Vector3 CurrentContactPos;
  SimRobotInner.GetWorldPosition( NonlinearOptimizerInfo::RobotLinkInfo[SwingLinkInfoIndex].AvgLocalContact,
                                  NonlinearOptimizerInfo::RobotLinkInfo[SwingLinkInfoIndex].LinkIndex, CurrentContactPos);
  Matrix JacMat;
  SimRobotInner.GetPositionJacobian(NonlinearOptimizerInfo::RobotLinkInfo[SwingLinkInfoIndex].AvgLocalContact,
                                    NonlinearOptimizerInfo::RobotLinkInfo[SwingLinkInfoIndex].LinkIndex, JacMat);
  Jac = JacMat;

  EndEffectorPathObj = EndEffectorPathInner;
  sCur = sVal;
  Vector3 CurrentsPos, CurrentContactsSlope;
  EndEffectorPathObj.PosNTang(sCur, CurrentsPos, CurrentContactsSlope);
  double Ks = 1e-2;
  d_q_s = CurrentContactsSlope + Ks * (SimParaObj.getContactGoal() - CurrentsPos);

  double delta_s;
  _NextConfig = TaskSpaceOptimazation(SimRobotInner, delta_s);
  _NextVelocity = NextVelocity;
  sNew = sVal + delta_s;
  return false;
}
