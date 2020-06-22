#include "CommonHeader.h"
#include "NonlinearOptimizerInfo.h"

static Robot SimRobotObj;
static int SwingLinkInfoIndex;
static std::vector<int> SwingLinkChain;
static Vector3 GoalPos;
static Vector3 GoalDir;
static std::vector<double> ReferenceConfig;
static SelfLinkGeoInfo SelfLinkGeoObj;
static double alignmentVal = 0.0;
static double alignmentTol= 1e-5;

struct OrientationOpt: public NonlinearOptimizerInfo
{
  OrientationOpt():NonlinearOptimizerInfo(){};

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

      std::vector<double> F_val = OrientationOptNCons(*n, *neF, x_vec);
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
  static std::vector<double> OrientationOptNCons(const int & nVar, const int & nObjNCons, const std::vector<double> & SwingLinkConfig)
  {
    // This funciton provides the constraint for the configuration variable
    std::vector<double> F(nObjNCons);
    for (int i = 0; i < nVar; i++)
        ReferenceConfig[SwingLinkChain[i]] = SwingLinkConfig[i];
    SimRobotObj.UpdateConfig(Config(ReferenceConfig));
    SimRobotObj.UpdateGeometry();
    RobotLink3D EndEffectorLink = SimRobotObj.links[RobotLinkInfo[SwingLinkInfoIndex].LinkIndex];
    Vector3 EndEffectorDir;
    EndEffectorDir.x = EndEffectorLink.T_World.R.data[2][0];
    EndEffectorDir.y = EndEffectorLink.T_World.R.data[2][1];
    EndEffectorDir.z = EndEffectorLink.T_World.R.data[2][2];

    double projValue = EndEffectorDir.dot(GoalDir);
    F[0] = (projValue - alignmentVal) * (projValue - alignmentVal);

    int ConstraintIndex = 1;
    for (int i = 0; i < RobotLinkInfo[SwingLinkInfoIndex].LocalContacts.size(); i++)
    {
      Vector3 LinkiPjPos;
      SimRobotObj.GetWorldPosition(RobotLinkInfo[SwingLinkInfoIndex].LocalContacts[i], RobotLinkInfo[SwingLinkInfoIndex].LinkIndex, LinkiPjPos);
      F[ConstraintIndex] = SDFInfo.SignedDistance(LinkiPjPos);
      ConstraintIndex = ConstraintIndex + 1;
    }
    return F;
  }
};

std::vector<double> OrientationOptimazation(const Robot & SimRobot, const std::vector<int> & _SwingLinkChain, const SimPara & SimParaObj, const double & _alignmentVal, const int & StageIndex){
  // This function is used to optimize robot's configuration such that a certain contact can be reached for that end effector.
  SimRobotObj = SimRobot;
  SwingLinkInfoIndex = SimParaObj.SwingLinkInfoIndex;
  SwingLinkChain = _SwingLinkChain;
  GoalDir = SimParaObj.DirectionGoal;
  ReferenceConfig = SimRobot.q;
  alignmentVal = _alignmentVal;

  OrientationOpt OrientationOptProblem;

  // Static Variable Substitution
  int n = 2;
  std::vector<double> SwingLinkChainGuess(n);

  // Cost function on the norm difference between the reference avg position and the modified contact position.
  int neF = 1;
  neF += NonlinearOptimizerInfo::RobotLinkInfo[SwingLinkInfoIndex].LocalContacts.size();

  OrientationOptProblem.InnerVariableInitialize(n, neF);
  /*
    Initialize the bounds of variables
  */
  std::vector<double> xlow_vec(n), xupp_vec(n);
  for (int i = 0; i < n; i++)
  {
    // Configuration
    xlow_vec[i] = SimRobot.qMin(SwingLinkChain[i]);
    xupp_vec[i] = SimRobot.qMax(SwingLinkChain[i]);
    SwingLinkChainGuess[i] = ReferenceConfig[SwingLinkChain[i]];
  }
  OrientationOptProblem.VariableBoundsUpdate(xlow_vec, xupp_vec);

  /*
    Initialize the bounds of variables
  */
  std::vector<double> Flow_vec(neF), Fupp_vec(neF);
  for (int i = 0; i < neF; i++)
  {
    Flow_vec[i] = 0;
    Fupp_vec[i] = 1e10;
  }
  OrientationOptProblem.ConstraintBoundsUpdate(Flow_vec, Fupp_vec);

  /*
    Initialize the seed guess
  */
  OrientationOptProblem.SeedGuessUpdate(SwingLinkChainGuess);

  /*
    Given a name of this problem for the output
  */
  OrientationOptProblem.ProblemNameUpdate("OrientationOptProblem", 0);

  // Here we would like allow much more time to be spent on IK
  OrientationOptProblem.NonlinearProb.setIntParameter("Iterations limit", 250);
  OrientationOptProblem.NonlinearProb.setIntParameter("Major iterations limit", 25);
  OrientationOptProblem.NonlinearProb.setIntParameter("Major print level", 0);
  OrientationOptProblem.NonlinearProb.setIntParameter("Minor print level", 0);
  /*
    ProblemOptions seting
  */
  // Solve with Finite-Difference
  OrientationOptProblem.ProblemOptionsUpdate(0, 3);
  OrientationOptProblem.Solve(SwingLinkChainGuess);

  std::vector<double> OptConfig = ReferenceConfig;

  for (int i = 0; i < n; i++)
  {
    OptConfig[SwingLinkChain[i]] = SwingLinkChainGuess[i];
  }
  SimRobotObj.UpdateConfig(Config(OptConfig));
  SimRobotObj.UpdateGeometry();
  OptConfig = YPRShifter(OptConfig);

  std::string ConfigPath = "/home/motion/Desktop/Whole-Body-Planning-for-Push-Recovery/build/";
  std::string OptConfigFile = "OptWholeBodyUpdatedConfig" + std::to_string(StageIndex) + ".config";
  RobotConfigWriter(OptConfig, ConfigPath, OptConfigFile);

  return OptConfig;
}
