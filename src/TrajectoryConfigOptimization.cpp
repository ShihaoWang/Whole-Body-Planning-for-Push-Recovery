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

struct TrajConfigOpt: public NonlinearOptimizerInfo
{
  TrajConfigOpt():NonlinearOptimizerInfo(){};

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

      std::vector<double> F_val = TrajConfigOptNCons(*n, *neF, x_vec);
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
  static std::vector<double> TrajConfigOptNCons(const int & nVar, const int & nObjNCons, const std::vector<double> & SwingLinkConfig)
  {
    // This funciton provides the constraint for the configuration variable
    std::vector<double> F(nObjNCons);
    for (int i = 0; i < SwingLinkChain.size(); i++)
        ReferenceConfig[SwingLinkChain[i]] = SwingLinkConfig[i];
    SimRobotObj.UpdateConfig(Config(ReferenceConfig));
    SimRobotObj.UpdateGeometry();
    Vector3 EndEffectorAvgPos;
    SimRobotObj.GetWorldPosition( NonlinearOptimizerInfo::RobotLinkInfo[SwingLinkInfoIndex].AvgLocalContact,
                                  NonlinearOptimizerInfo::RobotLinkInfo[SwingLinkInfoIndex].LinkIndex,
                                  EndEffectorAvgPos);
    Vector3 AvgDiff = EndEffectorAvgPos - GoalPos;
    F[0] = AvgDiff.normSquared();
    int ConstraintIndex = 1;

    // Self-collision constraint
    std::vector<double> SelfCollisionDistVec(SwingLinkChain.size()-3);
    for (int i = 0; i < SwingLinkChain.size() - 3; i++){
      Box3D Box3DObj = SimRobotObj.geometry[SwingLinkChain[i]]->GetBB();
      std::vector<Vector3> BoxVerticesVec = BoxVertices(Box3DObj);
      std::vector<double> DistVec(BoxVerticesVec.size());
      for (int j = 0; j < BoxVerticesVec.size(); j++)
        DistVec[j] = SelfLinkGeoObj.SelfCollisionDist(SwingLinkInfoIndex, BoxVerticesVec[j]);
      SelfCollisionDistVec[i] = *std::min_element(DistVec.begin(), DistVec.end());
    }
    F[ConstraintIndex] = *std::min_element(SelfCollisionDistVec.begin(), SelfCollisionDistVec.end());
    ConstraintIndex+=1;

    RobotLink3D EndEffectorLink = SimRobotObj.links[NonlinearOptimizerInfo::RobotLinkInfo[SwingLinkInfoIndex].LinkIndex];

    Vector3 EndEffectorDir;
    EndEffectorDir.x = EndEffectorLink.T_World.R.data[2][0];
    EndEffectorDir.y = EndEffectorLink.T_World.R.data[2][1];
    EndEffectorDir.z = EndEffectorLink.T_World.R.data[2][2];

    double projValue = EndEffectorDir.dot(GoalDir);
    F[ConstraintIndex] = projValue - alignmentVal ;
    ConstraintIndex+=1;

    return F;
  }
};

std::vector<double> TrajConfigOptimazation(const Robot & SimRobot, ReachabilityMap & RMObject, SelfLinkGeoInfo & _SelfLinkGeoObj, SimPara & SimParaObj, const double & _alignmentVal, const int & StageIndex){
  // This function is used to optimize robot's configuration such that a certain contact can be reached for that end effector.
  SimRobotObj = SimRobot;
  SwingLinkInfoIndex = SimParaObj.getSwingLinkInfoIndex();
  SwingLinkChain = RMObject.EndEffectorLink2Pivotal[SwingLinkInfoIndex];
  GoalPos = SimParaObj.getCurrentContactPos();
  GoalDir = SimParaObj.getGoalDirection();
  ReferenceConfig = SimRobot.q;
  alignmentVal = _alignmentVal;
  SelfLinkGeoObj = _SelfLinkGeoObj;

  TrajConfigOpt TrajConfigOptProblem;

  // Static Variable Substitution
  std::vector<double> SwingLinkChainGuess(SwingLinkChain.size());
  int n = SwingLinkChain.size();

  // Cost function on the norm difference between the reference avg position and the modified contact position.
  int neF = 1;
  neF += 1;                                                                                           // Self-Collision Avoidance
  neF += 2;                                                                                           // End effector alignment expressed as two inequality constrains
  TrajConfigOptProblem.InnerVariableInitialize(n, neF);

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
  TrajConfigOptProblem.VariableBoundsUpdate(xlow_vec, xupp_vec);

  /*
    Initialize the bounds of variables
  */
  std::vector<double> Flow_vec(neF), Fupp_vec(neF);
  for (int i = 0; i < neF; i++)
  {
    Flow_vec[i] = 0;
    Fupp_vec[i] = 1e10;
  }
  for (int i = neF-2; i < neF; i++)
  {
    Flow_vec[i] = -alignmentTol;
    Fupp_vec[i] = alignmentTol;
  }

  TrajConfigOptProblem.ConstraintBoundsUpdate(Flow_vec, Fupp_vec);

  /*
    Initialize the seed guess
  */
  TrajConfigOptProblem.SeedGuessUpdate(SwingLinkChainGuess);

  /*
    Given a name of this problem for the output
  */
  TrajConfigOptProblem.ProblemNameUpdate("TrajConfigOptProblem", 0);

  // Here we would like allow much more time to be spent on IK
  TrajConfigOptProblem.NonlinearProb.setIntParameter("Iterations limit", 250);
  TrajConfigOptProblem.NonlinearProb.setIntParameter("Major iterations limit", 25);
  TrajConfigOptProblem.NonlinearProb.setIntParameter("Major print level", 0);
  TrajConfigOptProblem.NonlinearProb.setIntParameter("Minor print level", 0);
  /*
    ProblemOptions seting
  */
  // Solve with Finite-Difference
  TrajConfigOptProblem.ProblemOptionsUpdate(0, 3);
  TrajConfigOptProblem.Solve(SwingLinkChainGuess);

  std::vector<double> OptConfig = ReferenceConfig;

  for (int i = 0; i < n; i++)
  {
    OptConfig[SwingLinkChain[i]] = SwingLinkChainGuess[i];
  }
  SimRobotObj.UpdateConfig(Config(OptConfig));
  SimRobotObj.UpdateGeometry();

  std::string ConfigPath = "/home/motion/Desktop/Online-Contact-Planning-for-Fall-Mitigation/user/hrp2/";
  string _OptConfigFile = to_string(StageIndex) + "InnerOptConfig.config";
  RobotConfigWriter(OptConfig, ConfigPath, _OptConfigFile);

  // Self-collision constraint numerical checker
  std::vector<double> SelfCollisionDistVec(SwingLinkChain.size()-3);
  for (int i = 0; i < SwingLinkChain.size()-3; i++)     // Due to the bounding box size of torso link
  {
    Box3D Box3DObj = SimRobotObj.geometry[SwingLinkChain[i]]->GetBB();
    std::vector<Vector3> BoxVerticesVec = BoxVertices(Box3DObj);
    std::vector<double> DistVec(BoxVerticesVec.size());
    for (int j = 0; j < BoxVerticesVec.size(); j++)
      DistVec[j] = _SelfLinkGeoObj.SelfCollisionDist(SwingLinkInfoIndex, BoxVerticesVec[j]);
    SelfCollisionDistVec[i] = *std::min_element(DistVec.begin(), DistVec.end());
  }
  double SelfCollisionDistTol = *std::min_element(SelfCollisionDistVec.begin(), SelfCollisionDistVec.end());

  bool OptFlag = true;
  if(SelfCollisionDistTol<-0.0025){
      std::printf("Transient Optimization Failure due to Self-collision for Link %d! \n", NonlinearOptimizerInfo::RobotLinkInfo[SwingLinkInfoIndex].LinkIndex);
      OptFlag = false;
  }

  Vector3 EndEffectorAvgPos;
  SimRobotObj.GetWorldPosition(NonlinearOptimizerInfo::RobotLinkInfo[SwingLinkInfoIndex].AvgLocalContact, NonlinearOptimizerInfo::RobotLinkInfo[SwingLinkInfoIndex].LinkIndex, EndEffectorAvgPos);
  Vector3 AvgDiff = EndEffectorAvgPos - GoalPos;
  double DistTestTol = 0.0225;
  double DistTest = AvgDiff.normSquared();
  if(DistTest>DistTestTol){
    std::printf("Transient Optimization Failure due to Goal Contact Non-reachability for Link %d! \n", NonlinearOptimizerInfo::RobotLinkInfo[SwingLinkInfoIndex].LinkIndex);
    OptFlag = false;
  }
  SimParaObj.setTrajConfigOptFlag(OptFlag);
  OptConfig = YPRShifter(OptConfig);
  return OptConfig;
}
