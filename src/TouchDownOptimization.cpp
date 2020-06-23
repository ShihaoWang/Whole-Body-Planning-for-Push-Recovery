#include "CommonHeader.h"
#include "NonlinearOptimizerInfo.h"

static Robot SimRobotObj;
static int SwingLinkInfoIndex;
static std::vector<int> SwingLinkChain;
static std::vector<double> ReferenceConfig;
static std::vector<double> CurrentConfig;
static SelfLinkGeoInfo SelfLinkGeoObj;
static double EndEffectorTol = 0.001;     // 1mm

struct TouchDownOptimization: public NonlinearOptimizerInfo
{
  TouchDownOptimization():NonlinearOptimizerInfo(){};

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
      {
        x_vec[i] = x[i];
      }
      std::vector<double> F_val = TouchDownOptimizationNCons(*n, *neF, x_vec);
      for (int i = 0; i < *neF; i++)
      {
        F[i] = F_val[i];
      }
    }
  void Solve(std::vector<double> &RobotConfig)
  {
    int StartType = 0;
    NonlinearProb.solve(StartType, neF, n, ObjAdd, ObjRow, ObjNConstraint,
      xlow, xupp, Flow, Fupp,
      x, xstate, xmul, F, Fstate, Fmul,
      nS, nInf, sumInf);
      for (int i = 0; i < n; i++)
      {
        RobotConfig[i] = x[i];
      }
      delete []x;      delete []xlow;   delete []xupp;
      delete []xmul;   delete []xstate;

      delete []F;      delete []Flow;   delete []Fupp;
      delete []Fmul;   delete []Fstate;
  }
  static std::vector<double> TouchDownOptimizationNCons(const int & nVar, const int & nObjNCons, const std::vector<double> & ActiveConfigOpt)
  {
    // This funciton provides the constraint for the configuration variable
    std::vector<double> F(nObjNCons);
    double ConfigDiff = 0.0;
    for (int i = 0; i < SwingLinkChain.size(); i++)
    {
      ReferenceConfig[SwingLinkChain[i]] = ActiveConfigOpt[i];
      double ConfigDiff_i = CurrentConfig[SwingLinkChain[i]] - ActiveConfigOpt[i];
      ConfigDiff+=ConfigDiff_i*ConfigDiff_i;
    }
    SimRobotObj.UpdateConfig(Config(ReferenceConfig));
    SimRobotObj.UpdateGeometry();
    F[0] = ConfigDiff;

    int ConstraintIndex = 1;
    // Self-collision constraint
    std::vector<double> SelfCollisionDistVec(SwingLinkChain.size()-3);
    for (int i = 0; i < SwingLinkChain.size()-3; i++)     // Due to the bounding box size of torso link
    {
      Box3D Box3DObj = SimRobotObj.geometry[SwingLinkChain[i]]->GetBB();
      std::vector<Vector3> BoxVerticesVec = BoxVertices(Box3DObj);
      std::vector<double> DistVec(BoxVerticesVec.size());
      for (int j = 0; j < BoxVerticesVec.size(); j++)
        DistVec[j] = SelfLinkGeoObj.SelfCollisionDist(SwingLinkInfoIndex, BoxVerticesVec[j]);

      SelfCollisionDistVec[i] = *std::min_element(DistVec.begin(), DistVec.end());
    }

    F[ConstraintIndex] = *std::min_element(SelfCollisionDistVec.begin(), SelfCollisionDistVec.end());
    ConstraintIndex+=1;

    for (int i = 0; i < NonlinearOptimizerInfo::RobotLinkInfo[SwingLinkInfoIndex].LocalContacts.size(); i++)
    {
      Vector3 LinkiPjPos;
      SimRobotObj.GetWorldPosition(NonlinearOptimizerInfo::RobotLinkInfo[SwingLinkInfoIndex].LocalContacts[i], NonlinearOptimizerInfo::RobotLinkInfo[SwingLinkInfoIndex].LinkIndex, LinkiPjPos);
      F[ConstraintIndex] = SDFInfo.SignedDistance(LinkiPjPos);
      ConstraintIndex+=1;
    }
    return F;
  }
};

std::vector<double> TouchDownConfigOptimization(const Robot & SimRobot, ReachabilityMap & RMObject, SelfLinkGeoInfo & _SelfLinkGeoObj, const std::vector<double> & RawqDes, const int & _SwingLinkInfoIndex){
  SimRobotObj = SimRobot;
  SwingLinkInfoIndex = _SwingLinkInfoIndex;
  SwingLinkChain = RMObject.EndEffectorLink2Pivotal[SwingLinkInfoIndex];
  ReferenceConfig = SimRobot.q;         // Used for internal optimization
  CurrentConfig = SimRobot.q;           // Used for configuration comparison
  SelfLinkGeoObj = _SelfLinkGeoObj;

  TouchDownOptimization TouchDownOptimizationProblem;

  // Static Variable Substitution
  int n = SwingLinkChain.size();
  std::vector<double> SwingLinkConfigs(n);

  int neF = 1;
  neF += 1;                                                                                           // Self-Collision
  neF = neF + NonlinearOptimizerInfo::RobotLinkInfo[SwingLinkInfoIndex].LocalContacts.size();         // Touch-down constraint
  TouchDownOptimizationProblem.InnerVariableInitialize(n, neF);

  /*
    Initialize the bounds of variables
  */
  std::vector<double> xlow_vec(n), xupp_vec(n);
  for (int i = 0; i < n; i++){
    // Configuration
    xlow_vec[i] = SimRobot.qMin(SwingLinkChain[i]);
    xupp_vec[i] = SimRobot.qMax(SwingLinkChain[i]);
    SwingLinkConfigs[i] = ReferenceConfig[SwingLinkChain[i]];
  }
  TouchDownOptimizationProblem.VariableBoundsUpdate(xlow_vec, xupp_vec);

  /*
    Initialize the bounds of variables
  */
  std::vector<double> Flow_vec(neF), Fupp_vec(neF);
  for (int i = 0; i < neF; i++) {
    Flow_vec[i] = 0;
    Fupp_vec[i] = 1e10;
  }
  for (int i = neF - NonlinearOptimizerInfo::RobotLinkInfo[SwingLinkInfoIndex].LocalContacts.size(); i < neF; i++) {
    Flow_vec[i] = 0;
    Fupp_vec[i] = EndEffectorTol;
  }

  TouchDownOptimizationProblem.ConstraintBoundsUpdate(Flow_vec, Fupp_vec);

  /*
    Initialize the seed guess
  */
  TouchDownOptimizationProblem.SeedGuessUpdate(SwingLinkConfigs);

  /*
    Given a name of this problem for the output
  */
  TouchDownOptimizationProblem.ProblemNameUpdate("TouchDownOptimizationProblem", 0);

  // Here we would like allow much more time to be spent on IK
  TouchDownOptimizationProblem.NonlinearProb.setIntParameter("Iterations limit", 5000);
  TouchDownOptimizationProblem.NonlinearProb.setIntParameter("Major iterations limit", 250);
  TouchDownOptimizationProblem.NonlinearProb.setIntParameter("Major print level", 0);
  TouchDownOptimizationProblem.NonlinearProb.setIntParameter("Minor print level", 0);
  /*
    ProblemOptions seting
  */
  // Solve with Finite-Difference
  TouchDownOptimizationProblem.ProblemOptionsUpdate(0, 3);
  TouchDownOptimizationProblem.Solve(SwingLinkConfigs);

  std::vector<double> OptConfig = ReferenceConfig;

  for (int i = 0; i < n; i++)
  {
    OptConfig[SwingLinkChain[i]] = SwingLinkConfigs[i];
  }
  SimRobotObj.UpdateConfig(Config(OptConfig));
  SimRobotObj.UpdateGeometry();

  std::string ConfigPath = "/home/motion/Desktop/Whole-Body-Planning-for-Push-Recovery/build/";
  string _OptConfigFile = "TouchDownConfig.config";
  RobotConfigWriter(OptConfig, ConfigPath, _OptConfigFile);

  // Self-collision constraint numerical checker
  std::vector<double> SelfCollisionDistVec(SwingLinkChain.size()-3);
  for (int i = 0; i < SwingLinkChain.size()-3; i++)     // Due to the bounding box size of torso link
  {
    Box3D Box3DObj = SimRobotObj.geometry[SwingLinkChain[i]]->GetBB();
    std::vector<Vector3> BoxVerticesVec = BoxVertices(Box3DObj);
    std::vector<double> DistVec(BoxVerticesVec.size());
    for (int j = 0; j < BoxVerticesVec.size(); j++)
    {
      DistVec[j] = SelfLinkGeoObj.SelfCollisionDist(SwingLinkInfoIndex, BoxVerticesVec[j]);
    }
    SelfCollisionDistVec[i] = *std::min_element(DistVec.begin(), DistVec.end());
  }
  double SelfCollisionDistTol = *std::min_element(SelfCollisionDistVec.begin(), SelfCollisionDistVec.end());

  if(SelfCollisionDistTol<-0.0025){
      std::printf("Touch Down Optimization Failure due to Self-collision for Link %d! \n", NonlinearOptimizerInfo::RobotLinkInfo[SwingLinkInfoIndex].LinkIndex);
      return RawqDes;
  }
  OptConfig = YPRShifter(OptConfig);
  return OptConfig;
}
