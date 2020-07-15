#include "CommonHeader.h"
#include "RobotInfo.h"
#include "NonlinearOptimizerInfo.h"

static double ReductionMagnitude(const double & cur_s, const double & phase_s, const double & reduction_ratio){
  if(cur_s<phase_s) return 1.0;
  else{
    if(cur_s<=1.0){
      return 1.0 - (1.0 - reduction_ratio)/(1.0 - phase_s) * (cur_s - phase_s);
    } else return reduction_ratio;
  }
}

ControlReferenceInfo TaskTrajectoryPlanning(Robot & SimRobotInner, const InvertedPendulumInfo & InvertedPendulumInner, ReachabilityMap & RMObject,SelfLinkGeoInfo & SelfLinkGeoObj,
                                            EndEffectorPathInfo & EndEffectorPathObj, SimPara & SimParaObj){

  // Here a task-space planning method is used to figure out robot's whole-body motion.
  InvertedPendulumInfo InvertedPendulumObj = InvertedPendulumInner;
  int SwingLinkInfoIndex = SimParaObj.getSwingLinkInfoIndex();
  Vector3 EndEffectorInitxDir, EndEffectorInityDir;   // Eventually these two directions should be orthgonal to goal direction.
  getEndEffectorXYAxes(SimRobotInner, SwingLinkInfoIndex, EndEffectorInitxDir, EndEffectorInityDir);
  Vector3 EndEffectorGoalDirection = SimParaObj.getGoalDirection();
  std::vector<int> SwingLinkChain = RMObject.EndEffectorLink2Pivotal[SwingLinkInfoIndex];

  // The main idea is that end effector will gradually move to be aligned with goal direction.
  double CurrentTime = 0.0;
  Config CurrentConfig = Config(SimRobotInner.q);
  Config CurrentVelocity = SimRobotInner.dq;
  Vector3 CurrentContactPos = SimParaObj.getContactInit();

  std::vector<double>   TimeTraj;
  std::vector<Config>   WholeBodyConfigTraj;
  std::vector<Config>   WholeBodyVelocityTraj;
  std::vector<Vector3>  PlannedEndEffectorTraj;

  TimeTraj.push_back(CurrentTime);
  WholeBodyConfigTraj.push_back(CurrentConfig);
  WholeBodyVelocityTraj.push_back(CurrentVelocity);
  PlannedEndEffectorTraj.push_back(CurrentContactPos);

  SimParaObj.setTrajConfigOptFlag(false);

  ControlReferenceInfo ControlReferenceObj;
  ControlReferenceObj.setReadyFlag(false);

  double  sVal = 0.0;
  double  sBoundary = SimParaObj.PhaseRatio;
  double  PhaseTimeStep = SimParaObj.PhaseTimeStep;
  bool    PenetrationFlag = false;

  Config PreVelocity = 0.0 * CurrentVelocity;           // Initialize to be zero
  std::vector<double> sTraj;
  sTraj.push_back(0.0);

  while ((sVal<1.0)&& (!PenetrationFlag)){
    std::vector<double> NextConfig, NextVelocity;
    double sNew;
    std::vector<double> CurrentBaseDelta = CurrentBaseDeltaCal(SimRobotInner, InvertedPendulumObj, PhaseTimeStep);
    double DampingRatio = ReductionMagnitude(sVal, SimParaObj.PhaseRatio, SimParaObj.ReductionRatio);
    TaskTrajectoryPlanningInner(      sVal, sNew,                       SimRobotInner,
                                      PreVelocity,                      CurrentBaseDelta,
                                      WholeBodyConfigTraj.back(),       WholeBodyVelocityTraj.back(),
                                      NextConfig,                       NextVelocity,
                                      EndEffectorPathObj,               SwingLinkChain,
                                      EndEffectorInitxDir,              EndEffectorInityDir,
                                      SimParaObj,
                                      PhaseTimeStep,                    DampingRatio);
    SimRobotInner.UpdateConfig(Config(NextConfig));
    Vector3 EndEffectorAvgPos;
    SimRobotInner.GetWorldPosition( NonlinearOptimizerInfo::RobotLinkInfo[SwingLinkInfoIndex].AvgLocalContact,
                                    NonlinearOptimizerInfo::RobotLinkInfo[SwingLinkInfoIndex].LinkIndex, EndEffectorAvgPos);
    sVal = EndEffectorPathObj.Pos2s(EndEffectorAvgPos);

    Config UpdatedConfig  = WholeBodyDynamicsIntegrator(SimRobotInner, InvertedPendulumObj, PhaseTimeStep);
    SimRobotInner.UpdateConfig(UpdatedConfig);

    //
    // Vector3 EndEffectorAvgPosPath;
    // EndEffectorPathObj.s2Pos(sVal, EndEffectorAvgPosPath);

    // sVal = sNew;
    // PenetrationFlag = PenetrationTester(SimRobotInner, SwingLinkInfoIndex);
    // if(PenetrationFlag)
    //   continue;

    CurrentTime+=PhaseTimeStep;
    CurrentConfig = UpdatedConfig;
    CurrentVelocity = Config(NextVelocity);

    PreVelocity = WholeBodyVelocityTraj.back();

    sTraj.push_back(sVal);
    TimeTraj.push_back(CurrentTime);
    WholeBodyConfigAppender(WholeBodyConfigTraj, CurrentConfig);
    WholeBodyVelocityTraj.push_back(Config(NextVelocity));
    PlannedEndEffectorTraj.push_back(CurrentContactPos);

  }
  std::cout<<"sTraj: "<<std::endl;
  for(int i = 0; i<sTraj.size(); i++){
    std::cout<<sTraj[i]<<std::endl;
  }
  std::cout<<"TimeTraj: "<<std::endl;
  for(int i = 0; i<TimeTraj.size(); i++){
    std::cout<<TimeTraj[i]<<std::endl;
  }
  for (int i = 0; i < WholeBodyConfigTraj.size(); i++) {
    std::string ConfigPath = "/home/motion/Desktop/Whole-Body-Planning-for-Push-Recovery/build/";
    std::string OptConfigFile = "InnerOpt" + std::to_string(i) + ".config";
    RobotConfigWriter(WholeBodyConfigTraj[i], ConfigPath, OptConfigFile);
  }
  for (int i = 0; i < WholeBodyConfigTraj.size(); i++) {
    std::string ConfigPath = "/home/motion/Desktop/Whole-Body-Planning-for-Push-Recovery/build/";
    std::string OptConfigFile = "InnerVel" + std::to_string(i) + ".config";
    RobotConfigWriter(WholeBodyVelocityTraj[i], ConfigPath, OptConfigFile);
  }
  std::cout<<"Config"<<std::endl;
  for (int i = 0; i < WholeBodyConfigTraj.size(); i++) {
    SwingLinkStatePrint(WholeBodyConfigTraj[i], SwingLinkChain);
  }
  std::cout<<"Velocity"<<std::endl;
  for (int i = 0; i < WholeBodyConfigTraj.size(); i++) {
    SwingLinkStatePrint(WholeBodyVelocityTraj[i], SwingLinkChain);
  }
  LinearPath WholeBodyConfigTrajPath(TimeTraj, WholeBodyConfigTraj);
  std::ofstream WholeBodyConfigTrajFile;
  const string  WholeBodyConfigTrajName = "WholeBodyConfigTraj.path";
  WholeBodyConfigTrajFile.open(WholeBodyConfigTrajName.c_str());
  WholeBodyConfigTrajPath.Save(WholeBodyConfigTrajFile);
  WholeBodyConfigTrajFile.close();

  if(SimParaObj.getTrajConfigOptFlag()){
    std::vector<ContactStatusInfo> GoalContactInfo = SimParaObj.FixedContactStatusInfo;
    for(int i = 0; i<GoalContactInfo[SwingLinkInfoIndex].LocalContactStatus.size(); i++)
      GoalContactInfo[SwingLinkInfoIndex].LocalContactStatus[i] = 1;
    ControlReferenceObj.SetInitContactStatus(SimParaObj.FixedContactStatusInfo);
    ControlReferenceObj.SetGoalContactStatus(GoalContactInfo);
    ControlReferenceObj.TrajectoryUpdate(TimeTraj, WholeBodyConfigTraj, PlannedEndEffectorTraj);
    double EstFailureMetric = EstimatedFailureMetric(SimRobotInner, GoalContactInfo, InvertedPendulumObj.COMPos, InvertedPendulumObj.COMVel);
    ControlReferenceObj.setFailueMetric(EstFailureMetric);
    ControlReferenceObj.setReadyFlag(true);
  }
  return ControlReferenceObj;
}
