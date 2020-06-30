#include "CommonHeader.h"
#include "NonlinearOptimizerInfo.h"

static std::vector<double> ProjectionLength(const Vector3 & InitDir, const Vector3 & GoalDir, const int & gridNo){
  double x = GoalDir.dot(InitDir);
  double y = (InitDir - x * GoalDir).length();
  double angle_rad = atan2 (y,x);
  double angle_unit = angle_rad/(1.0 * gridNo);
  std::vector<double> projVec(gridNo);
  for (int i = 0; i < gridNo; i++) {
    angle_rad-=angle_unit;
    projVec[i] = cos(angle_rad);
  }
  return projVec;
}

static std::vector<double> ReductionRatioVecGene(const double & ReductionRatio, const int & sRes){
  // A linear reduction is assumed at velocity bound.
  double ReductionUnit = (1.0 - ReductionRatio)/(1.0 * sRes);
  std::vector<double> ReductionRatioVec(sRes);
  double ReductionRatio_i = 1.0;
  for (int i = 0; i < sRes; i++){
    ReductionRatio_i-=ReductionUnit;
    ReductionRatioVec[i] = ReductionRatio_i;
  }
  return ReductionRatioVec;
}

static double PositivePosNVel(  double PosDiff,
                                const double & RawInitVelocity,  double & GoalVelocity,
                                const double & VelocityBound, const double & AccBound){
   // Two procedure => accelerate to velocity limit and remains there
   double InitVelocity = RawInitVelocity;
   double AccDist = (VelocityBound * VelocityBound - InitVelocity * InitVelocity)/(2.0 * AccBound);
   double AccTimeOffset = 0.0;
   if(AccDist<0.0){
     // This happens in velocity reduction stage where the initial velocity is larger than velocity bound.
     // As a result, the velocity needs to be dampened to be within boundary
     double AccDistOffset = -AccDist;
     PosDiff-=AccDistOffset;
     AccTimeOffset = (InitVelocity - VelocityBound)/AccBound;
     InitVelocity = VelocityBound;
     AccDist = 0.0;
   }
   if(AccDist>PosDiff){ // Does not need to accelerate to limit to reach GoalPos
     double AccTime = (-InitVelocity + sqrt(InitVelocity * InitVelocity + 2.0 * AccBound * PosDiff))/(1.0 * AccBound);
     GoalVelocity = InitVelocity + AccTime * AccBound;
     return AccTime + AccTimeOffset;
   }
   else { // This means that DOF accelerates to limit and then remain on maximum speed
     double AccTime1 = (VelocityBound - InitVelocity)/AccBound;
     double AccTime2 = (PosDiff - AccDist)/VelocityBound;
     double AccTime = AccTime1 + AccTime2;
     GoalVelocity = VelocityBound;
     return AccTime + AccTimeOffset;
   }
 }

 static double NegativePosNVel(  double PosDiff,
                                 const double & RawInitVelocity,  double & GoalVelocity,
                                 const double & VelocityBound, const double & AccBound){
  double InitVelocity = RawInitVelocity;
  double AccDist = (VelocityBound * VelocityBound - InitVelocity * InitVelocity)/(2.0 * AccBound);
  double AccTimeOffset = 0.0;
  if(AccDist<0.0){
    double AccDistOffset = -AccDist;
    PosDiff+=AccDistOffset;
    AccTimeOffset = -(InitVelocity - VelocityBound)/AccBound;
    InitVelocity = -VelocityBound;
    AccDist = 0.0;
  }
  if(AccDist>-PosDiff){
    double AccTime = (InitVelocity + sqrt(InitVelocity * InitVelocity - 2.0 * AccBound * PosDiff))/(1.0 * AccBound);
    GoalVelocity = InitVelocity - AccTime * AccBound;
    return AccTime + AccTimeOffset;
  }
  else { // This means that DOF accelerates to limit and then remain on maximum speed
    double AccTime1 = (VelocityBound + InitVelocity)/AccBound;
    double AccTime2 = -(PosDiff + AccDist)/VelocityBound;
    double AccTime = AccTime1 + AccTime2;
    GoalVelocity = -VelocityBound;
    return AccTime + AccTimeOffset;
  }
}

static double AccPhaseTimeInner(const double & PosDiff,
                                const double & InitVelocity,  double & GoalVelocity,
                                const double & VelocityBound, const double & AccBound){
  // The Velocity/Acceleartion bound is assumed to be bidirectional.
  if(PosDiff>0.0){
    if(InitVelocity>0.0){
      return PositivePosNVel(PosDiff, InitVelocity, GoalVelocity, VelocityBound, AccBound);
    }
    else {
      double AccTime1 = -InitVelocity/AccBound;
      double AccDist1 = InitVelocity * InitVelocity/(2.0 * AccBound);
      double AccTime2 = PositivePosNVel(PosDiff + AccDist1, 0.0, GoalVelocity, VelocityBound, AccBound);
      return AccTime1 + AccTime2;
    }
  }
  else{
    if(InitVelocity<0.0){
      return NegativePosNVel(PosDiff, InitVelocity, GoalVelocity, VelocityBound, AccBound);
    } else {
      double AccTime1 = InitVelocity/AccBound;
      double AccDist1 = InitVelocity * InitVelocity/(2.0 * AccBound);
      double AccTime2 = NegativePosNVel(PosDiff - AccDist1, 0.0, GoalVelocity, VelocityBound, AccBound);
      return AccTime1 + AccTime2;
    }
  }
}

static void Velocity2Pos(   const double & PosDiff,
                            const double & InitVelocity,  double & GoalVelocity,
                            const double & VelocityBound, const double & AccBound,
                            const double & DurationTime){
  // Figure out the velocity that this system reaches GoalPos with Duration time.
  if(PosDiff>0.0){
    double AccEst = 2.0 * (PosDiff - InitVelocity * DurationTime)/(DurationTime * DurationTime);  // This value is sure to be less than the bound.
    double VelocityEst = InitVelocity + AccEst * DurationTime;
    if(VelocityEst<VelocityBound) GoalVelocity = VelocityEst;
    else GoalVelocity = VelocityBound;
  }
  else{
    double AccEst = -2.0 * (PosDiff - InitVelocity * DurationTime)/(DurationTime * DurationTime);  // This value is sure to be less than the bound.
    double VelocityEst = InitVelocity - AccEst * DurationTime;
    if(VelocityEst>-VelocityBound) GoalVelocity = VelocityEst;
    else GoalVelocity = -VelocityBound;
  }
}

double AccPhaseTimePathMethod(  const std::vector<double> & CurConfig,      const std::vector<double> & NextConfig,
                                const std::vector<double> & CurVelocity,    std::vector<double> & NextVelocity,
                                const std::vector<double> & VelocityBound,  const std::vector<double> & AccelerationBound,
                                const std::vector<int> SwingLinkChain){
                       // This function solves for the time in acceleration phase.
  /*
    For each DOF, the following compuation procedure is conducted to figure the time.
    1. Determine whether should accelerate in postive or negative direction.
    2. Determine
  */
  std::vector<double> AccPhaseTimeTotal(SwingLinkChain.size());
  std::vector<double> SwingLinkChainVelocity(SwingLinkChain.size());
  for (int i = 0; i < SwingLinkChain.size(); i++) {
    double InitPos = CurConfig[SwingLinkChain[i]];
    double GoalPos = NextConfig[SwingLinkChain[i]];
    double PosDiff = GoalPos - InitPos;
    double InitVelocity = CurVelocity[SwingLinkChain[i]];
    double GoalVelocity;
    double VelBound = VelocityBound[SwingLinkChain[i]];
    double AccBound = AccelerationBound[SwingLinkChain[i]];
    double AccPhaseTime = AccPhaseTimeInner(PosDiff, InitVelocity, GoalVelocity, VelBound, AccBound);
    // printf("Link: %d PosDiff: %f,   InitVelocity: %f,   GoalVelocity: %f,   AccTime: %f and Valid: %d\n",
    //                               SwingLinkChain[i], PosDiff, InitVelocity, GoalVelocity, AccPhaseTime, GoalVelocity * PosDiff>0.0);
    AccPhaseTimeTotal[i] = AccPhaseTime;
    SwingLinkChainVelocity[i] = GoalVelocity;
  }
  double AccTime = *max_element(AccPhaseTimeTotal.begin(), AccPhaseTimeTotal.end());
  double AccLinkIndex = std::distance(AccPhaseTimeTotal.begin(), std::max_element(AccPhaseTimeTotal.begin(), AccPhaseTimeTotal.end()));
  for (int i = 0; i < SwingLinkChain.size(); i++) {
    double InitPos = CurConfig[SwingLinkChain[i]];
    double GoalPos = NextConfig[SwingLinkChain[i]];
    double PosDiff = GoalPos - InitPos;
    double InitVelocity = CurVelocity[SwingLinkChain[i]];
    double GoalVelocity;
    double VelBound = VelocityBound[SwingLinkChain[i]];
    double AccBound = AccelerationBound[SwingLinkChain[i]];
    Velocity2Pos(PosDiff, InitVelocity, GoalVelocity, VelBound, AccBound, AccTime);
    // printf("Updated Link: %d,   PosDiff: %f,  InitVelocity: %f,   GoalVelocity: %f and Valid: %d\n",
    //         SwingLinkChain[i], PosDiff, InitVelocity, GoalVelocity,  GoalVelocity * PosDiff>0.0);
    NextVelocity[SwingLinkChain[i]] = GoalVelocity;
  }
  return AccTime;
}

static double DecPhaseTimePathMethod(  const std::vector<double> & CurConfig,      const std::vector<double> & NextConfig,
                                const std::vector<double> & CurVelocity,    std::vector<double> & NextVelocity,
                                const std::vector<double> & VelocityBound,  const std::vector<double> & AccelerationBound,
                                const std::vector<int> SwingLinkChain, const double & ReductionRatio){
   // This function solves for the time in deceleration phase.
   std::vector<double> DampingVelocityBound(VelocityBound.size());
   for (int i = 0; i < DampingVelocityBound.size(); i++)
     DampingVelocityBound[i] = VelocityBound[i] * ReductionRatio;
   return AccPhaseTimePathMethod(CurConfig, NextConfig, CurVelocity, NextVelocity, DampingVelocityBound, AccelerationBound, SwingLinkChain);
}

ControlReferenceInfo TrajectoryPlanning(Robot & SimRobotInner, const InvertedPendulumInfo & InvertedPendulumInner, ReachabilityMap & RMObject,SelfLinkGeoInfo & SelfLinkGeoObj,
                                        EndEffectorPathInfo & EndEffectorPathObj, SimPara & SimParaObj){


  InvertedPendulumInfo InvertedPendulumObj = InvertedPendulumInner;
  RobotLink3D EndEffectorLink = SimRobotInner.links[NonlinearOptimizerInfo::RobotLinkInfo[SimParaObj.getSwingLinkInfoIndex()].LinkIndex];

  Vector3 EndEffectorInitDir;
  EndEffectorInitDir.x = EndEffectorLink.T_World.R.data[2][0];
  EndEffectorInitDir.y = EndEffectorLink.T_World.R.data[2][1];
  EndEffectorInitDir.z = EndEffectorLink.T_World.R.data[2][2];

  Vector3 EndEffectorGoalDir = SimParaObj.getGoalDirection();

  int SwingLinkInfoIndex = SimParaObj.getSwingLinkInfoIndex();
  std::vector<int> SwingLinkChain = RMObject.EndEffectorLink2Pivotal[SwingLinkInfoIndex];

  // The main idea is that end effector will gradually move to be aligned with goal direction.
  const int sNumber = 5;                 // (sNumber-1) segments
  double sDiff = 1.0/(1.0 * sNumber - 1.0);
  double sVal = 0.0;
  std::vector<double> projVec = ProjectionLength(EndEffectorInitDir, EndEffectorGoalDir, sNumber - 1);

  double CurrentTime = 0.0;
  Config CurrentConfig = Config(YPRShifter(SimRobotInner.q));
  Config CurrentVelocity = SimRobotInner.dq;
  Vector3 CurrentContactPos = SimParaObj.getContactInit();

  std::vector<double>   TimeTraj;                 TimeTraj.reserve(sNumber);
  std::vector<Config>   WholeBodyConfigTraj;      WholeBodyConfigTraj.reserve(sNumber);
  std::vector<Config>   WholeBodyVelocityTraj;    WholeBodyVelocityTraj.reserve(sNumber);
  std::vector<Vector3>  PlannedEndEffectorTraj;   PlannedEndEffectorTraj.reserve(sNumber);

  TimeTraj.push_back(CurrentTime);
  WholeBodyConfigTraj.push_back(CurrentConfig);
  WholeBodyVelocityTraj.push_back(CurrentVelocity);
  PlannedEndEffectorTraj.push_back(CurrentContactPos);

  SimParaObj.setTrajConfigOptFlag(false);

  Vector3 EndEffectorPosOffset; EndEffectorPosOffset.setZero();
  double sBoundary = round((sNumber - 1) * SimParaObj.PhaseRatio);
  double resS = sNumber - sBoundary;
  std::vector<double> ReducRatioVec = ReductionRatioVecGene(SimParaObj.PhaseRatio, resS);

  ControlReferenceInfo ControlReferenceObj;
  ControlReferenceObj.setReadyFlag(false);
  bool LastStageFlag = false;
  for (int sIndex = 1; sIndex < sNumber; sIndex++) {
    sVal = 1.0 * sIndex * sDiff;
    EndEffectorPathObj.s2Pos(sVal, CurrentContactPos);
    Vector3 PlannedCurrentContactPos = CurrentContactPos;
    SimParaObj.setCurrentContactPos(PlannedCurrentContactPos);
    switch (sIndex) {
      case 4: LastStageFlag = true;
      break;
      default:
      break;
    }
    std::vector<double> NextConfig = TrajConfigOptimazation(SimRobotInner, RMObject,
                                                            SelfLinkGeoObj, SimParaObj, sIndex);
    std::vector<double> NextVelocity = WholeBodyVelocityTraj[WholeBodyVelocityTraj.size()-1];
    Config UpdatedConfig;
    double StageTime;
    if(!SimParaObj.getTrajConfigOptFlag()) break;
    else {
      if(sIndex<sBoundary){
        StageTime = AccPhaseTimePathMethod( CurrentConfig, NextConfig,
                                            WholeBodyVelocityTraj[sIndex-1], NextVelocity,
                                            SimRobotInner.velMax, SimRobotInner.accMax,
                                            SwingLinkChain);
      }
      else {
        int ReducRatioVecIndex = sIndex - sBoundary;
        StageTime = DecPhaseTimePathMethod( CurrentConfig, NextConfig,
                                            WholeBodyVelocityTraj[sIndex-1], NextVelocity,
                                            SimRobotInner.velMax, SimRobotInner.accMax,
                                            SwingLinkChain, ReducRatioVec[ReducRatioVecIndex]);
      }
      SimRobotInner.UpdateConfig(Config(NextConfig));
      UpdatedConfig  = WholeBodyDynamicsIntegrator(SimRobotInner, InvertedPendulumObj, StageTime, sIndex);
      SimRobotInner.UpdateConfig(UpdatedConfig);
      if(!LastStageFlag)
        UpdatedConfig = OrientationOptimazation(SimRobotInner, SwingLinkChain, SimParaObj, projVec[sIndex], sIndex);
      else
        UpdatedConfig = LastStageConfigOptimazation(SimRobotInner, RMObject, SelfLinkGeoObj, SimParaObj, sIndex);
      SimRobotInner.UpdateConfig(UpdatedConfig);

      std::string ConfigPath = "/home/motion/Desktop/Whole-Body-Planning-for-Push-Recovery/build/";
      std::string OptConfigFile = "InnerOpt" + std::to_string(sIndex) + ".config";
      RobotConfigWriter(UpdatedConfig, ConfigPath, OptConfigFile);

      Vector3 EndEffectorContactPos;
      SimRobotInner.GetWorldPosition( NonlinearOptimizerInfo::RobotLinkInfo[SwingLinkInfoIndex].AvgLocalContact,
                                      NonlinearOptimizerInfo::RobotLinkInfo[SwingLinkInfoIndex].LinkIndex,
                                      EndEffectorContactPos);
      EndEffectorPosOffset = CurrentContactPos - EndEffectorContactPos;
      CurrentTime+=StageTime;
      CurrentConfig = UpdatedConfig;
      CurrentVelocity = Config(NextVelocity);
      // printf("CurrentTime: %f, %d th's StageTime %f\n", CurrentTime, sIndex, StageTime);

      TimeTraj.push_back(CurrentTime);
      WholeBodyConfigTraj.push_back(CurrentConfig);
      WholeBodyVelocityTraj.push_back(Config(NextVelocity));
      PlannedEndEffectorTraj.push_back(PlannedCurrentContactPos);
    }
  }
  for (int i = 0; i < WholeBodyConfigTraj.size(); i++) {
    std::string ConfigPath = "/home/motion/Desktop/Whole-Body-Planning-for-Push-Recovery/build/";
    std::string OptConfigFile = "InnerOpt" + std::to_string(i) + ".config";
    RobotConfigWriter(WholeBodyConfigTraj[i], ConfigPath, OptConfigFile);
  }
  if(SimParaObj.getTrajConfigOptFlag()){
    std::vector<ContactStatusInfo> GoalContactInfo = SimParaObj.FixedContactStatusInfo;
    for(int i = 0; i<GoalContactInfo[SwingLinkInfoIndex].LocalContactStatus.size(); i++)
      GoalContactInfo[SwingLinkInfoIndex].LocalContactStatus[i] = 1;
    ControlReferenceObj.SetInitContactStatus(SimParaObj.FixedContactStatusInfo);
    ControlReferenceObj.SetGoalContactStatus(GoalContactInfo);
    ControlReferenceObj.TrajectoryUpdate(TimeTraj, WholeBodyConfigTraj, PlannedEndEffectorTraj);
    ControlReferenceObj.setReadyFlag(true);
  }
  return ControlReferenceObj;
}
