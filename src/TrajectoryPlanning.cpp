#include "CommonHeader.h"
#include "RobotInfo.h"
#include "NonlinearOptimizerInfo.h"
#include "gurobi_c++.h"

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

static std::vector<double> EdgeProjVec(const Vector3 & InitxDir, const Vector3 & GoalDir, const int & gridNo){
  double xProj = InitxDir.dot(GoalDir);
  double xProjUnit = xProj/(1.0 * gridNo);
  std::vector<double> xProjVec(gridNo);
  double xProjVal = xProj - xProjUnit;
  for (int i = 0; i < gridNo; i++) {
    xProjVec[i] = xProjVal;
    xProjVal = xProjVal - xProjUnit;
  }
  return xProjVec;
}

static double EdgeProjMagnitude(const double & cur_s,  const Vector3 & InitxDir, const Vector3 & GoalDir){
  // This function calculates the projection based on current s'length.
  double xProj = InitxDir.dot(GoalDir);
  return (1.0 - cur_s) * xProj;
}

static double ReductionMagnitude(const double & cur_s, const double & phase_s, const double & reduction_ratio){
  if(cur_s<phase_s) return 1.0;
  else{
    if(cur_s<=1.0){
      return 1.0 - (1.0 - reduction_ratio)/(1.0 - phase_s) * (cur_s - phase_s);
    } else return reduction_ratio;
  }
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
    AccTimeOffset = -(InitVelocity + VelocityBound)/AccBound;
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

static void AccVioPosNVel2Time( double InitPos,       double & GoalPos,
                                double InitVelocity,  double & GoalVelocity,
                                double VelocityBound, double AccBound,
                                double DurationTime,  double AccEst){
  // Update GoalPos and GoalVelocity according to DurationTime when the violation is at acceleration level
  if(AccEst<=0){
    AccBound = -AccBound;
    VelocityBound = -VelocityBound;
  }
  double VelocityEst = InitVelocity + DurationTime * AccBound;
  if(VelocityEst * VelocityEst <=VelocityBound * VelocityBound){
    GoalPos = InitPos + (VelocityEst * VelocityEst - InitVelocity * InitVelocity)/(2.0 * AccBound);
    GoalVelocity = VelocityEst;
  }
  else{
      // Then the velocity first reaches maximum and then remains there for the rest of the time.
      double AccTime = (VelocityBound - InitVelocity)/AccBound;
      double RestTime = DurationTime - AccTime;
      double AccDist = (VelocityBound * VelocityBound - InitVelocity * InitVelocity)/(2.0 * AccBound);
      double FlatDist = GoalVelocity * RestTime;
      GoalPos = InitPos + AccDist + FlatDist;
      GoalVelocity = VelocityBound;
  }
  return;
}

static void VelVioPosNVel2Time( double InitPos,       double & GoalPos,
                                double InitVelocity,  double & GoalVelocity,
                                double VelocityBound, double AccBound,
                                double DurationTime,  double AccEst){
  // Update GoalPos and GoalVelocity according to DurationTime when the violation is at velocity level so AccEst should be within AccBound
  if(AccEst<=0){
    AccBound = -AccBound;
    VelocityBound = -VelocityBound; 
  }
  // First compare the time needed to accelerate to velocity bound
  double AccTime = (VelocityBound - InitVelocity)/AccBound;
  if(AccTime<=DurationTime){
    // First to maximum and then remain there.
    double AccDist = (VelocityBound * VelocityBound - InitVelocity * InitVelocity)/(2.0 * AccBound);
    double RestTime = DurationTime - AccTime;
    double RestDist = RestTime * VelocityBound;
    GoalPos = InitPos + AccDist + RestDist;
    GoalVelocity = VelocityBound;
  } else {
    // Just need to accelerate to DurationTime then
    GoalVelocity = InitVelocity + AccBound * DurationTime;
    GoalPos = InitPos + (GoalVelocity * GoalVelocity - InitVelocity * InitVelocity)/(2.0 * AccBound);
  }
}

static bool Velocity2Pos(     const double & InitPos,       double & GoalPos,
                              const double & InitVelocity,  double & GoalVelocity,
                              const double & VelocityBound, const double & AccBound,
                              const double & DurationTime, bool ModiFlag){
  bool FinderFlag = true;
  double AccTol = 1e-3;
  double PosDiff = GoalPos - InitPos;
  double AccEst = 2.0 * (PosDiff - InitVelocity * DurationTime)/(DurationTime * DurationTime);
  if(AccEst * AccEst<=AccBound * AccBound + AccTol){
    double VelocityEst = InitVelocity + AccEst * DurationTime;
    if(VelocityEst * VelocityEst<=VelocityBound * VelocityBound)
      GoalVelocity = VelocityEst;
    else{
      double MinDurationTime = AccPhaseTimeInner(PosDiff, InitVelocity, GoalVelocity, VelocityBound, AccBound);
      if(MinDurationTime<=DurationTime){
        if(VelocityEst>0) GoalVelocity = VelocityBound;
        else GoalVelocity = -VelocityBound;
      }
      else {
        if(ModiFlag)
        VelVioPosNVel2Time( InitPos, GoalPos, InitVelocity, GoalVelocity,
                            VelocityBound, AccBound, DurationTime, AccEst);
        FinderFlag = false;
      }
    }
  }
  else {
    // This indicates that the current time is too large or too small.
    if(ModiFlag)
    AccVioPosNVel2Time( InitPos, GoalPos, InitVelocity, GoalVelocity, 
                        VelocityBound, AccBound, DurationTime, AccEst);
    FinderFlag = false;
  }
  return FinderFlag;
}

static bool AccPhaseTimePathMethod( const std::vector<double> & CurConfig,          std::vector<double> & NextConfig,
                                    const std::vector<double> & CurVelocity,        std::vector<double> & NextVelocity,
                                    const std::vector<double> & VelocityBound,      const std::vector<double> & AccelerationBound,
                                    const std::vector<int> & SwingLinkChain,        double & AccTime){
  // This function solves for the time in acceleration phase.
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

  // std::string ConfigPath = "/home/motion/Desktop/Whole-Body-Planning-for-Push-Recovery/build/";
  // std::string OptConfigFile = "CurConfig.config";
  // RobotConfigWriter(CurConfig, ConfigPath, OptConfigFile);


  // ConfigPath = "/home/motion/Desktop/Whole-Body-Planning-for-Push-Recovery/build/";
  // OptConfigFile = "NextConfig.config";
  // RobotConfigWriter(NextConfig, ConfigPath, OptConfigFile);

  AccTime = *max_element(AccPhaseTimeTotal.begin(), AccPhaseTimeTotal.end());
  const int AccTimeGrid = 10;
  double AccTimeStep = AccTime/(1.0 * AccTimeGrid - 1.0);
  int AccTimeIndex = 0;
  std::vector<int> ReachableSizeVec(AccTimeGrid);
  while (AccTimeIndex<AccTimeGrid) {
    double AccTimeInner = AccTime - 1.0 * AccTimeIndex * AccTimeStep;
    int ReachableSize = 0;
    for (int i = 0; i < SwingLinkChain.size(); i++) {
      double InitPos = CurConfig[SwingLinkChain[i]];
      double GoalPos = NextConfig[SwingLinkChain[i]];
      double PosDiff = GoalPos - InitPos;
      double InitVelocity = CurVelocity[SwingLinkChain[i]];
      double GoalVelocity;
      double VelBound = VelocityBound[SwingLinkChain[i]];
      double AccBound = AccelerationBound[SwingLinkChain[i]];
      bool FinderFlag = Velocity2Pos(InitPos, GoalPos, InitVelocity, GoalVelocity, VelBound, AccBound, AccTimeInner, false);
      if(FinderFlag) ReachableSize++;
    }
    ReachableSizeVec[AccTimeIndex] = ReachableSize;
    AccTimeIndex++;
  }

  std::reverse(ReachableSizeVec.begin(), ReachableSizeVec.end());
  int ReachableSizeIndex = std::distance(ReachableSizeVec.begin(), std::max_element(ReachableSizeVec.begin(), ReachableSizeVec.end()));
  AccTime = AccTime - AccTimeStep * (1.0 * AccTimeGrid - 1.0 - 1.0 * ReachableSizeIndex);
  printf("AccTime: %f\n", AccTime);
  for (int i = 0; i < SwingLinkChain.size(); i++) {
    double InitPos = CurConfig[SwingLinkChain[i]];
    double GoalPos = NextConfig[SwingLinkChain[i]];
    double PosDiff = GoalPos - InitPos;
    double InitVelocity = CurVelocity[SwingLinkChain[i]];
    double GoalVelocity;
    double VelBound = VelocityBound[SwingLinkChain[i]];
    double AccBound = AccelerationBound[SwingLinkChain[i]];
    Velocity2Pos(InitPos, GoalPos, InitVelocity, GoalVelocity, VelBound, AccBound, AccTime, true);
    NextConfig[SwingLinkChain[i]] = GoalPos;
    NextVelocity[SwingLinkChain[i]] = GoalVelocity;
  }
  return true;
}

static bool DecPhaseTimePathMethod(   const std::vector<double> & CurConfig,      std::vector<double> & NextConfig,
                                        const std::vector<double> & CurVelocity,    std::vector<double> & NextVelocity,
                                        const std::vector<double> & VelocityBound,  const std::vector<double> & AccelerationBound,
                                        const std::vector<int> & SwingLinkChain,      const double & ReductionRatio, double & DecTime){
   // This function solves for the time in deceleration phase.
   std::vector<double> DampingVelocityBound(VelocityBound.size());
   for (int i = 0; i < DampingVelocityBound.size(); i++)
     DampingVelocityBound[i] = VelocityBound[i] * ReductionRatio;
   return AccPhaseTimePathMethod(CurConfig, NextConfig, CurVelocity, NextVelocity, DampingVelocityBound, AccelerationBound, SwingLinkChain, DecTime);
}

static bool StageStateOptimization( const double & sVal, const double & sUnit, double & sNew, Vector3 & CurrentContactPos,Robot & SimRobotInner,
                                    const Config & CurrentConfig, const Config & CurrentVelocity,
                                    std::vector<double> & NextConfig, std::vector<double> & NextVelocity,
                                    ReachabilityMap & RMObject, SelfLinkGeoInfo & SelfLinkGeoObj,
                                    EndEffectorPathInfo & EndEffectorPathObj, const std::vector<int> & SwingLinkChain,
                                    const Vector3 & EndEffectorInitxDir, const Vector3 & EndEffectorInityDir,
                                    SimPara & SimParaObj, double & StageTime){
  int TotalInter = 25;
  int IterInd = 0;
  Vector3 EndEffectorGoalDir = SimParaObj.getGoalDirection();
  double sCurUnit = sUnit;
  double scaleRatio = 1.25;
  while(IterInd<TotalInter){
    double sCur = sVal;
    sCur+=sCurUnit;
    sNew = sCur;
    EndEffectorPathObj.s2Pos(sCur, CurrentContactPos);
    Vector3 PlannedCurrentContactPos = CurrentContactPos;
    SimParaObj.setCurrentContactPos(PlannedCurrentContactPos);
    double EndEffectorProjx = EdgeProjMagnitude(sCur, EndEffectorInitxDir, EndEffectorGoalDir);
    double EndEffectorProjy = EdgeProjMagnitude(sCur, EndEffectorInityDir, EndEffectorGoalDir);
    NextConfig = TrajConfigOptimazation(SimRobotInner, RMObject, SelfLinkGeoObj, SimParaObj, EndEffectorProjx, EndEffectorProjy, sCur);
    NextVelocity = CurrentVelocity;
    Config UpdatedConfig;
    bool FinderFlag;
    if(!SimParaObj.getTrajConfigOptFlag()) return false;
    else
    {
      if(sCur<SimParaObj.PhaseRatio){
        FinderFlag = AccPhaseTimePathMethod( CurrentConfig, NextConfig,
                                            CurrentVelocity, NextVelocity,
                                            SimRobotInner.velMax, SimRobotInner.accMax,
                                            SwingLinkChain, StageTime);
      }
      else
      {
        double ReductionMag = ReductionMagnitude(sCur, SimParaObj.PhaseRatio, SimParaObj.ReductionRatio);
        FinderFlag = DecPhaseTimePathMethod( CurrentConfig, NextConfig,
                                            CurrentVelocity, NextVelocity,
                                            SimRobotInner.velMax, SimRobotInner.accMax,
                                            SwingLinkChain, ReductionMag, StageTime);
      }
      if(FinderFlag == true) return true;
    }
    sCurUnit = scaleRatio * sCurUnit;
  }
  return false;
}

ControlReferenceInfo TrajectoryPlanning(Robot & SimRobotInner, const InvertedPendulumInfo & InvertedPendulumInner, ReachabilityMap & RMObject,SelfLinkGeoInfo & SelfLinkGeoObj,
                                        EndEffectorPathInfo & EndEffectorPathObj, SimPara & SimParaObj){


  InvertedPendulumInfo InvertedPendulumObj = InvertedPendulumInner;
  int SwingLinkInfoIndex = SimParaObj.getSwingLinkInfoIndex();
  Vector3 EndEffectorInitxDir, EndEffectorInityDir;   // Eventually these two directions should be orthgonal to goal direction.
  getEndEffectorXYAxes(SimRobotInner, SwingLinkInfoIndex, EndEffectorInitxDir, EndEffectorInityDir);
  Vector3 EndEffectorGoalDir = SimParaObj.getGoalDirection();
  std::vector<int> SwingLinkChain = RMObject.EndEffectorLink2Pivotal[SwingLinkInfoIndex];

  // The main idea is that end effector will gradually move to be aligned with goal direction.
  double sVal = 0.0;
  double sUnit = 0.2;
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

  double sBoundary = SimParaObj.PhaseRatio;

  ControlReferenceInfo ControlReferenceObj;
  ControlReferenceObj.setReadyFlag(false);

  bool PenetrationFlag = false;
  int sIndex = 0;

  while ((sVal<1.0)&& (!PenetrationFlag)){
    std::vector<double> NextConfig, NextVelocity;
    double StageTime;
    double sNew;
    bool StageOptFlag = StageStateOptimization( sVal, sUnit, sNew, CurrentContactPos, SimRobotInner, CurrentConfig, WholeBodyVelocityTraj.back(),
                                                NextConfig, NextVelocity, RMObject, SelfLinkGeoObj, EndEffectorPathObj, SwingLinkChain,
                                                EndEffectorInitxDir, EndEffectorInityDir, SimParaObj, StageTime);
    sVal = sNew;
    if(!StageOptFlag) break;
    SimRobotInner.UpdateConfig(Config(NextConfig));
    Config UpdatedConfig  = WholeBodyDynamicsIntegrator(SimRobotInner, InvertedPendulumObj, StageTime);
    SimRobotInner.UpdateConfig(UpdatedConfig);

    if(sVal>=0.5){
      PenetrationFlag = PenetrationTester(SimRobotInner, SwingLinkInfoIndex);
      if(PenetrationFlag)
        UpdatedConfig = LastStageConfigOptimazation(SimRobotInner, RMObject, SelfLinkGeoObj, SimParaObj, -1);
    }

    CurrentTime+=StageTime;
    CurrentConfig = UpdatedConfig;
    CurrentVelocity = Config(NextVelocity);

    TimeTraj.push_back(CurrentTime);
    WholeBodyConfigAppender(WholeBodyConfigTraj, CurrentConfig);
    WholeBodyVelocityTraj.push_back(Config(NextVelocity));
    PlannedEndEffectorTraj.push_back(CurrentContactPos);

    SelfLinkGeoObj.LinkBBsUpdate(SimRobotInner);
    sIndex++;
  }

  // std::cout<<"TimeTraj: "<<std::endl;
  // for(int i = 0; i<TimeTraj.size(); i++){
  //   std::cout<<TimeTraj[i]<<std::endl;
  // }
  // for (int i = 0; i < WholeBodyConfigTraj.size(); i++) {
  //   std::string ConfigPath = "/home/motion/Desktop/Whole-Body-Planning-for-Push-Recovery/build/";
  //   std::string OptConfigFile = "InnerOpt" + std::to_string(i) + ".config";
  //   RobotConfigWriter(WholeBodyConfigTraj[i], ConfigPath, OptConfigFile);
  // }
  // for (int i = 0; i < WholeBodyConfigTraj.size(); i++) {
  //   std::string ConfigPath = "/home/motion/Desktop/Whole-Body-Planning-for-Push-Recovery/build/";
  //   std::string OptConfigFile = "InnerVel" + std::to_string(i) + ".config";
  //   RobotConfigWriter(WholeBodyVelocityTraj[i], ConfigPath, OptConfigFile);
  // }
  // std::cout<<"Config"<<std::endl;
  // for (int i = 0; i < WholeBodyConfigTraj.size(); i++) {
  //   SwingLinkStatePrint(WholeBodyConfigTraj[i], SwingLinkChain);
  // }
  // std::cout<<"Velocity"<<std::endl;
  // for (int i = 0; i < WholeBodyConfigTraj.size(); i++) {
  //   SwingLinkStatePrint(WholeBodyVelocityTraj[i], SwingLinkChain);
  // }
  // LinearPath WholeBodyConfigTrajPath(TimeTraj, WholeBodyConfigTraj);
  // std::ofstream WholeBodyConfigTrajFile;
  // const string  WholeBodyConfigTrajName = "WholeBodyConfigTraj.path";
  // WholeBodyConfigTrajFile.open(WholeBodyConfigTrajName.c_str());
  // WholeBodyConfigTrajPath.Save(WholeBodyConfigTrajFile);
  // WholeBodyConfigTrajFile.close();

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
