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

static double PositivePosNVel(  const double & PosDiff,
                                const double & InitVelocity,  double & GoalVelocity,
                                const double & VelocityBound, const double & AccBound){
   // Two procedure => accelerate to velocity limit and remains there
   double AccDist = (VelocityBound * VelocityBound - InitVelocity * InitVelocity)/2.0;
   if(AccDist>PosDiff){ // Does not need to accelerate to limit to reach GoalPos
     double AccTime = (-InitVelocity + sqrt(InitVelocity * InitVelocity + 2.0 * AccBound * PosDiff))/(1.0 * AccBound);
     GoalVelocity = InitVelocity + AccTime * AccBound;
     return AccTime;
   }
   else { // This means that DOF accelerates to limit and then remain on maximum speed
     double AccTime1 = (VelocityBound - InitVelocity)/AccBound;
     double AccTime2 = (PosDiff - AccDist)/VelocityBound;
     double AccTime = AccTime1 + AccTime2;
     GoalVelocity = VelocityBound;
     return AccTime;
   }
 }

 static double NegativePosNVel(  const double & PosDiff,
                                 const double & InitVelocity,  double & GoalVelocity,
                                 const double & VelocityBound, const double & AccBound){
  double AccDist = (VelocityBound * VelocityBound - InitVelocity * InitVelocity)/2.0;
  if(AccDist>-PosDiff){
    double AccTime = (-InitVelocity + sqrt(InitVelocity * InitVelocity - 2.0 * AccBound * PosDiff))/(1.0 * AccBound);
    GoalVelocity = InitVelocity - AccTime * AccBound;
    return AccTime;
  }
  else { // This means that DOF accelerates to limit and then remain on maximum speed
    double AccTime1 = (VelocityBound + InitVelocity)/AccBound;
    double AccTime2 = -(PosDiff + AccDist)/VelocityBound;
    double AccTime = AccTime1 + AccTime2;
    GoalVelocity = -VelocityBound;
    return AccTime;
  }
}

static double AccPhaseTimeInner(const double & InitPos,       const double & GoalPos,
                                const double & InitVelocity,  double & GoalVelocity,
                                const double & VelocityBound, const double & AccBound){
  // The Velocity/Acceleartion bound is assumed to be bidirectional.
  double PosDiff = GoalPos - InitPos;
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

static void Velocity2Pos(   const double & InitPos,       const double & GoalPos,
                            const double & InitVelocity,  double & GoalVelocity,
                            const double & VelocityBound, const double & AccBound,
                            const double & DurationTime){
  // Figure out the velocity that this system reaches GoalPos with Duration time.
  double PosDiff = GoalPos - InitPos;
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
    else GoalVelocity = VelocityBound;
  }
}

static double AccPhaseTime( const std::vector<double> & CurConfig,      const std::vector<double> & NextConfig,
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
    double InitVelocity = CurVelocity[SwingLinkChain[i]];
    double GoalVelocity;
    double VelBound = VelocityBound[SwingLinkChain[i]];
    double AccBound = AccelerationBound[SwingLinkChain[i]];
    double AccPhaseTime = AccPhaseTimeInner(InitPos, GoalPos, InitVelocity, GoalVelocity, VelBound, AccBound);
    AccPhaseTimeTotal[i] = AccPhaseTime;
    SwingLinkChainVelocity[i] = GoalVelocity;
  }
  double AccTime = *max_element(AccPhaseTimeTotal.begin(), AccPhaseTimeTotal.end());
  double AccLinkIndex = std::distance(AccPhaseTimeTotal.begin(), std::max_element(AccPhaseTimeTotal.begin(), AccPhaseTimeTotal.end()));
  for (int i = 0; i < SwingLinkChain.size(); i++) {
    double InitPos = CurConfig[SwingLinkChain[i]];
    double GoalPos = NextConfig[SwingLinkChain[i]];
    double InitVelocity = CurVelocity[SwingLinkChain[i]];
    double GoalVelocity;
    double VelBound = VelocityBound[SwingLinkChain[i]];
    double AccBound = AccelerationBound[SwingLinkChain[i]];
    Velocity2Pos(InitPos, GoalPos, InitVelocity, GoalVelocity, VelBound, AccBound, AccTime);
    NextVelocity[SwingLinkChain[i]] = GoalVelocity;
  }
  return AccTime;
}

ControlReferenceInfo TrajectoryPlanning(Robot & SimRobotInner, ReachabilityMap & RMObject,SelfLinkGeoInfo & SelfLinkGeoObj,
                                              EndEffectorPathInfo & EndEffectorPathObj, SimPara & SimParaObj){



  RobotLink3D EndEffectorLink = SimRobotInner.links[NonlinearOptimizerInfo::RobotLinkInfo[SimParaObj.getSwingLinkInfoIndex()].LinkIndex];

  Vector3 EndEffectorInitDir;
  EndEffectorInitDir.x = EndEffectorLink.T_World.R.data[2][0];
  EndEffectorInitDir.y = EndEffectorLink.T_World.R.data[2][1];
  EndEffectorInitDir.z = EndEffectorLink.T_World.R.data[2][2];

  Vector3 EndEffectorGoalDir = SimParaObj.getGoalDirection();

  std::vector<int> SwingLinkChain = RMObject.EndEffectorLink2Pivotal[SimParaObj.getSwingLinkInfoIndex()];

  // The main idea is that end effector will gradually move to be aligned with goal direction.
  const int sNumber = 11;                 // 11 data points => 10 segments
  double sDiff = 1.0/(1.0 * sNumber - 1.0);
  double sVal = 0.0;
  std::vector<double> projVec = ProjectionLength(EndEffectorInitDir, EndEffectorGoalDir, sNumber - 1);

  double CurrentTime = 0.0;
  Config CurrentConfig = Config(YPRShifter(SimRobotInner.q));
  Vector3 CurrentContactPos = SimParaObj.getContactInit();

  std::vector<double>   TimeTraj;                 TimeTraj.reserve(sNumber);
  std::vector<Config>   WholeBodyConfigTraj;      WholeBodyConfigTraj.reserve(sNumber);
  std::vector<Config>   WholeBodyVelocityTraj;    WholeBodyVelocityTraj.reserve(sNumber);
  std::vector<Vector3>  EndEffectorTraj;          EndEffectorTraj.reserve(sNumber);

  TimeTraj.push_back(CurrentTime);
  WholeBodyConfigTraj.push_back(CurrentConfig);
  WholeBodyVelocityTraj.push_back(SimRobotInner.dq);
  EndEffectorTraj.push_back(CurrentContactPos);

  SimParaObj.setTrajConfigOptFlag(false);

  double sBoundary = round((sNumber - 1) * SimParaObj.PhaseRatio);
  double resS = sNumber - sBoundary;

  for (int sIndex = 1; sIndex < sNumber; sIndex++) {
    sVal = 1.0 * sIndex * sDiff;
    EndEffectorPathObj.s2Pos(sVal, CurrentContactPos);
    SimParaObj.setCurrentContactPos(CurrentContactPos);
    std::vector<double> NextConfig = TrajConfigOptimazation(SimRobotInner, RMObject, SelfLinkGeoObj, SimParaObj, projVec[sIndex], sIndex);
    std::vector<double> NextVelocity(NextConfig.size());

    if(!SimParaObj.getTrajConfigOptFlag()) break;
    else {
      if(sIndex<=sBoundary){
        double StageTime = AccPhaseTime(CurrentConfig, NextConfig,
                                        WholeBodyVelocityTraj[sIndex-1], NextVelocity,
                                        SimRobotInner.velMax, SimRobotInner.accMax,
                                        SwingLinkChain);
      } else {

      }
    }
  }
}

//   bool OptFlag = true;
//   while((sIndex<sNumber)&&(OptFlag == true))
//   {
//     sVal = 1.0 * sIndex * sDiff;
//     EndEffectorPathObj.s2Pos(sVal, CurrentContactPos);
//     switch (sIndex)
//     {
//       case 4:   LastFlag = true;
//       break;
//       default:  LastFlag = false;
//       break;
//     }
//     std::vector<double> OptConfig = TransientOptFn(SimRobotInner, SwingLimbIndex, SelfLinkGeoObj, CurrentContactPos, RMObject, OptFlag, LastFlag);;
//     if(OptFlag)
//     {
//       // Minimum Time Estimation.
//       double CurrentTime_i = MinimumTimeEstimation(SimRobotInner, RMObject.EndEffectorLink2Pivotal[SwingLimbIndex], CurrentConfig, Config(OptConfig));
//
//       CurrentTime+=CurrentTime_i;
//       TimeTraj.push_back(CurrentTime);
//       WholeBodyConfigTraj.push_back(Config(OptConfig));
//       EndEffectorTraj.push_back(CurrentContactPos);
//
//       // Then we should update the robot's CurrentConfig based on CurrentTime_i.
//       Config UpdatedConfig  = WholeBodyDynamicsIntegrator(SimRobotInner, OptConfig, PIPObj, InvertedPendulumObj, CurrentTime_i, sIndex);
//       SimRobotInner.UpdateConfig(UpdatedConfig);
//       CurrentConfig = UpdatedConfig;
//     }
//     else
//     {
//       break;
//     }
//     sIndex++;
//   }
// }
//
// switch (OptimalContact.size())
// {
//   case 0:
//   {
//     return ControlReferenceObj;
//   }
//   break;
//   default:
//   {
//     // Now too early to assume that FailureFlag is true.
//     bool FeasiFlag;
//     int OptimalContactIndex = 0;
//     while(OptimalContactIndex<OptimalContact.size())
//     {
//       Robot SimRobotInner = SimRobot;
//       Vector3 ContactGoal = OptimalContact[OptimalContactIndex];
//       Vector3 ContactGoalGrad = NonlinearOptimizerInfo::SDFInfo.SignedDistanceNormal(ContactGoal);
//       SplineObj = TransientTrajGene(SimRobotInner, SwingLimbIndex, SelfLinkGeoObj, RobotLinkInfo, ContactInit, ContactGoal, RMObject, DataRecorderObj, FeasiFlag);
//       if(FeasiFlag)
//       {
//         EndPathInfo EndEffectorPathObj(SplineObj, SwingLimbIndex);
//         InvertedPendulumInfo InvertedPendulumObj(PIPObj.theta, PIPObj.thetadot, COMPos, COMVel);
//
//         /*
//           1. At each sampled waypoints along the end effector trajectory, an end effector position is evaluated from path.
//           2. Based on robot's current configuration, an IK problem is solved to get robot's swing limb configuration.
//           3. A time-optimal executation duration is computed.
//           4. Based on that time, robot's whole-body configuration is updated with inverted pendulum model.
//           5. The whole algorithm terminates when robot's self-collision has been triggered or no feasible IK solution can be found.
//         */
//         const int sNumber = 5;                 // 6 sampled points will be extracted from EndEffectorPathObj.
//         int sIndex = 1;
//         double sDiff = 1.0/(1.0 * sNumber - 1.0);
//         double sVal = 0.0;
//         Config CurrentConfig = SimRobotInner.q;
//         CurrentConfig = YPRShifter(CurrentConfig);
//         double CurrentTime = 0.0;
//         Vector3 CurrentContactPos = ContactInit;
//
//         std::vector<double> TimeTraj;
//         std::vector<Config> WholeBodyConfigTraj;
//         std::vector<Vector3> EndEffectorTraj;
//
//         TimeTraj.reserve(sNumber);
//         WholeBodyConfigTraj.reserve(sNumber);
//         EndEffectorTraj.reserve(sNumber);
//
//         TimeTraj.push_back(CurrentTime);
//         WholeBodyConfigTraj.push_back(CurrentConfig);
//         EndEffectorTraj.push_back(CurrentContactPos);
//
//         bool OptFlag = true;
//         while((sIndex<sNumber)&&(OptFlag == true))
//         {
//           sVal = 1.0 * sIndex * sDiff;
//           EndEffectorPathObj.s2Pos(sVal, CurrentContactPos);
//           bool LastFlag;
//           switch (sIndex)
//           {
//             case 4:   LastFlag = true;
//             break;
//             default:  LastFlag = false;
//             break;
//           }
//           std::vector<double> OptConfig = TransientOptFn(SimRobotInner, SwingLimbIndex, SelfLinkGeoObj, CurrentContactPos, RMObject, OptFlag, LastFlag);;
//           if(OptFlag)
//           {
//             // Minimum Time Estimation.
//             double CurrentTime_i = MinimumTimeEstimation(SimRobotInner, RMObject.EndEffectorLink2Pivotal[SwingLimbIndex], CurrentConfig, Config(OptConfig));
//
//             CurrentTime+=CurrentTime_i;
//             TimeTraj.push_back(CurrentTime);
//             WholeBodyConfigTraj.push_back(Config(OptConfig));
//             EndEffectorTraj.push_back(CurrentContactPos);
//
//             // Then we should update the robot's CurrentConfig based on CurrentTime_i.
//             Config UpdatedConfig  = WholeBodyDynamicsIntegrator(SimRobotInner, OptConfig, PIPObj, InvertedPendulumObj, CurrentTime_i, sIndex);
//             SimRobotInner.UpdateConfig(UpdatedConfig);
//             CurrentConfig = UpdatedConfig;
//           }
//           else
//           {
//             break;
//           }
//           sIndex++;
//         }
//         // Here the inner optimiztion loop has been finished!
//         if(OptFlag)
//         {
//           ControlReferenceObj.TrajectoryUpdate(TimeTraj, WholeBodyConfigTraj, EndEffectorTraj, ContactGoal, ContactGoalGrad, EndEffectorPathObj.TotalLength, Type);
//           std::vector<ContactStatusInfo> GoalContactInfo = FixedRobotContactInfo;
//           for(int i = 0; i<FixedRobotContactInfo[SwingLimbIndex].LocalContactStatus.size(); i++)
//           {
//             GoalContactInfo[SwingLimbIndex].LocalContactStatus[i] = 1;
//           }
//           ControlReferenceObj.InitContactInfo = FixedRobotContactInfo;
//           ControlReferenceObj.GoalContactInfo = GoalContactInfo;
//
//           DataRecorderObj.OptConfigs = WholeBodyConfigTraj;
//
//           return ControlReferenceObj;
//         }
//       }
//       OptimalContactIndex++;
//     }
//   }
//   break;
// }
