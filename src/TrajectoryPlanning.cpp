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
  const int sNumber = 6;                 // 6 data points => 5 segments
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

  for (int sIndex = 1; sIndex < sNumber; sIndex++) {
    sVal = 1.0 * sIndex * sDiff;
    EndEffectorPathObj.s2Pos(sVal, CurrentContactPos);
    Vector3 PlannedCurrentContactPos = CurrentContactPos + EndEffectorPosOffset;
    SimParaObj.setCurrentContactPos(PlannedCurrentContactPos);
    std::vector<double> NextConfig = TrajConfigOptimazation(SimRobotInner, RMObject,
                                                            SelfLinkGeoObj, SimParaObj,
                                                            projVec[sIndex], sIndex);

    std::vector<double> NextVelocity = WholeBodyVelocityTraj[WholeBodyVelocityTraj.size()-1];
    Config UpdatedConfig;
    for (int i = 0; i < SwingLinkChain.size(); i++) {
      printf("Link %d's Velocity is %f\n", SwingLinkChain[i], NextVelocity[SwingLinkChain[i]]);
    }

    double StageTime;
    if(!SimParaObj.getTrajConfigOptFlag()) break;
    else {
      if(sIndex<sBoundary)
        StageTime = AccPhaseTimePathMethod( CurrentConfig, NextConfig,
                                            WholeBodyVelocityTraj[sIndex-1], NextVelocity,
                                            SimRobotInner.velMax, SimRobotInner.accMax,
                                            SwingLinkChain);
      else {


      }

      printf("%d th's StageTime %f\n", sIndex, StageTime);
      SimRobotInner.UpdateConfig(Config(NextConfig));
      UpdatedConfig  = WholeBodyDynamicsIntegrator(SimRobotInner, InvertedPendulumObj, StageTime, sIndex);
      SimRobotInner.UpdateConfig(UpdatedConfig);

      Vector3 EndEffectorContactPos;
      SimRobotInner.GetWorldPosition( NonlinearOptimizerInfo::RobotLinkInfo[SwingLinkInfoIndex].AvgLocalContact,
                                      NonlinearOptimizerInfo::RobotLinkInfo[SwingLinkInfoIndex].LinkIndex,
                                      EndEffectorContactPos);
      EndEffectorPosOffset = CurrentContactPos - EndEffectorContactPos;

      CurrentTime+=StageTime;
      CurrentConfig = UpdatedConfig;
      CurrentVelocity = Config(NextVelocity);

      TimeTraj.push_back(CurrentTime);
      WholeBodyConfigTraj.push_back(CurrentConfig);
      WholeBodyVelocityTraj.push_back(Config(NextVelocity));
      PlannedEndEffectorTraj.push_back(PlannedCurrentContactPos);
    }
  }
}

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
//         std::vector<Vector3> PlannedEndEffectorTraj;
//
//         TimeTraj.reserve(sNumber);
//         WholeBodyConfigTraj.reserve(sNumber);
//         PlannedEndEffectorTraj.reserve(sNumber);
//
//         TimeTraj.push_back(CurrentTime);
//         WholeBodyConfigTraj.push_back(CurrentConfig);
//         PlannedEndEffectorTraj.push_back(CurrentContactPos);
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
//             PlannedEndEffectorTraj.push_back(CurrentContactPos);
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
//           ControlReferenceObj.TrajectoryUpdate(TimeTraj, WholeBodyConfigTraj, PlannedEndEffectorTraj, ContactGoal, ContactGoalGrad, EndEffectorPathObj.TotalLength, Type);
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
