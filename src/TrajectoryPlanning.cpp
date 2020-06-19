#include "CommonHeader.h"
#include "NonlinearOptimizerInfo.h"

ControlReferenceInfo TrajectoryPlanning(const Robot & SimRobotObj, ReachabilityMap & RMObject,SelfLinkGeoInfo & SelfLinkGeoObj,
                                              EndEffectorPathInfo & EndEffectorPathObj, SimPara & SimParaObj){

  const int sNumber = 11;                 // 11 data points => 10 segments
  double sDiff = 1.0/(1.0 * sNumber - 1.0);
  double sVal = 0.0;

  Robot SimRobotInner = SimRobotObj;

  double CurrentTime = 0.0;
  // std::vector<double> CurrentConfigVec = SimRobotInner.q;
  Config CurrentConfig = Config(YPRShifter(SimRobotInner.q));
  Vector3 CurrentContactPos = SimParaObj.getContactInit();

  std::vector<double>   TimeTraj;                 TimeTraj.reserve(sNumber);
  std::vector<Config>   WholeBodyConfigTraj;      WholeBodyConfigTraj.reserve(sNumber);
  std::vector<Vector3>  EndEffectorTraj;          EndEffectorTraj.reserve(sNumber);

  TimeTraj.push_back(CurrentTime);
  WholeBodyConfigTraj.push_back(CurrentConfig);
  EndEffectorTraj.push_back(CurrentContactPos);

  for (int sIndex = 1; sIndex < sNumber; sIndex++) {
    sVal = 1.0 * sIndex * sDiff;
    EndEffectorPathObj.s2Pos(sVal, CurrentContactPos);
    bool LastStageFlag = false;
    if(sIndex == sNumber - 1) LastStageFlag = true;
    SimParaObj.setCurrentContactPos(CurrentContactPos);
    std::vector<double> OptConfig = TrajConfigOptimazation(SimRobotInner, RMObject, SelfLinkGeoObj, SimParaObj, LastStageFlag);
    if(!SimParaObj.getTrajConfigOptFlag()) break;
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
