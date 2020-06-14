#include "CommonHeader.h"
#include "NonlinearOptimizerInfo.h"
#include <ctime>

static AllContactStatus ContactStatusExplorer(Robot & SimRobot, const std::vector<ContactStatusInfo> & RobotContactInfo){
  // First part is for contact modification.
  // Second part is for contact addition.
  Vector3 COMPos(0.0, 0.0, 0.0), COMVel(0.0, 0.0, 0.0);
  getCentroidalState(SimRobot, COMPos, COMVel);

  AllContactStatus AllContactStatusObj;
  int LegActNo = RobotContactInfo[0].LocalContactStatus[0] + RobotContactInfo[1].LocalContactStatus[0];
  int AllActNo = LegActNo + RobotContactInfo[2].LocalContactStatus[0] + RobotContactInfo[3].LocalContactStatus[0];

  // Contact Modification
  switch (AllActNo)
  {
    case 0:
    {
      std::cerr<<"No Active Contact!"<<endl;
      exit(1);
    }
    break;
    case 1:
      // No contact modification
    break;
    default:
    {
      // Algorithms for contact modification
      switch (LegActNo)
      {
        case 0:
        {
          std::cerr<<"No Active Foot Contact!"<<endl;
          exit(1);
        }
        break;
        case 1:
        {
          // In this case, the contact modification can only be conducted for hand contact if there exists.
          switch (RobotContactInfo[2].LocalContactStatus[0])
          {
            case 1:
            {
              std::vector<ContactStatusInfo> RobotContactInfoModi = RobotContactInfo;
              RobotContactInfoModi[2].StatusSwitch(0);
              AllContactStatusObj.ContactStatusAppender(RobotContactInfoModi, 2, 0);
            }
            break;
            default:
            break;
          }
          switch (RobotContactInfo[3].LocalContactStatus[0])
          {
            case 1:
            {
              std::vector<ContactStatusInfo> RobotContactInfoModi = RobotContactInfo;
              RobotContactInfoModi[3].StatusSwitch(0);
              AllContactStatusObj.ContactStatusAppender(RobotContactInfoModi, 3, 0);
            }
            break;
            default:
            break;
          }
        }
        break;
        default:
        {
          // More general case where two feet contact is shown while hands may be involved.
          for (int i = 0; i < 4; i++)
          {
            switch (RobotContactInfo[i].LocalContactStatus[0])
            {
              case 1:
              {
                std::vector<ContactStatusInfo> RobotContactInfoModi = RobotContactInfo;
                RobotContactInfoModi[i].StatusSwitch(0);
                AllContactStatusObj.ContactStatusAppender(RobotContactInfoModi, i, 0);
              }
              break;
              default:
              break;
            }
          }
        }
        break;
      }
    }
    break;
  }

  // Contact Addition
  for (int i = 0; i < RobotContactInfo.size(); i++)
  {
    if(!RobotContactInfo[i].LocalContactStatus[0]){
      std::vector<ContactStatusInfo> RobotContactInfoTemp = RobotContactInfo;
      AllContactStatusObj.ContactStatusAppender(RobotContactInfoTemp, i, 1);
    }
  }
  return AllContactStatusObj;
}

// static ControlReferenceInfo ControlReferenceGeneInner(const Robot & SimRobot, const PIPInfo & PIPObj, const Vector3 & PredictedCOMPos, ReachabilityMap & RMObject, SelfLinkGeoInfo & SelfLinkGeoObj, const std::vector<LinkInfo> & RobotLinkInfo, const std::vector<ContactStatusInfo> & FixedRobotContactInfo, const int & SwingLimbIndex, const int & Type, const double & RefFailureMetric, DataRecorderInfo & DataRecorderObj)
// {
//   ControlReferenceInfo ControlReferenceObj;
//   Vector3 ContactInit;       // This is the position of the reference contact for robot's active end effector.
//   SimRobot.GetWorldPosition(RobotLinkInfo[SwingLimbIndex].AvgLocalContact, RobotLinkInfo[SwingLimbIndex].LinkIndex, ContactInit);
//   Vector3 COMPos(0.0, 0.0, 0.0), COMVel(0.0, 0.0, 0.0);
//   CentroidalState(SimRobot, COMPos, COMVel);
//   std::vector<Vector3> OptimalContact = OptimalContactSearcher(SimRobot, PIPObj, PredictedCOMPos, RMObject, RobotLinkInfo, FixedRobotContactInfo, SwingLimbIndex, RefFailureMetric, DataRecorderObj);
//   switch (OptimalContact.size())
//   {
//     case 0:
//     {
//       return ControlReferenceObj;
//     }
//     break;
//     default:
//     {
//       // Now too early to assume that FailureFlag is true.
//       bool FeasiFlag;
//       std::vector<SplineLib::cSpline3> SplineObj;
//       int OptimalContactIndex = 0;
//       while(OptimalContactIndex<OptimalContact.size())
//       {
//         Robot SimRobotInner = SimRobot;
//         Vector3 ContactGoal = OptimalContact[OptimalContactIndex];
//         Vector3 ContactGoalGrad = NonlinearOptimizerInfo::SDFInfo.SignedDistanceNormal(ContactGoal);
//         SplineObj = TransientTrajGene(SimRobotInner, SwingLimbIndex, SelfLinkGeoObj, RobotLinkInfo, ContactInit, ContactGoal, RMObject, DataRecorderObj, FeasiFlag);
//         if(FeasiFlag)
//         {
//           EndPathInfo EndPathObj(SplineObj, SwingLimbIndex);
//           InvertedPendulumInfo InvertedPendulumObj(PIPObj.theta, PIPObj.thetadot, COMPos, COMVel);
//
//           /*
//             1. At each sampled waypoints along the end effector trajectory, an end effector position is evaluated from path.
//             2. Based on robot's current configuration, an IK problem is solved to get robot's swing limb configuration.
//             3. A time-optimal executation duration is computed.
//             4. Based on that time, robot's whole-body configuration is updated with inverted pendulum model.
//             5. The whole algorithm terminates when robot's self-collision has been triggered or no feasible IK solution can be found.
//           */
//           const int sNumber = 5;                 // 6 sampled points will be extracted from EndPathObj.
//           int sIndex = 1;
//           double sDiff = 1.0/(1.0 * sNumber - 1.0);
//           double sVal = 0.0;
//           Config CurrentConfig = SimRobotInner.q;
//           CurrentConfig = YPRShifter(CurrentConfig);
//           double CurrentTime = 0.0;
//           Vector3 CurrentContactPos = ContactInit;
//
//           std::vector<double> TimeTraj;
//           std::vector<Config> ConfigTraj;
//           std::vector<Vector3> EndEffectorTraj;
//
//           TimeTraj.reserve(sNumber);
//           ConfigTraj.reserve(sNumber);
//           EndEffectorTraj.reserve(sNumber);
//
//           TimeTraj.push_back(CurrentTime);
//           ConfigTraj.push_back(CurrentConfig);
//           EndEffectorTraj.push_back(CurrentContactPos);
//
//           bool OptFlag = true;
//           while((sIndex<sNumber)&&(OptFlag == true))
//           {
//             sVal = 1.0 * sIndex * sDiff;
//             EndPathObj.s2Pos(sVal, CurrentContactPos);
//             bool LastFlag;
//             switch (sIndex)
//             {
//               case 4:   LastFlag = true;
//               break;
//               default:  LastFlag = false;
//               break;
//             }
//             std::vector<double> OptConfig = TransientOptFn(SimRobotInner, SwingLimbIndex, SelfLinkGeoObj, CurrentContactPos, RMObject, OptFlag, LastFlag);;
//             if(OptFlag)
//             {
//               // Minimum Time Estimation.
//               double CurrentTime_i = MinimumTimeEstimation(SimRobotInner, RMObject.EndEffectorLink2Pivotal[SwingLimbIndex], CurrentConfig, Config(OptConfig));
//
//               CurrentTime+=CurrentTime_i;
//               TimeTraj.push_back(CurrentTime);
//               ConfigTraj.push_back(Config(OptConfig));
//               EndEffectorTraj.push_back(CurrentContactPos);
//
//               // Then we should update the robot's CurrentConfig based on CurrentTime_i.
//               Config UpdatedConfig  = WholeBodyDynamicsIntegrator(SimRobotInner, OptConfig, PIPObj, InvertedPendulumObj, CurrentTime_i, sIndex);
//               SimRobotInner.UpdateConfig(UpdatedConfig);
//               CurrentConfig = UpdatedConfig;
//             }
//             else
//             {
//               break;
//             }
//             sIndex++;
//           }
//           // Here the inner optimiztion loop has been finished!
//           if(OptFlag)
//           {
//             ControlReferenceObj.TrajectoryUpdate(TimeTraj, ConfigTraj, EndEffectorTraj, ContactGoal, ContactGoalGrad, EndPathObj.TotalLength, Type);
//             std::vector<ContactStatusInfo> GoalContactInfo = FixedRobotContactInfo;
//             for(int i = 0; i<FixedRobotContactInfo[SwingLimbIndex].LocalContactStatus.size(); i++)
//             {
//               GoalContactInfo[SwingLimbIndex].LocalContactStatus[i] = 1;
//             }
//             ControlReferenceObj.InitContactInfo = FixedRobotContactInfo;
//             ControlReferenceObj.GoalContactInfo = GoalContactInfo;
//
//             DataRecorderObj.OptConfigs = ConfigTraj;
//
//             return ControlReferenceObj;
//           }
//         }
//         OptimalContactIndex++;
//       }
//     }
//     break;
//   }
//   return ControlReferenceObj;
// }

ControlReferenceInfo ControlReferenceGene(Robot & SimRobot,
                                          const std::vector<ContactStatusInfo> & RobotContactInfo,
                                          ReachabilityMap & RMObject,
                                          SelfLinkGeoInfo & SelfLinkGeoObj,
                                          const SimPara & SimParaObj){
  AllContactStatus AllContactStatusObj = ContactStatusExplorer(SimRobot, RobotContactInfo);
  Vector3 COMPos, COMVel;
  getCentroidalState(SimRobot, COMPos, COMVel);
  std::clock_t start_time = std::clock();          // get current time
  for (int i = 0; i < AllContactStatusObj.ContactTypeVec.size(); i++) {
    std::vector<ContactStatusInfo> curContactInfo = AllContactStatusObj.ContactStatusInfoVec[i];
    std::vector<Vector3> ActContactPos = ActiveContactFinder(SimRobot, curContactInfo);
    PIPInfo TipOverPIP = TipOverPIPGenerator(ActContactPos, COMPos, COMVel);

    // ControlReferenceInfo ControlReferenceGeneInner(SimRobot, const PIPInfo & PIPObj, const Vector3 & PredictedCOMPos, RMObject, SelfLinkGeoObj, const std::vector<ContactStatusInfo> & FixedRobotContactInfo, const int & SwingLimbIndex, const int & Type, DataRecorderInfo & DataRecorderObj);
  }
  // int ContactStatusOption = 0;
  // int LimbSuccessNo = 0;
  // std::vector<ControlReferenceInfo> RobotTrajVec;
  // ControlReferenceInfo RobotTraj;
  // std::vector<double> ExeTimeVec;
  // std::vector<int> ContactStatusOptionVec;
  // std::vector<double> EndPathLengthVec;
  //
  // for(ContactStatusOption = 0; ContactStatusOption< AllContactStatusObj.ContactStatusInfoVec.size(); ContactStatusOption++){
  //   if(ContactStatusOptionRef!=-1){
  //     if(ContactStatusOption!=ContactStatusOptionRef){
  //       continue;
  //     }
  //   }
  //   std::vector<ContactStatusInfo> RobotContactInfo = AllContactStatusObj.ContactStatusInfoVec[ContactStatusOption];
  //   int SwingLimbIndex = AllContactStatusObj.ContactStatusExplorer[ContactStatusOption];
  //   std::vector<Vector3> ActContactPos = ContactPositionFinder(SimRobot, NonlinearOptimizerInfo::RobotLinkInfo, RobotContactInfo);    // From ContactInfoActive
  //
  //   DataRecorderInfo DataRecorderObj;
  //

  //
  //   RobotTraj = ControlReferenceGenerationInner(SimRobot, TipOverPIP, PredictedCOMPos, RMObject, SelfLinkGeoObj, NonlinearOptimizerInfo::RobotLinkInfo, RobotContactInfo, SwingLimbIndex, AllContactStatusObj.ContactTypeVec[ContactStatusOption], RefFailureMetric, DataRecorderObj);
  //   RobotTraj.SwingLimbIndex = SwingLimbIndex;
  //   RobotTraj.ContactStatusOptionIndex = ContactStatusOption;
  //   double duration_time = (std::clock() - start_time)/(double)CLOCKS_PER_SEC;
  //   std::printf("Planning takes: %f ms\n", 1000.0 * duration_time);
  //   start_time = std::clock();          // get current time
  //   RobotTraj.PlanningTime = 1000.0 * duration_time;      // The unit is ms.
  //   PlanTime+= 1000.0 * duration_time;
  //   if(RobotTraj.ControlReferenceFlag)
  //   {
  //     DataRecorderObj.PlanningNo = PlanningSteps;
  //     DataRecorderObj.LimbNo = LimbSuccessNo;
  //     DataRecorderObj.DataRecorder(SpecificPath);
  //     for (int i = 0; i < DataRecorderObj.OptConfigs.size(); i++)
  //     {
  //       const string OptConfigFile = std::to_string(PlanningSteps) + "_" + std::to_string(LimbSuccessNo) + "_" + "OptConfig" + std::to_string(i) + ".config";
  //       RobotConfigWriter(DataRecorderObj.OptConfigs[i], SpecificPath, OptConfigFile);
  //     }
  //     LimbSuccessNo++;
  //
  //     RobotTrajVec.push_back(RobotTraj);
  //     ExeTimeVec.push_back(RobotTraj.PlanStateTraj.EndTime());
  //     ContactStatusOptionVec.push_back(ContactStatusOption);
  //     EndPathLengthVec.push_back(RobotTraj.PathLength);
  //   }
  // }
  // // Based on the value of the impulse, let's select the one with the lowest impulse.
  // switch (RobotTrajVec.size())
  // {
  //   case 0:
  //   {
  //     std::printf("Planning fails to find a feasible solution! \n");
  //     // Planning fails to find a feasible solution!
  //     return RobotTraj;
  //   }
  //   break;
  //   default:
  //   {
  //     PlanningInfoFileAppender(PlanningSteps, RobotTrajVec.size()-1, SpecificPath, CurTime);
  //     int RobotTrajIndex = EndEffectorSelector(ExeTimeVec, EndPathLengthVec, ContactStatusOptionVec, PreviousContactStatusIndex);
  //     std::printf("Planning successfully finds a feasible solution! \nRobot Limb Index: %d\n", NonlinearOptimizerInfo::RobotLinkInfo[RobotTrajVec[RobotTrajIndex].SwingLimbIndex].LinkIndex);
  //     RobotTraj = RobotTrajVec[RobotTrajIndex];
  //     if(PreviousContactStatusIndex!=RobotTraj.ContactStatusOptionIndex) PreviousContactStatusIndex = RobotTraj.ContactStatusOptionIndex;
  //   }
  //   break;
  // }
  // return RobotTraj;
}
