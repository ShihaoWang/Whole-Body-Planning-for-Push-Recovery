#include "CommonHeader.h"
#include "NonlinearOptimizerInfo.h"
#include <ctime>
#include <queue>
#include <algorithm>    // std::min

typedef std::pair<double, Vector3> qEle;
static int LastEndEffectorIndex = -1;

static std::vector<ContactForm> ContactStatusExplorer(Robot & SimRobot, const std::vector<ContactStatusInfo> & RobotContactInfo){
  // First part is for contact modification.
  // Second part is for contact addition.
  Vector3 COMPos(0.0, 0.0, 0.0), COMVel(0.0, 0.0, 0.0);
  getCentroidalState(SimRobot, COMPos, COMVel);

  std::vector<ContactForm> ContactFormVec;
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
              ContactFormVec.push_back(ContactForm(RobotContactInfoModi, 2, 0));
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
              ContactFormVec.push_back(ContactForm(RobotContactInfoModi, 3, 0));
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
                ContactFormVec.push_back(ContactForm(RobotContactInfoModi, i, 0));
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
      ContactFormVec.push_back(ContactForm(RobotContactInfoTemp, i, 1));
    }
  }
  return ContactFormVec;
}

static std::vector<std::pair<Vector3, double>> ContactFreeInfoFn(const Robot & SimRobot, const std::vector<ContactStatusInfo> & RobotContactInfo, ReachabilityMap & RMObject){
  std::vector<std::pair<Vector3, double>> ContactFreeInfo;
  for (int i = 0; i < RobotContactInfo.size(); i++){
    if(RobotContactInfo[i].LocalContactStatus[0]){
      Vector3 LinkiPjPos;
      SimRobot.GetWorldPosition(NonlinearOptimizerInfo::RobotLinkInfo[i].AvgLocalContact,
        NonlinearOptimizerInfo::RobotLinkInfo[i].LinkIndex,
        LinkiPjPos);
        double Radius = RMObject.EndEffectorCollisionRadius[i];
        auto ContactFreeInfo_i = std::make_pair (LinkiPjPos, Radius);
        ContactFreeInfo.push_back(ContactFreeInfo_i);
      }
    }
    return ContactFreeInfo;
  }

static std::vector<Vector3> SupportContactFinder(const Vector3 & COMPos, const PIPInfo & PIPObj, const std::vector<Vector3> & ContactFreeContact){
  Vector3 EdgeA = PIPObj.edge_a;
  Vector3 EdgeB = PIPObj.edge_b;
  Vector3 EdgeDir = EdgeB - EdgeA;
  std::vector<Vector3> SupportContact;
  SupportContact.reserve(ContactFreeContact.size());
  for (int i = 0; i < ContactFreeContact.size(); i++){
    Vector3 ContactFreePoint = ContactFreeContact[i];
    Vector3 PointNormDir = NonlinearOptimizerInfo::SDFInfo.SignedDistanceNormal(ContactFreePoint);
    Vector3 rPos2COMPos = ContactFreePoint - COMPos;
    Vector3 InducedMomentum = cross(rPos2COMPos, PointNormDir);
    double ProjMomentumVal = InducedMomentum.x * EdgeDir.x + InducedMomentum.y * EdgeDir.y + InducedMomentum.z * EdgeDir.z;
    if(ProjMomentumVal<=0) SupportContact.push_back(ContactFreePoint);
  }
  return SupportContact;
}

static std::vector<Vector3> OptimalContactFinder(const std::vector<Vector3> & SupportContact, const std::vector<Vector3> & FixedContacts, const Vector3 & COMPos, const Vector3 & COMVel, int CutOffNo, SimPara & SimParaObj){
  // This function selects the optimal contact given support contact.
  std::priority_queue<qEle, std::vector<qEle>, less<qEle> > OptimalContactQueue;
  std::vector<Vector3> CandidateContacts;
  std::vector<Vector3> CandidateContactWeights;
  std::vector<double> ContactFailureMetric(SupportContact.size());
  Vector3 CurContact = SimParaObj.getContactInit();
  const int ActContactNo = FixedContacts.size() + 1;
  for (int i = 0; i < SupportContact.size(); i++){
    std::vector<Vector3> ActContacts = FixedContacts;
    ActContacts.push_back(SupportContact[i]);
    std::vector<PIPInfo> PIPTotal = PIPGenerator(ActContacts, COMPos, COMVel);
    ContactFailureMetric[i] = FailureMetricEval(PIPTotal);
    if(ContactFailureMetric[i]>0.0){
      printf("ContactFailureMetric[%d]: %f\n", OptimalContactQueue.size(), ContactFailureMetric[i]);
      Vector3 ContactDiff = CurContact - SupportContact[i];
      double ContactDiffDist = ContactDiff.normSquared();
      double ContactDistCost = 1.0 * exp(-1.25 * ContactDiffDist);
      ContactFailureMetric[i]*= ContactDistCost;
      OptimalContactQueue.push(std::make_pair(ContactFailureMetric[i], SupportContact[i]));
      CandidateContacts.push_back(SupportContact[i]);
      CandidateContactWeights.push_back(ContactFailureMetric[i] * NonlinearOptimizerInfo::SDFInfo.SignedDistanceNormal(SupportContact[i]));
    }
  }

  std::vector<Vector3> SelectedContacts;
  if(!CandidateContacts.size()){
    int OptiIndex = std::distance(ContactFailureMetric.begin(), std::max_element(ContactFailureMetric.begin(), ContactFailureMetric.end()));
    CandidateContacts.push_back(SupportContact[OptiIndex]);
    CandidateContactWeights.push_back(0.0);
    SelectedContacts.push_back(SupportContact[OptiIndex]);
  } else {
    int OptimalContactNumber = CandidateContacts.size();
    int OptEleNo = std::min(OptimalContactNumber, CutOffNo);
    for (int i = 0; i < OptEleNo; i++) {
      SelectedContacts.push_back(OptimalContactQueue.top().second);
      OptimalContactQueue.pop();
    }
  }
  // Vector3Writer(CandidateContacts, "OptimalContact");
  // Vector3Writer(CandidateContactWeights, "OptimalContactWeights");
  SimParaObj.DataRecorderObj.setCCSData(CandidateContacts, CandidateContactWeights, SelectedContacts);
  return SelectedContacts;
}

static std::vector<Vector3> OptimalContactSearcher( Robot SimRobot,     const PIPInfo & PIPObj,
  ReachabilityMap & RMObject, const ContactForm & ContactFormObj,
  SimPara & SimParaObj){
    std::vector<Vector3> OptimalContact;
    std::vector<Vector3> FixedContactPos;
    for (int i = 0; i < NonlinearOptimizerInfo::RobotLinkInfo.size(); i++){
      for (int j = 0; j < NonlinearOptimizerInfo::RobotLinkInfo[i].LocalContacts.size(); j++){
        if(ContactFormObj.FixedContactStatusInfo[i].LocalContactStatus[j]){
          Vector3 LinkiPjPos;
          SimRobot.GetWorldPosition(NonlinearOptimizerInfo::RobotLinkInfo[i].LocalContacts[j],
            NonlinearOptimizerInfo::RobotLinkInfo[i].LinkIndex,
            LinkiPjPos);
            FixedContactPos.push_back(LinkiPjPos);
          }
        }
      }
      std::vector<std::pair<Vector3, double>> ContactFreeInfoVec = ContactFreeInfoFn(SimRobot, ContactFormObj.FixedContactStatusInfo, RMObject);

      Vector3 COMPos, COMVel;
      getCentroidalState(SimRobot, COMPos, COMVel);
      InvertedPendulumInfo InvertedPendulumObj(PIPObj.L, PIPObj.g, PIPObj.theta, PIPObj.thetadot, COMPos, COMVel);
      InvertedPendulumObj.setEdges(PIPObj.edge_a, PIPObj.edge_b);

      Config UpdatedConfig  = WholeBodyDynamicsIntegrator(SimRobot, InvertedPendulumObj, SimParaObj.ForwardDuration);
      SimRobot.UpdateConfig(UpdatedConfig);

      COMPos = InvertedPendulumObj.COMPos;
      COMVel = InvertedPendulumObj.COMVel;

      //  0. Reachable with respect to the pivotal joint
      std::vector<Vector3> ReachableContacts = RMObject.ReachablePointsFinder(SimRobot, ContactFormObj.SwingLinkInfoIndex, NonlinearOptimizerInfo::SDFInfo);
      if(!ReachableContacts.size()) return OptimalContact;
      // 1. Self-collision from other end effectors
      std::vector<Vector3> CollisionFreeContacts = RMObject.ContactFreePointsFinder(RMObject.EndEffectorCollisionRadius[ContactFormObj.SwingLinkInfoIndex], ReachableContacts, ContactFreeInfoVec);
      if(!CollisionFreeContacts.size()) return OptimalContact;
      // 2. Supportive
      std::vector<Vector3> SupportiveContacts = SupportContactFinder(COMPos, PIPObj, CollisionFreeContacts);
      if(!SupportiveContacts.size()) return OptimalContact;

      SimParaObj.DataRecorderObj.setRCSData(ReachableContacts, CollisionFreeContacts, SupportiveContacts);

      // 3. Optimal Contact
      int CutOffNo = 5;
      OptimalContact = OptimalContactFinder(SupportiveContacts, FixedContactPos, COMPos, COMVel, CutOffNo, SimParaObj);
      if(!OptimalContact.size()) return OptimalContact;

      return OptimalContact;
    }

static ControlReferenceInfo ControlReferenceGeneInner(const Robot & SimRobot, const PIPInfo & TipOverPIPObj, ReachabilityMap & RMObject, SelfLinkGeoInfo & SelfLinkGeoObj, const ContactForm & ContactFormObj, SimPara & SimParaObj){
  ControlReferenceInfo ControlReferenceObj;
  Vector3 ContactInit;       // This is the position of the reference contact for robot's active end effector.
  SimRobot.GetWorldPosition(NonlinearOptimizerInfo::RobotLinkInfo[SimParaObj.SwingLinkInfoIndex].AvgLocalContact,
  NonlinearOptimizerInfo::RobotLinkInfo[SimParaObj.SwingLinkInfoIndex].LinkIndex, ContactInit);
  SimParaObj.setContactInit(ContactInit);
  std::vector<Vector3> OptimalContact = OptimalContactSearcher(SimRobot, TipOverPIPObj, RMObject, ContactFormObj, SimParaObj);
  if(!OptimalContact.size()) return ControlReferenceObj;
  SimParaObj.setFixedContactStatusInfo(ContactFormObj.FixedContactStatusInfo);

  Vector3 COMPos, COMVel;
  getCentroidalState(SimRobot, COMPos, COMVel);
  InvertedPendulumInfo InvertedPendulumObj( TipOverPIPObj.L, TipOverPIPObj.g,
                                            TipOverPIPObj.theta, TipOverPIPObj.thetadot,
                                            COMPos, COMVel);
  InvertedPendulumObj.setEdges(TipOverPIPObj.edge_a, TipOverPIPObj.edge_b);
  for (int i = 0; i < OptimalContact.size(); i++) {
    Robot SimRobotInner = SimRobot;
    SimParaObj.setContactGoal(OptimalContact[i]);
    SimParaObj.setTransPathFeasiFlag(false);
    CubicSplineInfo CubicSplineInfoObj = TransientPathGene(SimRobotInner, SelfLinkGeoObj, SimParaObj);
    if(SimParaObj.getTransPathFeasiFlag()){
      EndEffectorPathInfo EndEffectorPathObj(CubicSplineInfoObj);
      // Here two methods will be conducted for comparison purpose.
      /*
      1. At each sampled waypoints along the end effector trajectory, an end effector position is evaluated from path.
      2. Based on robot's current configuration, an IK problem is solved to get robot's swing limb configuration.
      3. A time-optimal executation duration is computed.
      4. Based on that time, robot's whole-body configuration is updated with inverted pendulum model.
      5. The whole algorithm terminates when robot's self-collision has been triggered or no feasible IK solution can be found.
      */
      ControlReferenceObj = TrajectoryPlanning(SimRobotInner, InvertedPendulumObj, RMObject, SelfLinkGeoObj,
        EndEffectorPathObj, SimParaObj);
        if(ControlReferenceObj.getReadyFlag()) break;
      }
    }
    return ControlReferenceObj;
  }

ControlReferenceInfo ControlReferenceGene(Robot & SimRobot,
                                          const std::vector<ContactStatusInfo> & RobotContactInfo,
                                          ReachabilityMap & RMObject,
                                          SelfLinkGeoInfo & SelfLinkGeoObj,
                                          SimPara & SimParaObj){
 std::vector<ContactForm> ContactFormVec = ContactStatusExplorer(SimRobot, RobotContactInfo);
 Vector3 COMPos, COMVel;
 getCentroidalState(SimRobot, COMPos, COMVel);
 ControlReferenceInfo ControlReferenceInfoObj;
 ControlReferenceInfoObj.setReadyFlag(false);
 std::vector<ControlReferenceInfo> ControlReferenceObjVec;
 std::vector<double> ExecutionTimeVec;
 std::vector<double> EstFailureMetricVec;
 int PlanEndEffectorIndex = 0;
 std::clock_t start_time = std::clock();
 double stage_planning_time = 0.0;
 for (auto & ContactFormObj : ContactFormVec) {
   std::vector<ContactStatusInfo> curContactInfo = ContactFormObj.FixedContactStatusInfo;
   std::vector<Vector3> ActContactPos = ActiveContactFinder(SimRobot, curContactInfo);
   bool ValidFlag;
   PIPInfo TipOverPIP = TipOverPIPGenerator(ActContactPos, COMPos, COMVel, ValidFlag);
   if(!ValidFlag) continue;
   SimParaObj.setSwingLinkInfoIndex(ContactFormObj.SwingLinkInfoIndex);
   SimParaObj.setPlanEndEffectorIndex(PlanEndEffectorIndex);
   SimParaObj.DataRecorderObj.setPlanStageIndexNLinkNo(SimParaObj.getPlanStageIndex(), PlanEndEffectorIndex);
   ControlReferenceInfo ControlReferenceObj =  ControlReferenceGeneInner(SimRobot, TipOverPIP, RMObject, SelfLinkGeoObj, ContactFormObj, SimParaObj);
   double planning_time = (std::clock() - start_time)/(double)CLOCKS_PER_SEC;
   std::printf("Planning takes: %f ms\n", 1000.0 * planning_time);
   stage_planning_time+=planning_time;
   start_time = std::clock();
   ControlReferenceObj.setSwingLinkInfoIndex(ContactFormObj.SwingLinkInfoIndex);
   ControlReferenceObj.ControlReferenceType = ContactFormObj.ContactType;
   if(ControlReferenceObj.getReadyFlag()){
     ControlReferenceObjVec.push_back(ControlReferenceObj);
     ExecutionTimeVec.push_back(ControlReferenceObj.TimeTraj.back());
     EstFailureMetricVec.push_back(ControlReferenceObj.getFailueMetric());
     SimParaObj.DataRecorderObj.PlannedConfigTraj = ControlReferenceObj.PlannedConfigTraj;
     SimParaObj.DataRecorderObj.EndEffectorTraj = ControlReferenceObj.EndEffectorTraj;
     SimParaObj.DataRecorderObj.Write2File(SimParaObj.getCurrentCasePath());
     PlanEndEffectorIndex++;
   }
 }
 // // Select Executation Time
 // if(ExecutionTimeVec.size()){
 //   int ObjIndex = std::distance(ExecutionTimeVec.begin(), std::min_element(ExecutionTimeVec.begin(), ExecutionTimeVec.end()));
 //   ControlReferenceInfoObj = ControlReferenceObjVec[ObjIndex];
 //   PlanTimeRecorder(stage_planning_time, SimParaObj.getCurrentCasePath());
 //   PlanningInfoFileAppender(SimParaObj.getPlanStageIndex(), ExecutionTimeVec.size()-1, SimParaObj.getCurrentCasePath(), SimParaObj.getSimTime());
 // }
 // Select Estimated Failure Metric
 if(EstFailureMetricVec.size()){
   int ObjIndex;
   if(LastEndEffectorIndex == -1){
     ObjIndex = std::distance(EstFailureMetricVec.begin(), std::max_element(EstFailureMetricVec.begin(), EstFailureMetricVec.end()));
    //  ObjIndex = 1;
     ControlReferenceInfoObj = ControlReferenceObjVec[ObjIndex];
     LastEndEffectorIndex = ControlReferenceInfoObj.getSwingLinkInfoIndex();
   }
   else
   {
     ObjIndex = std::distance(EstFailureMetricVec.begin(), std::max_element(EstFailureMetricVec.begin(), EstFailureMetricVec.end()));
     ControlReferenceInfoObj = ControlReferenceObjVec[ObjIndex];
     if(ControlReferenceInfoObj.getSwingLinkInfoIndex() == LastEndEffectorIndex){
       EstFailureMetricVec[ObjIndex] = -1.0;
       ObjIndex = std::distance(EstFailureMetricVec.begin(), std::max_element(EstFailureMetricVec.begin(), EstFailureMetricVec.end()));
       ControlReferenceInfoObj = ControlReferenceObjVec[ObjIndex];
       LastEndEffectorIndex = ControlReferenceInfoObj.getSwingLinkInfoIndex();
     } else
        LastEndEffectorIndex = ControlReferenceInfoObj.getSwingLinkInfoIndex();
   }
   PlanTimeRecorder(stage_planning_time, SimParaObj.getCurrentCasePath());
   PlanningInfoFileAppender(SimParaObj.getPlanStageIndex(), EstFailureMetricVec.size()-1, SimParaObj.getCurrentCasePath(), SimParaObj.getSimTime());
 }
 return ControlReferenceInfoObj;
}
