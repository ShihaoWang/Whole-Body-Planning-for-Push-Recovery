#include "CommonHeader.h"
#include "NonlinearOptimizerInfo.h"
#include "Control/PathController.h"
#include "Simulation/WorldSimulation.h"
#include <ode/ode.h>
#include "RobotInfo.h"

std::vector<LinkInfo>   NonlinearOptimizerInfo::RobotLinkInfo;
SignedDistanceFieldInfo NonlinearOptimizerInfo::SDFInfo;
AnyCollisionGeometry3D  NonlinearOptimizerInfo::TerrColGeom;

static void mainInner(string ExperimentFolderPath, int FileIndex, string Type, int RandIter, string ForceMax, ReachabilityMap & RMObject, SelfLinkGeoInfo & SelfLinkGeoObj, SimPara & SimParaObj){
  
  RobotWorld world;
  SimGUIBackend Backend(&world);
  WorldSimulation& Sim = Backend.sim;
  
  string XMLFileStr =  ExperimentFolderPath + "Environment.xml";
  const char* XMLFile = XMLFileStr.c_str();    // Here we must give abstract path to the file
  if(!Backend.LoadAndInitSim(XMLFile)){
    std::cerr<<XMLFileStr<< "file does not exist in that path!" << endl;
    return;
  }

  Robot SimRobot = *world.robots[0];
  string CurrentCasePath = ExperimentFolderPath + std::to_string(FileIndex) + "/";
  RobotConfigLoader(SimRobot, CurrentCasePath, "InitConfig.config");
  const std::string ContactStatusPath = CurrentCasePath + "ContactStatus.txt";
  std::vector<ContactStatusInfo> InitContactInfo = ContactStatusInfoLoader(ContactStatusPath);

  if(RandIter<0){
    CurrentCasePath+= Type + "/" + ForceMax + "/";
    // struct stat buffer; 
    // if(stat (CurrentCasePath.c_str(), &buffer) != 0){
    string mkdir_call = "mkdir " + CurrentCasePath;
    std::system(mkdir_call.c_str()); // note the slash after accounts!
    // }
  }
  else {
    string CurrentCasePathTemp = CurrentCasePath + Type + "/" + ForceMax + "/" + to_string(RandIter) + "/";
    // struct stat buffer; 
    // if(stat (CurrentCasePathTemp.c_str(), &buffer) != 0){
    string mkdir_call = "mkdir " + CurrentCasePath +  Type + "/" + ForceMax + "/";
    std::system(mkdir_call.c_str()); // note the slash after accounts!
    CurrentCasePath+= Type + "/" + ForceMax + "/" + to_string(RandIter) + "/";
    mkdir_call = "mkdir " + CurrentCasePath;
    std::system(mkdir_call.c_str()); // note the slash after accounts!
    // }
  }

  SimParaObj.CurrentCasePathUpdate(CurrentCasePath);
  FilePathManager(SimParaObj.CurrentCasePath);

  std::vector<double> InitConfig(SimRobot.q);
  std::vector<double> InitVelocity(SimRobot.q.size(), 0.0);
  std::vector<double> RobotConfigRef = InitVelocity;

  Sim.world->robots[0]->UpdateConfig(Config(InitConfig));
  Sim.world->robots[0]->dq = InitVelocity;
  Sim.controlSimulators[0].oderobot->SetConfig(Config(InitConfig));
  Sim.controlSimulators[0].oderobot->SetVelocities(Config(InitVelocity));
  
  // Enable terrain robot links contact
  const int NumberOfTerrains = Sim.world->terrains.size();
  for (int i = 0; i < NumberOfTerrains; i++){
    int terrainid = Sim.world->TerrainID(i);
    for (int j = 0; j < Sim.world->robots[0]->links.size(); j++){
      Sim.EnableContactFeedback(terrainid, world.RobotLinkID(0, j));    
  }
  }
  int SimRes = SimulationTest(Sim, InitContactInfo, RMObject, SelfLinkGeoObj, SimParaObj);
  PlanResWriter(CurrentCasePath, SimRes);
  return;
}

int main(){
  // A high-level simulation program where a number of parameters can be set.

   /* 1. Load the Contact Link file */
  const std::string UserFilePath = "../user/";
  const std::string ContactLinkPath = UserFilePath + "ContactLink.txt";
  int NumberOfContactPoints;
  NonlinearOptimizerInfo::RobotLinkInfo = ContactInfoLoader(ContactLinkPath, NumberOfContactPoints);
  const std::string TorsoLinkFilePath = UserFilePath + "TorsoLink.txt";
  std::vector<int> TorsoLink = TorsoLinkReader(TorsoLinkFilePath);
  const std::string SelfCollisionFreeLinkFilePath = UserFilePath + "SelfCollisionFreeLink.txt";
  std::vector<int> SelfCollisionFreeLink = TorsoLinkReader(SelfCollisionFreeLinkFilePath);

  /* 
    Two main types of experimentation: 
      1.  Grad: Given a fixed impulse direction,  we gradually increase initial disturbance forces to check out how algorithm behaves.
      2.  Random: Given a fixed force magnitude,    the impulsive force is imposed in randomized direction.
  */

  string Type = "Grad";
  // string Type = "Random";
    
  std::vector<int> DisturbForceVec1Contact{ 1500, 2000, 2500, 3000, 3500, 4000, 4500}; 
  std::vector<int> DisturbForceVec2Contact{ 3500, 4000, 4500, 5000, 0500, 6000, 6500}; 
  
  const int RandTotal     = 10;    // How many randomized impulses will be tested
  int Rand1Contact     = 3500;
  int Rand2Contact     = 4500;

  double PushDuration     = 0.25;
  double DetectionWait    = 0.25;

  // Three inner variables
  double TimeStep         = 0.025;
  double InitDuration     = 2.0;
  double TotalDuration    = 5.0;                     // Simulation lasts for 5s after initial duration

  double ForwardDuration  = 0.5;
  double PhaseRatio       = 0.6;
  double ReductionRatio   = 0.6;

  // Three horizontal directions to be chosen from.
  Vector3 FixedDirectionX(1.0, 0.0, 0.0);
  Vector3 FixedDirectionY(0.0, 1.0, 0.0);
  Vector3 FixedDirectionXY(1.0, 1.0, 0.0);    FixedDirectionXY.setNormalized(FixedDirectionXY);

  // std::vector<std::string> ScenarioVec = { "flat_1Contact", "flat_2Contact", "uneven_1Contact", "uneven_2Contact"};  
  std::vector<std::string> ScenarioVec = {"flat_2Contact", "uneven_1Contact", "uneven_2Contact"};  

  const int GradTotal = 25;

  for (auto Scenario: ScenarioVec){
    const std::string ExperimentFolderPath = "/home/motion/Desktop/Whole-Body-Planning-for-Push-Recovery-Data/result/" + Scenario + "/";
    RobotWorld worldObj;
    SimGUIBackend BackendObj(&worldObj);
    WorldSimulation& SimObj = BackendObj.sim;
    string XMLFileStrObj =  ExperimentFolderPath + "Environment.xml";
    const char* XMLFileObj = XMLFileStrObj.c_str();    // Here we must give abstract path to the file
    if(!BackendObj.LoadAndInitSim(XMLFileObj)){
      std::cerr<< XMLFileStrObj<<" file does not exist in that path!"<<endl;
      return -1;
    }


    const int GridsNo = 251;
    struct stat buffer;   // This is used to check whether "SDFSpecs.bin" exists or not.
    const string SDFPath = ExperimentFolderPath + "SDFs/";
    const string SDFSpecsName = SDFPath + "SDFSpecs.bin";
    if(stat (SDFSpecsName.c_str(), &buffer) == 0)
      NonlinearOptimizerInfo::SDFInfo = SignedDistanceFieldLoader(SDFPath, GridsNo);
    else
      NonlinearOptimizerInfo::SDFInfo = SignedDistanceFieldGene(SDFPath, worldObj, GridsNo);

    ReachabilityMap RMObject = ReachabilityMapGenerator(*worldObj.robots[0], NonlinearOptimizerInfo::RobotLinkInfo, TorsoLink);
    const int NumberOfTerrains = worldObj.terrains.size();
    std::shared_ptr<Terrain> Terrain_ptr = std::make_shared<Terrain>(*worldObj.terrains[0]);
    Meshing::TriMesh EnviTriMesh  = Terrain_ptr->geometry->AsTriangleMesh();

    for (int i = 0; i < NumberOfTerrains-1; i++){
      std::shared_ptr<Terrain> Terrain_ptr = std::make_shared<Terrain>(*worldObj.terrains[i+1]);
      Meshing::TriMesh EnviTriMesh_i  = Terrain_ptr->geometry->AsTriangleMesh();
      EnviTriMesh.MergeWith(EnviTriMesh_i);
    }
    AnyCollisionGeometry3D TerrColGeom(EnviTriMesh);
    NonlinearOptimizerInfo::TerrColGeom = TerrColGeom;

    if(Type == "Grad"){
      int FileIndex = FileIndexFinder(false);
      while(FileIndex<=GradTotal){
      std::vector<int> DisturbForceVec;
        if (Scenario.find('1') != string::npos) DisturbForceVec = DisturbForceVec1Contact;
        else DisturbForceVec = DisturbForceVec2Contact;
        for (int i = 0; i < DisturbForceVec.size(); i++){
          SimPara SimParaObj( 1.0 * DisturbForceVec[i], PushDuration, DetectionWait, 
                              TimeStep, InitDuration, TotalDuration, 
                              ForwardDuration, PhaseRatio, ReductionRatio);
          SimParaObj.setImpulseForceMax(FixedDirectionX);
          SelfLinkGeoInfo SelfLinkGeoObj(*worldObj.robots[0], RMObject.EndEffectorLink2Pivotal, SelfCollisionFreeLink);
          mainInner(ExperimentFolderPath, FileIndex, Type, -1, std::to_string(DisturbForceVec[i]), RMObject, SelfLinkGeoObj, SimParaObj);
        }
        FileIndex = FileIndexFinder(true);
        FileIndex++;  
      }
    } else {
      int FileIndex = FileIndexFinder(false);
      while(FileIndex<=GradTotal){
        int  DisturbForce = 0.0;
        if (Scenario.find('1') != string::npos) DisturbForce = Rand1Contact;
        else DisturbForce = Rand2Contact;
        for (int RandIter = 0; RandIter < RandTotal; RandIter++){
          SimPara SimParaObj( 1.0 * DisturbForce, PushDuration, DetectionWait, 
                              TimeStep, InitDuration, TotalDuration, 
                              ForwardDuration, PhaseRatio, ReductionRatio);
          SimParaObj.setImpulseForceMax(FlatRandomDirection());
          SelfLinkGeoInfo SelfLinkGeoObj(*worldObj.robots[0], RMObject.EndEffectorLink2Pivotal, SelfCollisionFreeLink);
          mainInner(ExperimentFolderPath, FileIndex, Type, RandIter, std::to_string(DisturbForce), RMObject, SelfLinkGeoObj, SimParaObj);
        }
        FileIndex = FileIndexFinder(true);
        FileIndex++;  
      }
    }
  }
  //   switch (Type){
  //   case "Grad":{
  //     SimParaObj.setImpulseForceMax(FixedDirectionX);
  //     int FileIndex = FileIndexFinder(false);
  //     while(FileIndex<=GridsNo){
  //     std::vector<int> DisturbForceVec;
  //       if (Scenario.find('1') != string::npos) DisturbForceVec = DisturbForceVec1Contact;
  //       else DisturbForceVec = DisturbForceVec2Contact;
  //       for (int i = 0; i < DisturbForceVec.size(); i++){
  //         SimPara SimParaObj( 1.0 * stod(DisturbForceVec[i]), PushDuration, DetectionWait, 
  //                             TimeStep, InitDuration, TotalDuration, 
  //                             FowardDuration, PhaseRatio, ReductionRatio);
  //         SelfLinkGeoInfo SelfLinkGeoObj(*worldObj.robots[0], RMObject.EndEffectorLink2Pivotal, SelfCollisionFreeLink);
  //         mainInner(ExperimentFolderPath, FileIndex, Type, -1, DisturbForceVec[i], RMObject, SelfLinkGeoObj, SimParaObj);
  //       }
  //       FileIndex = FileIndexFinder(true);
  //       FileIndex++;
  //     }      
  //   } break;
  //   default:{
  //     if (Scenario.find('1') != string::npos) {
  //         // 1 Contact Case
  //         ForceMax = Rand1Contact;
  //       } else {
  //         // 2 Contact Case
  //         ForceMax = Rand2Contact;
  //       }
  //   } break;
  //   }
  // }
  return 1;
}