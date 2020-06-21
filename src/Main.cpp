#include <ctime>
#include "CommonHeader.h"
#include "NonlinearOptimizerInfo.h"
#include <omp.h>
#include "Control/PathController.h"
#include "Simulation/WorldSimulation.h"
#include <ode/ode.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include "RobotInfo.h"

std::vector<LinkInfo>   NonlinearOptimizerInfo::RobotLinkInfo;
SignedDistanceFieldInfo NonlinearOptimizerInfo::SDFInfo;
AnyCollisionGeometry3D  NonlinearOptimizerInfo::TerrColGeom;

int main(){
  /* 1. Load the Contact Link file */
  const std::string UserFilePath = "../user/";
  const std::string ContactLinkPath = UserFilePath + "ContactLink.txt";
  int NumberOfContactPoints;
  NonlinearOptimizerInfo::RobotLinkInfo = ContactInfoLoader(ContactLinkPath, NumberOfContactPoints);
  const std::string TorsoLinkFilePath = UserFilePath + "TorsoLink.txt";
  std::vector<int> TorsoLink = TorsoLinkReader(TorsoLinkFilePath);
  const std::string SelfCollisionFreeLinkFilePath = UserFilePath + "SelfCollisionFreeLink.txt";
  std::vector<int> SelfCollisionFreeLink = TorsoLinkReader(SelfCollisionFreeLinkFilePath);

  /* 2. Load the Envi and Contact Status file */
  std::ifstream FolderPathFile("../user/CaseSpecs.txt");
  std::string ExpName;
  std::getline(FolderPathFile, ExpName);
  FolderPathFile.close();

  const std::string ExperimentFolderPath = "/home/motion/Desktop/Whole-Body-Planning-for-Push-Recovery-Data/result/" + ExpName + "/";

  std::ifstream SettingsFile(ExperimentFolderPath + "Settings.txt");
  std::string ForceMaxStr, PushDurationStr, DetectionWaitStr;
  std::getline(SettingsFile, ForceMaxStr);
  std::getline(SettingsFile, PushDurationStr);
  std::getline(SettingsFile, DetectionWaitStr);
  SettingsFile.close();
  double ForceMax       = std::stod(ForceMaxStr);
  double PushDuration   = std::stod(PushDurationStr);
  double DetectionWait  = std::stod(DetectionWaitStr);
  double TimeStep       = 0.025;
  double InitDuration    = 2.0;
  double TotalDuration   = 5.0;                 // Simulation lasts for 5s after initial duration
  double ForwardDuartion = 0.75;                // Used to optimal contact point planning
  double PhaseRatio     = 0.75;
  double PhaseTimeStep  = 0.05;
  double ReductionRatio = 0.2;
  SimPara SimParaObj(ForceMax, PushDuration, DetectionWait, TimeStep, InitDuration, TotalDuration, ForwardDuartion, PhaseRatio, PhaseTimeStep, ReductionRatio);

  RobotWorld worldObj;
  SimGUIBackend BackendObj(&worldObj);
  WorldSimulation& SimObj = BackendObj.sim;
  string XMLFileStrObj =  ExperimentFolderPath + "Environment.xml";
  const char* XMLFileObj = XMLFileStrObj.c_str();    // Here we must give abstract path to the file
  if(!BackendObj.LoadAndInitSim(XMLFileObj)){
    std::cerr<< XMLFileStrObj<<" file does not exist in that path!"<<endl;
    return -1;
  }

  /* 3. Environment Geometry, Reachability Map and BBs for Self-collision Test */
  const int GridsNo = 251;
  struct stat buffer;   // This is used to check whether "SDFSpecs.bin" exists or not.
  const string SDFPath = ExperimentFolderPath + "SDFs/";
  const string SDFSpecsName = SDFPath + "SDFSpecs.bin";
  if(stat (SDFSpecsName.c_str(), &buffer) == 0){
    NonlinearOptimizerInfo::SDFInfo = SignedDistanceFieldLoader(SDFPath, GridsNo);
  } else {
      NonlinearOptimizerInfo::SDFInfo = SignedDistanceFieldGene(SDFPath, worldObj, GridsNo);
  }
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
  SelfLinkGeoInfo SelfLinkGeoObj(*worldObj.robots[0], RMObject.EndEffectorLink2Pivotal, SelfCollisionFreeLink);

  /* 5. Internal Experimentation Loop */
  int TotalNumber = 100;
  int FileIndex = FileIndexFinder(false);  FileIndex = 1;
  while(FileIndex<=TotalNumber){
    RobotWorld world;
    SimGUIBackend Backend(&world);
    WorldSimulation& Sim = Backend.sim;

    string XMLFileStr =  ExperimentFolderPath + "Environment.xml";
    const char* XMLFile = XMLFileStr.c_str();    // Here we must give abstract path to the file
    if(!Backend.LoadAndInitSim(XMLFile)){
      std::cerr<<XMLFileStr<< "file does not exist in that path!" << endl;
      return -1;
    }
    Robot SimRobot = *world.robots[0];
    string CurrentCasePath = ExperimentFolderPath + std::to_string(FileIndex) + "/";
    RobotConfigLoader(SimRobot, CurrentCasePath, "InitConfig.config");
    const std::string ContactStatusPath = CurrentCasePath + "ContactStatus.txt";
    std::vector<ContactStatusInfo> InitContactInfo = ContactStatusInfoLoader(ContactStatusPath);
    SimParaObj.CurrentCasePathUpdate(CurrentCasePath);

    std::vector<double> InitConfig(SimRobot.q);
    std::vector<double> InitVelocity(SimRobot.q.size(), 0.0);
    std::vector<double> RobotConfigRef = InitVelocity;

    //  Given the optimized result to be the initial state
    Sim.world->robots[0]->UpdateConfig(Config(InitConfig));
    Sim.world->robots[0]->dq = InitVelocity;
    Sim.controlSimulators[0].oderobot->SetConfig(Config(InitConfig));
    Sim.controlSimulators[0].oderobot->SetVelocities(Config(InitVelocity));

    Vector3 ImpulseDirection = ImpulseDirectionGene(*Sim.world->robots[0], InitContactInfo, 1);
    SimParaObj.setImpulseForceMax(ImpulseDirection);
    FilePathManager(SimParaObj.CurrentCasePath);
    int SimRes = SimulationTest(Sim, InitContactInfo, RMObject, SelfLinkGeoObj, SimParaObj);
    // int SimulationTest(Sim, InitContactInfo, RMObject, SelfLinkGeoObj)

    // CurrentCasePath+= PlanningType + "/";
    // SimulationTest(Sim, InitContactInfo, SelfLinkGeoObj, SimParaObj);
    // if(FailureFlag)
    // {
    //   PlanResWriter(CurrentCasePath, PushRecovFlag);
    //   FileIndex = FileIndexFinder(true);
    //   FileIndex++;
    // }
  }
  return 1;
}
