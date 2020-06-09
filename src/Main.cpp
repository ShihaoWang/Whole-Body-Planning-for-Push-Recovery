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

SignedDistanceFieldInfo NonlinearOptimizerInfo::SDFInfo;
std::vector<LinkInfo>   NonlinearOptimizerInfo::RobotLinkInfo;
static bool SDFFlag = false;

int main()
{
  /* 1. Load the Contact Link file */
  const std::string UserFilePath = "../user/hrp2/";
  const std::string ContactLinkPath = UserFilePath + "ContactLink.txt";
  int NumberOfContactPoints;
  NonlinearOptimizerInfo::RobotLinkInfo = ContactInfoLoader(ContactLinkPath, NumberOfContactPoints);
  const std::string TorsoLinkFilePath = UserFilePath + "TorsoLink.txt";
  std::vector<int> TorsoLink = TorsoLinkReader(TorsoLinkFilePath);
  const std::string SelfCollisionFreeLinkFilePath = UserFilePath + "SelfCollisionFreeLink.txt";
  std::vector<int> SelfCollisionFreeLink = TorsoLinkReader(SelfCollisionFreeLinkFilePath);

  /* 2. Load the Envi and Contact Status file */
  std::ifstream FolderPathFile("./Specs/EnviSpecs.txt");
  std::string EnviName;
  std::string ExpName;
  std::string ContactType;
  std::getline(FolderPathFile, EnviName);
  std::getline(FolderPathFile, ExpName);
  std::getline(FolderPathFile, ContactType);
  FolderPathFile.close();

  const std::string ExperimentPath = "../result/" + ExpName + "/" + ContactType;

  std::ifstream SettingsFile(ExperimentPath + "/" + "Settings.txt");
  std::string ForceMaxStr, PushDurationStr, DetectionWaitStr;
  std::getline(SettingsFile, ForceMaxStr);
  std::getline(SettingsFile, PushDurationStr);
  std::getline(SettingsFile, DetectionWaitStr);
  SettingsFile.close();
  double ForceMax = std::stod(ForceMaxStr);
  double PushDuration = std::stod(PushDurationStr);
  double DetectionWait = std::stod(DetectionWaitStr);

  RobotWorld worldObj;
  SimGUIBackend BackendObj(&worldObj);
  WorldSimulation& SimObj = BackendObj.sim;
  string XMLFileStrObj =  "../" + EnviName;
  const char* XMLFileObj = XMLFileStrObj.c_str();    // Here we must give abstract path to the file
  if(!BackendObj.LoadAndInitSim(XMLFileObj))
  {
    std::cerr<< EnviName<<" file does not exist in that path!"<<endl;
    return -1;
  }

  /* 3. Environment Geometry and Reachability Map */
  const int GridsNo = 251;
  struct stat buffer;   // This is used to check whether "SDFSpecs.bin" exists or not.
  const string SDFSpecsName = "./SDFs/SDFSpecs.bin";
  if(stat (SDFSpecsName.c_str(), &buffer) == 0)
  {
    NonlinearOptimizerInfo::SDFInfo = SignedDistanceFieldLoader(GridsNo);
  }
  else
  {
    if(!SDFFlag)
    {
      NonlinearOptimizerInfo::SDFInfo = SignedDistanceFieldGene(worldObj, GridsNo);
      SDFFlag = true;
    }
  }
  // ReachabilityMap RMObject = ReachabilityMapGenerator(*worldObj.robots[0], NonlinearOptimizerInfo::RobotLinkInfo, TorsoLink);
  //
  // const int NumberOfTerrains = worldObj.terrains.size();
  // std::shared_ptr<Terrain> Terrain_ptr = std::make_shared<Terrain>(*worldObj.terrains[0]);
  // Meshing::TriMesh EnviTriMesh  = Terrain_ptr->geometry->AsTriangleMesh();
  // for (int i = 0; i < NumberOfTerrains-1; i++)
  // {
  //   std::shared_ptr<Terrain> Terrain_ptr = std::make_shared<Terrain>(*worldObj.terrains[i+1]);
  //   Meshing::TriMesh EnviTriMesh_i  = Terrain_ptr->geometry->AsTriangleMesh();
  //   EnviTriMesh.MergeWith(EnviTriMesh_i);
  // }
  // AnyCollisionGeometry3D TerrColGeom(EnviTriMesh);
  //
  // /* 4. BBs for Robot's Links */
  // SelfLinkGeoInfo SelfLinkGeoObj(*worldObj.robots[0], RMObject.EndEffectorLink2Pivotal, SelfCollisionFreeLink);      // Here SelfLinkGeoObj is instantiated with robot at zero configuration.
  //
  // int FileIndex = FileIndexFinder(false);
  // int TotalNumber = 100;
  //
  // string PlanningType = "RHP";
  // // string PlanningType = "OLP";
  //
  // /* 5. Internal Experimentation Loop */
  // while(FileIndex<=TotalNumber)
  // {
  //   string SpecificPath = ExperimentPath + "/" + std::to_string(FileIndex) + "/";
  //
  //   RobotWorld world;
  //   SimGUIBackend Backend(&world);
  //   WorldSimulation& Sim = Backend.sim;
  //
  //   string XMLFileStr = "../" + EnviName;
  //   const char* XMLFile = XMLFileStr.c_str();    // Here we must give abstract path to the file
  //   if(!Backend.LoadAndInitSim(XMLFile))
  //   {
  //     std::cerr<< EnviName << " file does not exist in that path!" << endl;
  //     return -1;
  //   }
  //   Robot SimRobot = *world.robots[0];
  //   RobotConfigLoader(SimRobot, SpecificPath, "InitConfig.config");
  //
  //   const std::string ContactStatusPath = SpecificPath + "ContactStatus.txt";
  //   std::vector<ContactStatusInfo> RobotContactInfo = ContactStatusInfoLoader(ContactStatusPath);
  //
  //   std::vector<double> InitRobotConfig(SimRobot.q);
  //   std::vector<double> InitRobotVelocity(SimRobot.q.size(), 0.0);
  //   std::vector<double> RobotConfigRef = InitRobotVelocity;
  //
  //   //  Given the optimized result to be the initial state
  //   Sim.world->robots[0]->UpdateConfig(Config(InitRobotConfig));
  //   Sim.world->robots[0]->dq = InitRobotVelocity;
  //
  //   Sim.controlSimulators[0].oderobot->SetConfig(Config(InitRobotConfig));
  //   Sim.controlSimulators[0].oderobot->SetVelocities(Config(InitRobotVelocity));
  //
  //   int PushRecovFlag = 0;
  //   int FailureFlag = 0;
  //   SpecificPath+= PlanningType + "/";
  //   FilePathManager(SpecificPath);
  //
  //   SimulationTest(Sim, NonlinearOptimizerInfo::RobotLinkInfo, RobotContactInfo, RMObject, TerrColGeom, SelfLinkGeoObj, SpecificPath, ForceMax, PushDuration, DetectionWait, PushRecovFlag, FailureFlag, PlanningType);
  //   if(FailureFlag)
  //   {
  //     PlanResWriter(SpecificPath, PushRecovFlag);
  //     FileIndex = FileIndexFinder(true);
  //     FileIndex++;
  //   }
  // }
  return 1;
}
