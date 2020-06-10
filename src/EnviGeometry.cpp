#include "CommonHeader.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>
#include <map>
#include <random>


static void SignedDistanceFieldWriter(const string & SDFPath, const std::vector<double> & SDFVector, const std::vector<double> & SDFSpecs)
{
  // This function will write the computed SDFTensor and SDFSpecs into file
  FILE * SDFTensorFile = NULL;
  const std::string SDFTensorFileNameStr = SDFPath + "SDFTensor.bin";
  const char* SDFTensorFileName = SDFTensorFileNameStr.c_str() ;
  SDFTensorFile = fopen(SDFTensorFileName, "wb");
  fwrite(&SDFVector[0], sizeof(double), SDFVector.size(), SDFTensorFile);
  fclose(SDFTensorFile);

  FILE * SDFSpecsFile = NULL;
  const std::string SDFSpecsFileNameStr = SDFPath + "SDFSpecs.bin";
  const char* SDFSpecsFileName = SDFSpecsFileNameStr.c_str();
  SDFSpecsFile = fopen(SDFSpecsFileName, "wb");
  fwrite(&SDFSpecs[0], sizeof(double), SDFSpecs.size(), SDFSpecsFile);
  fclose(SDFSpecsFile);
}

SignedDistanceFieldInfo SignedDistanceFieldGene(const string & SDFPath, const RobotWorld& WorldObj, const int& GridsNo)
{
  double resolution = 0.025;

  const int NumberOfTerrains = WorldObj.terrains.size();

  std::shared_ptr<Terrain> Terrain_ptr = std::make_shared<Terrain>(*WorldObj.terrains[0]);
  Meshing::TriMesh EnviTriMesh  = Terrain_ptr->geometry->AsTriangleMesh();

  // This step is used to merge the meshes into a single one.
  for (int i = 0; i < NumberOfTerrains-1; i++)
  {
    std::shared_ptr<Terrain> Terrain_ptr = std::make_shared<Terrain>(*WorldObj.terrains[i+1]);
    Meshing::TriMesh EnviTriMesh_i  = Terrain_ptr->geometry->AsTriangleMesh();
    EnviTriMesh.MergeWith(EnviTriMesh_i);
  }

  Meshing::VolumeGrid SDFGrid;
  CollisionMesh EnviTriMeshTopology(EnviTriMesh);
  EnviTriMeshTopology.InitCollisions();
  EnviTriMeshTopology.CalcTriNeighbors();
  MeshToImplicitSurface_FMM(EnviTriMeshTopology, SDFGrid, resolution);

  cout<<"FMM grid bounding box "<<SDFGrid.bb<<endl;

  // Now it is time to calculate SignedDistanceFieldInfo struct obj from SDFGrid

  // The estimated sizes of the environment
  double BB_x_min = SDFGrid.bb.bmin[0];
  double BB_x_max = SDFGrid.bb.bmax[0];

  double BB_y_min = SDFGrid.bb.bmin[1];
  double BB_y_max = SDFGrid.bb.bmax[1];

  double BB_z_min = SDFGrid.bb.bmin[2];
  double BB_z_max = SDFGrid.bb.bmax[2];

  double BB_x_length = BB_x_max - BB_x_min;
  double BB_y_length = BB_y_max - BB_y_min;
  double BB_z_length = BB_z_max - BB_z_min;

  double ExtCoeff = 0.0;

  // The estimated sizes of the environment
  double Envi_x_min = BB_x_min - ExtCoeff * BB_x_length;
  double Envi_x_max = BB_x_max + ExtCoeff * BB_x_length;

  double Envi_y_min = BB_y_min - ExtCoeff * BB_y_length;
  double Envi_y_max = BB_y_max + ExtCoeff * BB_y_length;

  double Envi_z_min = BB_z_min - ExtCoeff * BB_z_length;
  double Envi_z_max = BB_z_max + ExtCoeff * BB_z_length;

  double Envi_x_length = Envi_x_max - Envi_x_min;
  double Envi_y_length = Envi_y_max - Envi_y_min;
  double Envi_z_length = Envi_z_max - Envi_z_min;

  double Envi_x_unit = Envi_x_length/(1.0* GridsNo - 1.0);
  double Envi_y_unit = Envi_y_length/(1.0* GridsNo - 1.0);
  double Envi_z_unit = Envi_z_length/(1.0* GridsNo - 1.0);

  std::vector<double> Envi_x_coor(GridsNo), Envi_y_coor(GridsNo), Envi_z_coor(GridsNo);

  // Here Envi_x/y/z_coor is the actual grids from the given environment
  for (int i = 0; i < GridsNo; i++)
  {
    Envi_x_coor[i] = Envi_x_min + (1.0 * i) * Envi_x_unit;
    Envi_y_coor[i] = Envi_y_min + (1.0 * i) * Envi_y_unit;
    Envi_z_coor[i] = Envi_z_min + (1.0 * i) * Envi_z_unit;
  }

  // Generation of the SDFTensor structure
  Eigen::Tensor<double,3> SDFTensor(GridsNo, GridsNo, GridsNo);
  SDFTensor.setZero();

  Vector3 GridPoint;
  double GridPointDist, GridPoint_x, GridPoint_y, GridPoint_z;
  std::vector<double> SDFVector;
  SDFVector.reserve(GridsNo * GridsNo * GridsNo);
  for (int i = 0; i < GridsNo; i++)
  {
    GridPoint_x = Envi_x_coor[i];
    for (int j = 0; j < GridsNo; j++)
    {
      GridPoint_y = Envi_y_coor[j];
      for (int k = 0; k < GridsNo; k++)
      {
        GridPoint_z = Envi_z_coor[k];
        GridPoint.set(GridPoint_x, GridPoint_y, GridPoint_z);
        GridPointDist = SDFGrid.TrilinearInterpolate(GridPoint);
        SDFTensor(i,j,k) = GridPointDist;
        SDFVector.push_back(GridPointDist);
      }
    }
  }

  std::vector<double> SDFSpecs;
  SDFSpecs.push_back(Envi_x_min);             SDFSpecs.push_back(Envi_x_max);
  SDFSpecs.push_back(Envi_y_min);             SDFSpecs.push_back(Envi_y_max);
  SDFSpecs.push_back(Envi_z_min);             SDFSpecs.push_back(Envi_z_max);

  SDFSpecs.push_back(Envi_x_unit);            SDFSpecs.push_back(Envi_y_unit);            SDFSpecs.push_back(Envi_z_unit);
  SDFSpecs.push_back(Envi_x_length);          SDFSpecs.push_back(Envi_y_length);          SDFSpecs.push_back(Envi_z_length);
  SDFSpecs.push_back(GridsNo);

  SignedDistanceFieldWriter(SDFPath, SDFVector, SDFSpecs);
  SignedDistanceFieldInfo SDFInfo(SDFTensor, SDFSpecs);

  return SDFInfo;
}

SignedDistanceFieldInfo SignedDistanceFieldLoader(const string & SDFPath, const int GridsNo)
{
  // This function will read in the computed SDF_File into a Eigen::Tensor
  const std::string SDF_FileNameStr = SDFPath + "SDFTensor.bin";
  const char* SDF_FileName = SDF_FileNameStr.c_str();
  FILE* SDF_File = fopen(SDF_FileName, "rb");
  std::vector<double> SDFVector(GridsNo * GridsNo * GridsNo);
  fread(&SDFVector[0], sizeof(double), GridsNo * GridsNo * GridsNo, SDF_File);
  fclose(SDF_File);

  // The next job is to write to the Eigen::Tensor
  Eigen::Tensor<double,3> SDFTensor(GridsNo,GridsNo, GridsNo);
  int SDFIndex = 0;
  for (int i = 0; i < GridsNo; i++)
  {
    for (int j = 0; j < GridsNo; j++)
    {
      for (int k = 0; k < GridsNo; k++)
      {
        SDFTensor(i,j,k) = SDFVector[SDFIndex];
        SDFIndex +=1;
      }
    }
  }
  FILE * SDFSpecsFile = NULL;
  std::vector<double> SDFSpecs(13);
  const std::string SDFSpecsFileNameStr = SDFPath + "SDFSpecs.bin";
  const char* SDFSpecsFileName = SDFSpecsFileNameStr.c_str();
  SDFSpecsFile = fopen(SDFSpecsFileName, "rb");
  fread(&SDFSpecs[0], sizeof(double), 13, SDFSpecsFile);
  fclose(SDFSpecsFile);

  SignedDistanceFieldInfo SDFInfo(SDFTensor, SDFSpecs);

  return SDFInfo;
}

static std::vector<RMPoint> RMLayerGene(const double & r, const int & n){
  std::vector<RMPoint> RMLayer;
  RMLayer.reserve(n);
  double phi = M_PI * (3.0 - sqrt(5.0));      // golden angle in radians
  const double gr=(sqrt(5.0) + 1.0) / 2.0;    // golden ratio = 1.6180339887498948482
  const double ga=(2.0 - gr) * (2.0*M_PI);    // golden angle = 2.39996322972865332
  for (int i = 1; i <= n; i++) {
    const double lat = asin(-1.0 + 2.0 * double(i) / (n+1));
    const double lon = ga * i;

    const double x = cos(lon)*cos(lat);
    const double y = sin(lon)*cos(lat);
    const double z = sin(lat);

    double PointNorm = std::sqrt(x * x + y * y + z * z);
    double Ratio = PointNorm/r;
    Vector3 RMPointPos(x/Ratio, y/Ratio, z/Ratio);
    RMPoint RMPoint_i(r, RMPointPos);
    RMLayer.push_back(RMPoint_i);
  }
  return RMLayer;
}

// Here this function is used to generate a sufficiently dense reachability map for end effector(s)
ReachabilityMap ReachabilityMapGenerator(Robot & SimRobot, const std::vector<LinkInfo> & RobotLinkInfo, const std::vector<int> & TorsoLink){
  double MaxRadius = 0.9;
  int LayerNumber = 61;
  int PointNumberOnInner = 1;
  double LayerDiff = MaxRadius/(LayerNumber * 1.0);
  double MinRadius = LayerDiff;

  // Three uniform distributions
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> Xdis(-MaxRadius, MaxRadius);
  std::uniform_real_distribution<> Ydis(-MaxRadius, MaxRadius);
  std::uniform_real_distribution<> Zdis(-MaxRadius, MaxRadius);

  std::map<int, std::vector<RMPoint>> RMObj;

  int TotalPoint = 0;
  for (int i = 0; i < LayerNumber; i++)
  {
    double Radius = (1.0 * i + 1.0) * LayerDiff;
    int LayerPointLimit = (i + 1) * (i + 1) * PointNumberOnInner;
    int LayerPointNumber = 0;
    std::vector<RMPoint> RMLayer;
    RMLayer = RMLayerGene(Radius, LayerPointLimit);
    RMObj[i] = RMLayer;
  }
  ReachabilityMap RMObject(RMObj);
  RMObject.ReachabilityMapPara(MaxRadius, LayerNumber, PointNumberOnInner, LayerDiff, MinRadius);

  // The next step is to estimate the end effector radius.
  int DOF = SimRobot.q.size();
  std::vector<double> MeasureConfiguration(DOF);
  for (int i = 0; i < DOF; i++){
    MeasureConfiguration[i] = 0.0;
  }

  // These two values are modified to ensure the arms are straight.
  MeasureConfiguration[23] = SimRobot.qMin[23];
  MeasureConfiguration[30] = SimRobot.qMax[30];
  SimRobot.UpdateConfig(Config(MeasureConfiguration));
  SimRobot.UpdateGeometry();
  std::vector<double> EndEffectorRadius(RobotLinkInfo.size());
  std::vector<double> EndEffectorCollisionRadius(RobotLinkInfo.size());
  std::vector<int> EndEffectorLinkIndex(RobotLinkInfo.size());
  std::vector<int> EndEffectorPivotalIndex(RobotLinkInfo.size());
  std::map<int, std::vector<int>> EndEffectorLink2Pivotal;

  std::map<int, int> EndEffectorIndices;

  for (int i = 0; i < RobotLinkInfo.size(); i++){
    int ParentIndex = -1;
    int CurrentIndex = RobotLinkInfo[i].LinkIndex;
    EndEffectorIndices[CurrentIndex] = 1;
    EndEffectorLinkIndex[i] = RobotLinkInfo[i].LinkIndex;
    std::vector<int> EndEffectorLink2PivotalIndex;
    EndEffectorLink2PivotalIndex.push_back(CurrentIndex);
    while(std::find(TorsoLink.begin(), TorsoLink.end(), ParentIndex)==TorsoLink.end()) {
      ParentIndex = SimRobot.parents[CurrentIndex];
      CurrentIndex = ParentIndex;
      EndEffectorLink2PivotalIndex.push_back(CurrentIndex);
    }
    EndEffectorLink2PivotalIndex.pop_back();
    Vector3 PivotalPos, EndPos, PivotalRef(0.0, 0.0, 0.0);
    SimRobot.GetWorldPosition(PivotalRef, EndEffectorLink2PivotalIndex[EndEffectorLink2PivotalIndex.size()-1], PivotalPos);
    SimRobot.GetWorldPosition(RobotLinkInfo[i].AvgLocalContact, RobotLinkInfo[i].LinkIndex, EndPos);
    EndEffectorPivotalIndex[i] = EndEffectorLink2PivotalIndex[EndEffectorLink2PivotalIndex.size()-1];
    EndEffectorLink2Pivotal[i] = EndEffectorLink2PivotalIndex;

    std::vector<double> CollisionRadius(RobotLinkInfo[i].LocalContacts.size());
    for (int j = 0; j < RobotLinkInfo[i].LocalContacts.size(); j++){
      Vector3 Center2Edge = RobotLinkInfo[i].LocalContacts[j] - RobotLinkInfo[i].AvgLocalContact;
      CollisionRadius[j] = sqrt(Center2Edge.x * Center2Edge.x + Center2Edge.y * Center2Edge.y + Center2Edge.z * Center2Edge.z);
    }
    EndEffectorCollisionRadius[i] = *std::max_element(CollisionRadius.begin(), CollisionRadius.end());
    Vector3 Pivotal2End = PivotalPos - EndPos;
    double Pivotal2EndRadius = sqrt(Pivotal2End.x * Pivotal2End.x + Pivotal2End.y * Pivotal2End.y + Pivotal2End.z * Pivotal2End.z);
    EndEffectorRadius[i] = Pivotal2EndRadius;
  }
  RMObject.EndEffectorRadius = EndEffectorRadius;
  RMObject.EndEffectorLinkIndex = EndEffectorLinkIndex;
  RMObject.EndEffectorCollisionRadius = EndEffectorCollisionRadius;
  RMObject.EndEffectorPivotalIndex = EndEffectorPivotalIndex;
  RMObject.EndEffectorLink2Pivotal = EndEffectorLink2Pivotal;
  RMObject.EndEffectorIndices = EndEffectorIndices;

  RMObject.TotalPoint = TotalPoint;
  return RMObject;
}
