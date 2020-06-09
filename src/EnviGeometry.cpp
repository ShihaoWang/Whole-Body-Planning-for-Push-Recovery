#include "CommonHeader.h"
#include <vector>

static void SignedDistanceFieldWriter(const std::vector<double> & SDFVector, const std::vector<double> & SDFSpecs)
{
  // This function will write the computed SDFTensor and SDFSpecs into file
  FILE * SDFTensorFile = NULL;
  SDFTensorFile = fopen("./SDFs/SDFTensor.bin", "wb");
  fwrite(&SDFVector[0], sizeof(double), SDFVector.size(), SDFTensorFile);
  fclose(SDFTensorFile);

  FILE * SDFSpecsFile = NULL;
  SDFSpecsFile = fopen("./SDFs/SDFSpecs.bin", "wb");
  fwrite(&SDFSpecs[0], sizeof(double), SDFSpecs.size(), SDFSpecsFile);
  fclose(SDFSpecsFile);
}

SignedDistanceFieldInfo SignedDistanceFieldGene(const RobotWorld& WorldObj, const int& GridsNo)
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

  SignedDistanceFieldWriter(SDFVector, SDFSpecs);
  SignedDistanceFieldInfo SDFInfo(SDFTensor, SDFSpecs);

  return SDFInfo;
}

SignedDistanceFieldInfo SignedDistanceFieldLoader(const int GridsNo)
{
  // This function will read in the computed SDF_File into a Eigen::Tensor
  FILE* SDF_File = fopen("./SDFs/SDFTensor.bin", "rb");
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
  SDFSpecsFile = fopen("./SDFs/SDFSpecs.bin", "rb");
  fread(&SDFSpecs[0], sizeof(double), 13, SDFSpecsFile);
  fclose(SDFSpecsFile);

  SignedDistanceFieldInfo SDFInfo(SDFTensor, SDFSpecs);

  return SDFInfo;
}
