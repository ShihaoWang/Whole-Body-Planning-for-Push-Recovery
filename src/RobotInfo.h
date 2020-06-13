#ifndef ROBOTINFO_H
#define ROBOTINFO_H
#include <KrisLibrary/robotics/RobotDynamics3D.h>
#include <KrisLibrary/robotics/NewtonEuler.h>
#include <Interface/SimulationGUI.h>
#include <unsupported/Eigen/CXX11/Tensor>
#include <KrisLibrary/meshing/TriMeshTopology.h>
#include <KrisLibrary/math3d/geometry3d.h>
#include <KrisLibrary/geometry/CollisionMesh.h>
#include <KrisLibrary/geometry/PQP/src/PQP.h>
#include "Modeling/Paths.h"

struct LinkInfo {
  LinkInfo(){LinkIndex = -1;}
  LinkInfo(const int &link_index){ LinkIndex = link_index; }
  void AddLocalConact(const Vector3& local_contact){ LocalContacts.push_back(local_contact); }
  void AvgContactUpdate(){
    switch (LocalContacts.size()){
      case 0: throw std::invalid_argument( "LocalContacts should have been initialized!" );
      break;
      default:{
        Vector3 SumLocalContacts(0.0, 0.0, 0.0);
        for (int i = 0; i < LocalContacts.size(); i++){
          SumLocalContacts.x = SumLocalContacts.x + LocalContacts[i].x;
          SumLocalContacts.y = SumLocalContacts.y + LocalContacts[i].y;
          SumLocalContacts.z = SumLocalContacts.z + LocalContacts[i].z;
        }
        AvgLocalContact.x = SumLocalContacts.x/LocalContacts.size();
        AvgLocalContact.y = SumLocalContacts.y/LocalContacts.size();
        AvgLocalContact.z = SumLocalContacts.z/LocalContacts.size();
      }
      break;
    }
  }
  int LinkIndex;
  std::vector<Vector3> LocalContacts;
  std::vector<Vector3> ContactPositions;
  std::vector<Vector3> ContactVelocities;
  std::vector<double> ContactDists;
  Vector3 AvgLocalContact;
};

struct ContactStatusInfo{
  ContactStatusInfo(){ LinkIndex = -1;}
  ContactStatusInfo(const int & link_index){ LinkIndex = link_index; }
  void AddLocalConactStatus(const int & _contactstatus){ LocalContactStatus.push_back(_contactstatus); }
  void StatusSwitch(const int & Val){
    for(int i = 0; i<LocalContactStatus.size(); i++){
      LocalContactStatus[i] = Val;
    }
  }
  int LinkIndex;
  std::vector<int> LocalContactStatus;
};

struct TerrainInfo{
  TerrainInfo(){num_tris = -1;}
  int num_tris;
  std::vector<Tri> Tris;
  std::vector<Vector3> TriNormals;
  std::vector<int> Indices;
};

struct SignedDistanceFieldInfo{
  SignedDistanceFieldInfo(){
    Envi_x_min = 0;           Envi_x_max = 0;
    Envi_y_min = 0;           Envi_y_max = 0;
    Envi_z_min = 0;           Envi_z_max = 0;
    Envi_x_unit = 0;          Envi_y_unit = 0;          Envi_z_unit = 0;
    Envi_x_length = 0;        Envi_y_length = 0;        Envi_z_length = 0;
    GridNo = 0;
  }
  SignedDistanceFieldInfo(const Eigen::Tensor<double, 3>& _SDFTensor, const std::vector<double> &_SDFSpecs){
    SDFTensor = _SDFTensor;
    Envi_x_min = _SDFSpecs[0];          Envi_x_max = _SDFSpecs[1];
    Envi_y_min = _SDFSpecs[2];          Envi_y_max = _SDFSpecs[3];
    Envi_z_min = _SDFSpecs[4];          Envi_z_max = _SDFSpecs[5];
    Envi_x_unit = _SDFSpecs[6];         Envi_y_unit = _SDFSpecs[7];         Envi_z_unit = _SDFSpecs[8];
    Envi_x_length = _SDFSpecs[9];       Envi_y_length = _SDFSpecs[10];      Envi_z_length = _SDFSpecs[11];
    GridNo = (int)_SDFSpecs[12];
  }
  double SignedDistance(const Vector3 &Point) const{
    double x_FloatIndex = (Point.x - Envi_x_min)/Envi_x_unit * 1.0;
    double y_FloatIndex = (Point.y - Envi_y_min)/Envi_y_unit * 1.0;
    double z_FloatIndex = (Point.z - Envi_z_min)/Envi_z_unit * 1.0;

    int x_leftindex = std::floor(x_FloatIndex);
    int y_leftindex = std::floor(y_FloatIndex);
    int z_leftindex = std::floor(z_FloatIndex);

    if(x_leftindex<0){
      x_leftindex = 0;
    } else {
      if(x_leftindex>GridNo-2) {
        x_leftindex = GridNo-2;
      }
    }

    if(y_leftindex<0) {
      y_leftindex = 0;
    } else {
      if(y_leftindex>GridNo-2) {
        y_leftindex = GridNo-2;
      }
    }

    if(z_leftindex<0) {
      z_leftindex = 0;
    } else {
      if(z_leftindex>GridNo-2) {
        z_leftindex = GridNo-2;
      }
    }

    int x_rightindex = x_leftindex + 1;
    int y_rightindex = y_leftindex + 1;
    int z_rightindex = z_leftindex + 1;

    double valA = SDFTensor(x_leftindex, y_leftindex, z_leftindex);
    double valB = SDFTensor(x_rightindex, y_leftindex, z_leftindex);
    double valC = SDFTensor(x_rightindex, y_rightindex, z_leftindex);
    double valD = SDFTensor(x_leftindex, y_rightindex, z_leftindex);

    double valE = SDFTensor(x_leftindex, y_leftindex, z_rightindex);
    double valF = SDFTensor(x_rightindex, y_leftindex, z_rightindex);
    double valG = SDFTensor(x_rightindex, y_rightindex, z_rightindex);
    double valH = SDFTensor(x_leftindex, y_rightindex, z_rightindex);

    // Since this is a tri-linear interpolation, there are three ways to do the interpolation

    // Type1:
    // Along x-direction
    double valMAB = (x_FloatIndex - x_leftindex*1.0) * (valB - valA) + valA;
    double valMDC = (x_FloatIndex - x_leftindex*1.0) * (valC - valD) + valD;
    double valMEF = (x_FloatIndex - x_leftindex*1.0) * (valF - valE) + valE;
    double valMHG = (x_FloatIndex - x_leftindex*1.0) * (valG - valH) + valH;

    // Along y-drection
    double valMABDC = (y_FloatIndex - y_leftindex*1.0) * (valMDC - valMAB) + valMAB;
    double valMEFGH = (y_FloatIndex - y_leftindex*1.0) * (valMHG - valMEF) + valMEF;

    // Along z-direction
    double valMABCDEFHG = (z_FloatIndex - z_leftindex*1.0) * (valMEFGH - valMABDC) + valMABDC;

    // // Type2:
    // // Along y-drection
    // double valMAD = (y_FloatIndex - y_leftindex*1.0) * (valD - valA) + valA;
    // double valMBC = (y_FloatIndex - y_leftindex*1.0) * (valC - valB) + valB;
    // double valMEH = (y_FloatIndex - y_leftindex*1.0) * (valH - valE) + valE;
    // double valMFG = (y_FloatIndex - y_leftindex*1.0) * (valG - valF) + valF;
    //
    // //Along x-direction
    // double valMADBC = (x_FloatIndex - x_leftindex * 1.0) * (valMBC - valMAD) + valMAD;
    // double valMEHFG = (x_FloatIndex - x_leftindex * 1.0) * (valMFG - valMEH) + valMEH;
    //
    // //Along z-direction
    // double valMADBCEHFG = (z_FloatIndex - z_leftindex*1.0) * (valMEHFG - valMADBC) + valMADBC;

    return valMABCDEFHG;
  }
  Vector3 SignedDistanceNormal(const Vector3 &Point) const {
    // This function is used to calculate the (1 x 3) Jacobian matrix given the current Position
    // This function is used to compute the distance from a 3D point to the environment terrain
    // The first job is to figure out the nearest neighbours of the Points

    double x_FloatIndex = (Point.x - Envi_x_min)/Envi_x_unit * 1.0;
    double y_FloatIndex = (Point.y - Envi_y_min)/Envi_y_unit * 1.0;
    double z_FloatIndex = (Point.z - Envi_z_min)/Envi_z_unit * 1.0;

    int x_leftindex = std::floor(x_FloatIndex);
    int y_leftindex = std::floor(y_FloatIndex);
    int z_leftindex = std::floor(z_FloatIndex);

    if(x_leftindex<0) {
      x_leftindex = 0;
    } else {
      if(x_leftindex>GridNo-2) {
        x_leftindex = GridNo-2;
      }
    }

    if(y_leftindex<0) {
      y_leftindex = 0;
    } else {
      if(y_leftindex>GridNo-2) {
        y_leftindex = GridNo-2;
      }
    }

    if(z_leftindex<0) {
      z_leftindex = 0;
    } else {
      if(z_leftindex>GridNo-2) {
        z_leftindex = GridNo-2;
      }
    }

    int x_rightindex = x_leftindex + 1;
    int y_rightindex = y_leftindex + 1;
    int z_rightindex = z_leftindex + 1;

    double valA = SDFTensor(x_leftindex, y_leftindex, z_leftindex);
    double valB = SDFTensor(x_rightindex, y_leftindex, z_leftindex);
    double valC = SDFTensor(x_rightindex, y_rightindex, z_leftindex);
    double valD = SDFTensor(x_leftindex, y_rightindex, z_leftindex);

    double valE = SDFTensor(x_leftindex, y_leftindex, z_rightindex);
    double valF = SDFTensor(x_rightindex, y_leftindex, z_rightindex);
    double valG = SDFTensor(x_rightindex, y_rightindex, z_rightindex);
    double valH = SDFTensor(x_leftindex, y_rightindex, z_rightindex);

    // Since this is a tri-linear interpolation, there are three ways to do the interpolation

    /*
      Jacobian matrix in the x-direction
    */

    // Cut the point with a plane orthgonal to z axis
    double valMAE =  (z_FloatIndex - z_leftindex*1.0) * (valE - valA) + valA;
    double valMBF =  (z_FloatIndex - z_leftindex*1.0) * (valF - valB) + valB;
    double valMDH =  (z_FloatIndex - z_leftindex*1.0) * (valH - valD) + valD;
    double valMCG =  (z_FloatIndex - z_leftindex*1.0) * (valG - valC) + valC;
    // Cut the point with a plane orthgonal to y axis
    double valMAEDH = (y_FloatIndex - y_leftindex*1.0) * (valMDH - valMAE) + valMAE;
    double valMBFCG = (y_FloatIndex - y_leftindex*1.0) * (valMCG - valMBF) + valMBF;
    // The values at the edge give the jacobian to x
    double JacDistTo_x = (valMBFCG - valMAEDH)/Envi_x_unit;

    /*
      Jacobian matrix in the y-direction
    */
    double valMABEF = (x_FloatIndex - x_leftindex*1.0) * (valMBF - valMAE) + valMAE;
    double valMDHCG = (x_FloatIndex - x_leftindex*1.0) * (valMCG - valMDH) + valMDH;

    double JacDistTo_y = (valMBFCG - valMAEDH)/Envi_y_unit;

    /*
      Jacobian matrix in the z-direction
    */
    // Cut the point with a plane orthgonal to x axis
    double valMAB = (x_FloatIndex - x_leftindex*1.0) * (valB - valA) + valA;
    double valMDC = (x_FloatIndex - x_leftindex*1.0) * (valC - valD) + valD;
    double valMEF = (x_FloatIndex - x_leftindex*1.0) * (valF - valE) + valE;
    double valMHG = (x_FloatIndex - x_leftindex*1.0) * (valG - valH) + valH;
    // Cut the point with a plane orthgonal to y axis
    double valMABDC = (y_FloatIndex - y_leftindex*1.0) * (valMDC - valMAB) + valMAB;
    double valMEFHG = (y_FloatIndex - y_leftindex*1.0) * (valMHG - valMEF) + valMHG;
    // The values at the edge give the jacobian to z
    double JacDistTo_z = (valMEFHG - valMABDC)/Envi_z_unit;

    Vector3 RawNormal(JacDistTo_x, JacDistTo_y, JacDistTo_z);
    RawNormal.setNormalized(RawNormal);
    return RawNormal;
  }
  Eigen::Tensor<double, 3> SDFTensor;
  double Envi_x_min, Envi_x_max;
  double Envi_y_min, Envi_y_max;
  double Envi_z_min, Envi_z_max;
  double Envi_x_unit, Envi_y_unit, Envi_z_unit;
  double Envi_x_length, Envi_y_length, Envi_z_length;
  int GridNo;
};

struct RMPoint {
  // This struct is used to save the information of a point in Reachability Map.
  RMPoint();
  RMPoint(const double & r, const Vector3 & Pos) {
    // Constructor
    Radius = r;
    Position = Pos;
    Direction = Position;
    double PositionLength = std::sqrt(Position.x * Position.x + Position.y * Position.y + Position.z * Position.z);
    Direction.x = Direction.x/PositionLength;
    Direction.y = Direction.y/PositionLength;
    Direction.z = Direction.z/PositionLength;
  }
  double Radius;
  Vector3 Position;
  Vector3 Direction;
};

struct ReachabilityMap {
  // This struct is used to save the reachability map
  ReachabilityMap();
  ReachabilityMap(const std::map<int, std::vector<RMPoint>> & _RMLayers):RMLayers(_RMLayers){};
  void ReachabilityMapPara(const double _MaxRadius, const int _LayerNumber, const int _PointNumberOnInner, const double _LayerDiff, const double _MinRadius){
    MaxRadius = _MaxRadius;
    LayerNumber = _LayerNumber;
    PointNumberOnInner = _PointNumberOnInner;
    LayerDiff = _LayerDiff;
    MinRadius = _MinRadius;
  }
  std::vector<Vector3> IdealReachablePointsFinder(const Robot & SimRobot, const int & LinkInfoIndex){
    std::vector<Vector3> ReachablePoints;
    ReachablePoints.reserve(TotalPoint);
    double PivotalLinkIndex = EndEffectorPivotalIndex[LinkInfoIndex];
    Vector3 RefPoint, ZeroPos(0.0, 0.0, 0.0);
    SimRobot.GetWorldPosition(ZeroPos, PivotalLinkIndex, RefPoint);
    for (int i = 0; i < LayerNumber; i++){
      std::vector<RMPoint> RMLayer_i = RMLayers[i];
      for (int j = 0; j < RMLayer_i.size(); j++){
        Vector3 RMPointPos = RMLayer_i[j].Position + RefPoint;
        ReachablePoints.push_back(RMPointPos);
      }
    }
    return ReachablePoints;
  }
  std::vector<Vector3> ReachablePointsFinder(const Robot & SimRobot, const int & LinkInfoIndex, SignedDistanceFieldInfo & SDFInfo, const Vector3 & COMVel){
    double Radius = EndEffectorRadius[LinkInfoIndex];
    double PivotalLinkIndex = EndEffectorPivotalIndex[LinkInfoIndex];

    Vector3 RefPoint, ZeroPos(0.0, 0.0, 0.0);
    SimRobot.GetWorldPosition(ZeroPos, PivotalLinkIndex, RefPoint);
    int ReachablePointNo = 0;
    return ReachablePointsGene(RefPoint, Radius, SDFInfo, ReachablePointNo, COMVel);
  }
  std::vector<Vector3> ReachablePointsFinder(const Robot & SimRobot, const int & LinkInfoIndex, SignedDistanceFieldInfo & SDFInfo){
    double Radius = EndEffectorRadius[LinkInfoIndex];
    double PivotalLinkIndex = EndEffectorPivotalIndex[LinkInfoIndex];

    Vector3 RefPoint, ZeroPos(0.0, 0.0, 0.0);
    SimRobot.GetWorldPosition(ZeroPos, PivotalLinkIndex, RefPoint);
    Vector3 COMVel(0.0, 0.0, 0.0);
    int ReachablePointNo = 0;
    return ReachablePointsGene(RefPoint, Radius, SDFInfo, ReachablePointNo, COMVel);
  }
  std::vector<Vector3> ReachablePointsGene(const Vector3 & RefPoint, const double & Radius, SignedDistanceFieldInfo & SDFInfo, int & ReachablePointNo, const Vector3 & COMVel){
    // Here MinRadius is used for comparison while Radius is the maximum allowed radius of robot's end effector.
    // This function is used to get presumably active ReachablePointsGene() from all sampled points.
    const double DisTol = 0.01;        // 1cm as a signed distance tolerance.
    std::vector<Vector3> ReachablePoints;
    ReachablePoints.reserve(TotalPoint);
    ReachablePointNo = 0;
    int LayerIndex = 0;
    double LayerRadius = MinRadius;
    if(Radius>MaxRadius){
      LayerIndex = LayerNumber-1;
    } else{
      while (LayerRadius<Radius){
        LayerRadius+=LayerDiff;
        LayerIndex++;
      }
      LayerRadius-=LayerDiff;
    }

    for (int i = 0; i < LayerIndex; i++){
      std::vector<RMPoint> RMLayer_i = RMLayers[i];
      for (int j = 0; j < RMLayer_i.size(); j++){
        Vector3 RMPointPos = RMLayer_i[j].Position + RefPoint;
        double CurrentDist = SDFInfo.SignedDistance(RMPointPos);
        if(CurrentDist*CurrentDist<DisTol*DisTol){
          ReachablePoints.push_back(RMPointPos);
          ReachablePointNo++;
        }
      }
    }
    return ReachablePoints;
  }
  std::vector<Vector3> ContactFreePointsFinder(const double & radius, const std::vector<Vector3> & ReachablePoints,const std::vector<std::pair<Vector3, double>> & ContactFreeInfo){
    // This function can only be called after ReachablePointsFinder() to reduce the extra point further.
    std::vector<Vector3> ContactFreePoints;
    ContactFreePoints.reserve(ReachablePoints.size());
    int ContactFreeNo = 0;
    for (int i = 0; i < ReachablePoints.size(); i++){
      Vector3 ReachablePoint = ReachablePoints[i];
      bool ContactFreeFlag = true;
      int ContactFreeInfoIndex = 0;
      while (ContactFreeInfoIndex<ContactFreeInfo.size()){
        Vector3 RefPoint = ContactFreeInfo[ContactFreeInfoIndex].first;
        double Radius = ContactFreeInfo[ContactFreeInfoIndex].second;
        Vector3 PosDiff = ReachablePoint - RefPoint;
        double PosDiffDis = sqrt(PosDiff.x * PosDiff.x + PosDiff.y * PosDiff.y + PosDiff.z * PosDiff.z);
        if(PosDiffDis<=(Radius + radius)){
          ContactFreeFlag = false;
          break;
        }
        ContactFreeInfoIndex++;
      }
      if(ContactFreeFlag){
        ContactFreePoints.push_back(ReachablePoints[i]);
        ContactFreeNo++;
      }
    }
    return ContactFreePoints;
  };
  std::map<int, std::vector<RMPoint>> RMLayers;
  std::vector<double> EndEffectorCollisionRadius;
  std::vector<double> EndEffectorRadius;
  std::vector<int> EndEffectorLinkIndex;              //
  std::vector<int> EndEffectorPivotalIndex;
  std::map<int, std::vector<int>> EndEffectorLink2Pivotal;       // This map saves intermediate joint from End Effector Joint to Pivotal Joint.
  std::map<int, int> EndEffectorIndices;
  double MaxRadius;
  int LayerNumber;
  int PointNumberOnInner;
  double LayerDiff;
  double MinRadius;
  int TotalPoint;
};

struct SelfLinkGeoInfo{
  SelfLinkGeoInfo();
  SelfLinkGeoInfo(const Robot & SimRobot, const std::map<int, std::vector<int>> & EndEffectorLink2Pivotal, const std::vector<int> & SelfCollisionFreeLink){
    for (int i = 5; i < SimRobot.q.size(); i++){
      AABB3D AABB3D_i = SimRobot.geometry[i]->GetAABBTight();
      Frame3D LinkTransforms_i = SimRobot.links[i].T_World;
      LinkBBs.push_back(AABB3D_i);
      LinkTransforms.push_back(LinkTransforms_i);
    }

    for (int i = 0; i < EndEffectorLink2Pivotal.size(); i++){
      std::vector<int> EndLink2PivIndices = EndEffectorLink2Pivotal.at(i);
      for (int j = 5; j < SimRobot.q.size(); j++){
        if(std::find(EndLink2PivIndices.begin(), EndLink2PivIndices.end(), j) == EndLink2PivIndices.end()){
          if(std::find(SelfCollisionFreeLink.begin(), SelfCollisionFreeLink.end(), j) == SelfCollisionFreeLink.end()){
            SelfCollisionLinkMap[i].push_back(j-5);
          }
        }
      }
    }
  };
  void LinkBBsUpdate(const Robot& SimRobot){
    for (int i = 5; i < SimRobot.q.size(); i++){
      AABB3D AABB3D_i = SimRobot.geometry[i]->GetAABBTight();
      LinkBBs[i-5] = AABB3D_i;
      Frame3D LinkTransforms_i = SimRobot.links[i].T_World;
      LinkTransforms[i-5] = LinkTransforms_i;
    }
  }
  double SingleLinkDist(const int & LinkCountIndex, const Vector3 & GlobalPoint){
    return LinkBBs[LinkCountIndex].signedDistance(GlobalPoint);
  }

  void SingleLinkDistNGrad(const int & LinkCountIndex, const Vector3 & GlobalPoint, double & Dist, Vector3 & DistGrad){
    const int GridNo = 100;
    double dx = LinkBBs[LinkCountIndex].size().x/(1.0 * GridNo);
    double dy = LinkBBs[LinkCountIndex].size().y/(1.0 * GridNo);
    double dz = LinkBBs[LinkCountIndex].size().z/(1.0 * GridNo);

    Dist = LinkBBs[LinkCountIndex].signedDistance(GlobalPoint);

    Vector3 GlobalPointx = GlobalPoint;
    GlobalPointx.x += dx;
    double Distx = LinkBBs[LinkCountIndex].signedDistance(GlobalPointx);

    Vector3 GlobalPointy = GlobalPoint;
    GlobalPointy.y += dy;
    double Disty = LinkBBs[LinkCountIndex].signedDistance(GlobalPointy);

    Vector3 GlobalPointz = GlobalPoint;
    GlobalPointz.z += dz;
    double Distz = LinkBBs[LinkCountIndex].signedDistance(GlobalPointz);

    DistGrad.x = (Distx - Dist)/dx;
    DistGrad.y = (Disty - Dist)/dy;
    DistGrad.z = (Distz - Dist)/dz;
    DistGrad.getNormalized(DistGrad);
  }

  std::vector<Vector3> BBVertices(const int & BBIndex){
    std::vector<Vector3> Vertices(8);
    // This function calculates the vertices for current bounding box.
    AABB3D CurBB = LinkBBs[BBIndex];
    Vector3 CurBBSize = CurBB.size();
    RigidTransform CurBBtoWorld = LinkTransforms[BBIndex];

    Vector3 bmin = CurBB.bmin;
    Vector3 bmax = CurBB.bmax;

    Vector3 bmin2bmax = bmax - bmin;

    Vector3 xAxis(CurBBtoWorld.R.data[0][0], CurBBtoWorld.R.data[0][1], CurBBtoWorld.R.data[0][2]);
    Vector3 yAxis(CurBBtoWorld.R.data[1][0], CurBBtoWorld.R.data[1][1], CurBBtoWorld.R.data[1][2]);
    Vector3 zAxis(CurBBtoWorld.R.data[2][0], CurBBtoWorld.R.data[2][1], CurBBtoWorld.R.data[2][2]);

    double xbmin2bmax = xAxis.dot(bmin2bmax);
    double ybmin2bmax = yAxis.dot(bmin2bmax);
    double zbmin2bmax = zAxis.dot(bmin2bmax);

    // Now we can have eight points explicitly.
    Vertices[0] = bmin;
    Vertices[1] = bmin + xbmin2bmax * xAxis;
    Vertices[2] = bmin + zbmin2bmax * zAxis;
    Vertices[3] = bmin + xbmin2bmax * xAxis + zbmin2bmax * zAxis;

    Vertices[4] = bmin + ybmin2bmax * yAxis;
    Vertices[5] = bmin + ybmin2bmax * yAxis + xbmin2bmax * xAxis;
    Vertices[6] = bmin + ybmin2bmax * yAxis + zbmin2bmax * zAxis;
    Vertices[7] = bmin + ybmin2bmax * yAxis + xbmin2bmax * xAxis + zbmin2bmax * zAxis;
    return Vertices;
  }

  void Vector3Writer(const std::vector<Vector3> & ContactPoints, const std::string &ContactPointFileName){
    if(ContactPoints.size()==0){
      std::cerr<< " ContactPoints has zero element!\n" << endl;
    }
    int NumberOfContactPoints = ContactPoints.size();
    std::vector<double> FlatContactPoints(3 * NumberOfContactPoints);
    int FlatContactPointIndex = 0;
    for (int i = 0; i < NumberOfContactPoints; i++){
      FlatContactPoints[FlatContactPointIndex] = ContactPoints[i].x;
      FlatContactPointIndex++;
      FlatContactPoints[FlatContactPointIndex] = ContactPoints[i].y;
      FlatContactPointIndex++;
      FlatContactPoints[FlatContactPointIndex] = ContactPoints[i].z;
      FlatContactPointIndex++;
    }
    FILE * FlatContactPointsFile = NULL;
    string ContactPointFile = ContactPointFileName + ".bin";
    const char *ContactPointFile_Name = ContactPointFile.c_str();
    FlatContactPointsFile = fopen(ContactPointFile_Name, "wb");
    fwrite(&FlatContactPoints[0], sizeof(double), FlatContactPoints.size(), FlatContactPointsFile);
    fclose(FlatContactPointsFile);
    return;
  }

  double SelfCollisionDist(const int & LinkIndex, const Vector3 & GlobalPoint){
    const int ActLinkNo = SelfCollisionLinkMap[LinkIndex].size();
    std::vector<double> DistVec;
    DistVec.reserve(ActLinkNo);
    for (int i = 0; i < ActLinkNo; i++){
      int SelfLinkIndex = SelfCollisionLinkMap[LinkIndex][i];
      double Dist_i = LinkBBs[SelfLinkIndex].signedDistance(GlobalPoint);
      DistVec.push_back(Dist_i);
    }
    double Dist = *std::min_element(DistVec.begin(), DistVec.end());
    return Dist;
  }

  void SelfCollisionDistNGrad(const int & LinkIndex, const Vector3 & GlobalPoint, double & Dist, Vector3 & Grad){
    // This function is used to calculate robot's self-collision distance given a point
    const int ActLinkNo = SelfCollisionLinkMap[LinkIndex].size();
    std::vector<double> DistVec;
    DistVec.reserve(ActLinkNo);
    std::vector<Vector3> GradVec;
    GradVec.reserve(ActLinkNo);
    std::vector<double> DistWeights;
    DistWeights.reserve(ActLinkNo);
    double Dist_i;
    Vector3 Grad_i;
    for (int i = 0; i < ActLinkNo; i++){
      int SelfLinkIndex = SelfCollisionLinkMap[LinkIndex][i];
      SingleLinkDistNGrad(SelfLinkIndex, GlobalPoint, Dist_i, Grad_i);
      DistVec.push_back(Dist_i);
      GradVec.push_back(Grad_i);
    }
    Dist = *std::min_element(DistVec.begin(), DistVec.end());
    double Scale = abs(Dist);
    for (int i = 0; i < ActLinkNo; i++){
      DistWeights.push_back(exp(-1.0 * DistVec[i]/Scale));
    }
    // Set its value to be zero!
    Grad.x = 0.0;
    Grad.y = 0.0;
    Grad.z = 0.0;
    for (int i = 0; i < ActLinkNo; i++){
      Grad+=DistWeights[i] * GradVec[i];
    }
    Grad.getNormalized(Grad);
  }
  std::vector<AABB3D> LinkBBs;
  std::vector<RigidTransform> LinkTransforms;
  std::map<int, std::vector<int>> SelfCollisionLinkMap;       // This map saves intermediate joint from End Effector Joint to Pivotal Joint.
};

struct SimPara{
  SimPara();
  SimPara(const double & _ForceMax,
          const double & _PushDuration,
          const double & _DetectionWait,
          const double & _TimeStep,
          const double & _InitDuration,
          const double & _TotalDuration,
          const double & _PhaseRatio):    ForceMax(_ForceMax),
                                          PushDuration(_PushDuration),
                                          DetectionWait(_DetectionWait),
                                          TimeStep(_TimeStep),
                                          InitDuration(_InitDuration),
                                          TotalDuration(_TotalDuration),
                                          PhaseRatio(_PhaseRatio){}
  void CurrentCasePathUpdate(const string _CurrentCasePath){
    CurrentCasePath = _CurrentCasePath;

    string fEdgeAFile = CurrentCasePath + "EdgeATraj.txt";
    // const char *fEdgeAFile_Name = fEdgeAFile.c_str();
    string fEdgeBFile = CurrentCasePath + "EdgeBTraj.txt";
    // const char *fEdgeBFile_Name = fEdgeBFile.c_str();
    string fEdgeCOMFile = CurrentCasePath + "EdgeCOMTraj.txt";
    // const char *fEdgeCOMFile_Name = fEdgeCOMFile.c_str();
    string fEdgexTrajFile = CurrentCasePath + "EdgexTraj.txt";
    // const char *fEdgexTrajFile_Name = fEdgexTrajFile.c_str();
    string fEdgeyTrajFile = CurrentCasePath + "EdgeyTraj.txt";
    // const char *fEdgeyTrajFile_Name = fEdgeyTrajFile.c_str();
    string fEdgezTrajFile = CurrentCasePath + "EdgezTraj.txt";
    // const char *fEdgezTrajFile_Name = fEdgezTrajFile.c_str();
    string fVertexTrajFile = CurrentCasePath + "EdgeVertexTraj.txt";

    EdgeFileNames.push_back(fEdgeAFile);
    EdgeFileNames.push_back(fEdgeBFile);
    EdgeFileNames.push_back(fEdgeCOMFile);
    EdgeFileNames.push_back(fEdgexTrajFile);
    EdgeFileNames.push_back(fEdgeyTrajFile);
    EdgeFileNames.push_back(fEdgezTrajFile);
    EdgeFileNames.push_back(fVertexTrajFile);

    FailureStateTrajStr =  CurrentCasePath + "FailureStateTraj.path";
    // const char *FailureStateTrajStr_Name = FailureStateTrajStr.c_str();
    CtrlStateTrajStr    =  CurrentCasePath + "CtrlStateTraj.path";
    // const char *CtrlStateTrajStr_Name = CtrlStateTrajStr.c_str();
    PlanStateTrajFileStr = CurrentCasePath + "PlanStateTraj.path";
    // const char *PlanStateTrajStr_Name = PlanStateTrajFileStr.c_str();
  }
  double ForceMax;
  double PushDuration;
  double DetectionWait;
  double TimeStep;
  double InitDuration;
  double TotalDuration;
  double PhaseRatio;            // This ratio determines the boundary between acceleration and deceleration.
  std::string CurrentCasePath;
  std::vector<string >EdgeFileNames;
  string FailureStateTrajStr, CtrlStateTrajStr, PlanStateTrajFileStr;
};

struct ControlReferenceInfo{
  ControlReferenceInfo(){
    ReadyFlag = false;
    ControlReferenceType = -1;
    LinkInfoIndex = -1;           // Used for RobotLinkInfo
    ContactStatusInfoIndex = -1;
    GoalContactPos.setZero();
    GoalContactGrad.setZero();
    }
  bool getReadyFlag(){ return ReadyFlag;}
  void setReadyFlag(const bool & _ReadyFlag ){ ReadyFlag = _ReadyFlag; }
  int  getControlReferenceType(){ return ControlReferenceType; }
  void setControlReferenceType(const int &_ControlReferenceType) { ControlReferenceType = _ControlReferenceType; }
  void setGoalContactPosNGrad(const Vector3 & _GoalContactPos, const Vector3 & _GoalContactGrad){
    GoalContactPos = _GoalContactPos;
    GoalContactGrad = _GoalContactGrad;
  }
  Vector3 getGoalContactPos(){  return GoalContactPos; }
  Vector3 getGoalContactGrad(){ return GoalContactGrad;}

  void SetInitContactStatus(const std::vector<ContactStatusInfo> &_InitContactStatus){ InitContactStatus = _InitContactStatus; }
  std::vector<ContactStatusInfo> getInitContactStatus(){return InitContactStatus;}

  void SetGoalContactStatus(const std::vector<ContactStatusInfo> & _GoalContactStatus) {GoalContactStatus = _GoalContactStatus;}
  std::vector<ContactStatusInfo> getGoalContactStatus() {return GoalContactStatus;}

  void TrajectoryUpdate(const std::vector<double> & timeTraj, const std::vector<Config> & configTraj, const std::vector<Vector3> & endeffectorTraj){
    PlannedConfigTraj = LinearPath(timeTraj, configTraj);
    std::vector<Vector> endeffectorPath;
    for (Vector3 EndEffectorPos: endeffectorTraj){
      Vector EndEffectorPosVec;
      EndEffectorPosVec.resize(3);
      EndEffectorPosVec[0] = EndEffectorPos[0];
      EndEffectorPosVec[1] = EndEffectorPos[1];
      EndEffectorPosVec[2] = EndEffectorPos[2];
      endeffectorPath.push_back(EndEffectorPosVec);
    }
    EndEffectorTraj = LinearPath(timeTraj, endeffectorPath);
  }

  bool  ReadyFlag;
  int   ControlReferenceType;
  int   LinkInfoIndex;
  int   ContactStatusInfoIndex;

  Vector3 GoalContactPos;
  Vector3 GoalContactGrad;

  LinearPath PlannedConfigTraj;
  LinearPath EndEffectorTraj;

  std::vector<ContactStatusInfo> InitContactStatus;
  std::vector<ContactStatusInfo> GoalContactStatus;
};

struct FailureStateInfo{
  FailureStateInfo(){
    FailureTime = 0.0;
    FailureStateFlag = false;
  };
  Config getFailureStateConfig() { return FailureConfig; }
  Config getFailureStateVelocity(){ return FailureVelocity; }
  bool   getFailureStateFlag() { return FailureStateFlag; }
  void FailureStateUpdate(const double & _FailureTime, const Config & _FailureConfig, const Config & _FailureVelocity){
    FailureTime = _FailureTime;
    FailureConfig = _FailureConfig;
    FailureVelocity = _FailureVelocity;
    FailureStateFlag = true;
  };
  double  FailureTime;
  Config  FailureConfig;
  Config  FailureVelocity;
  bool    FailureStateFlag;
};

struct FacetInfo{
  FacetInfo(){};
  void setFacetEdges(const std::vector<std::pair<Vector3, Vector3>> & _FacetEdges) { FacetEdges = _FacetEdges; }
  void setFacetNorm(const Vector3& _FacetNorm){ FacetNorm = _FacetNorm;}
  double ProjPoint2EdgeDist(const Vector3& _Point){
    std::vector<double> ProjPoint2Edge_vec(EdgeNorms.size());
    Vector3 Vertex2Point = _Point - FacetEdges[0].first;
    double Point2Facet = Vertex2Point.dot(FacetNorm);
    Vector3 Facet2Point = Point2Facet * FacetNorm;
    for (int i = 0; i < EdgeNorms.size(); i++){
      Vertex2Point = _Point - FacetEdges[i].first;
      Vector3 Vertex2ProjPoint = Vertex2Point - Facet2Point;
      double ProjPoint2Edge_i = Vertex2ProjPoint.dot(EdgeNorms[i]);
      ProjPoint2Edge_vec[i] = ProjPoint2Edge_i;
    }
    return *min_element(ProjPoint2Edge_vec.begin(), ProjPoint2Edge_vec.end());
  }
  std::vector<double> ProjPoint2EdgeDistVec(const Vector3& _Point){
    std::vector<double> ProjPoint2Edge_vec(EdgeNorms.size());
    Vector3 Vertex2Point = _Point - FacetEdges[0].first;
    double Point2Facet = Vertex2Point.dot(FacetNorm);
    Vector3 Facet2Point = Point2Facet * FacetNorm;
    for (int i = 0; i < EdgeNorms.size(); i++){
      Vertex2Point = _Point - FacetEdges[i].first;
      Vector3 Vertex2ProjPoint = Vertex2Point - Facet2Point;
      double ProjPoint2Edge_i = Vertex2ProjPoint.dot(EdgeNorms[i]);
      ProjPoint2Edge_vec[i] = ProjPoint2Edge_i;
    }
    return ProjPoint2Edge_vec;
  }
  void EdgesUpdate(){
    Edges.reserve(EdgeNorms.size());
    EdgesDirection.reserve(EdgeNorms.size());
    for (int i = 0; i < EdgeNorms.size(); i++){
      Vector3 Edge_i = FacetEdges[i].second - FacetEdges[i].first;
      Edges.push_back(Edge_i);
      Vector3 Edge_i_normalized;
      Edge_i.getNormalized(Edge_i_normalized);
      EdgesDirection.push_back(Edge_i_normalized);
    }
  }
  std::vector<std::pair<Vector3, Vector3>> FacetEdges;
  std::vector<Vector3> EdgeNorms;
  Vector3 FacetNorm;
  std::vector<Vector3> Edges;
  std::vector<Vector3> EdgesDirection;
};


#endif
