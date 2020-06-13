#include "RobotInfo.h"
#include "CommonHeader.h"
#include <set>

#include <cilantro/spatial/convex_polytope.hpp>
#include <cilantro/utilities/point_cloud.hpp>
#include <cilantro/visualization/visualizer.hpp>
#include <cilantro/visualization/common_renderables.hpp>
#include <cilantro/spatial/flat_convex_hull_3d.hpp>
static double eps_dist = 0.001;        // 1 mm

FacetInfo FlatConvexHullGeneration(const std::vector<Vector3> & ContactPoints){
  FacetInfo FacetInfoObj;
  std::vector<Eigen::Vector3f> Vertices;
  Vertices.reserve(ContactPoints.size());
  for (int i = 0; i < ContactPoints.size(); i++)
      Vertices.emplace_back(ContactPoints[i].x,ContactPoints[i].y,ContactPoints[i].z);
  cilantro::PointCloud3f datacloud(Vertices);
  cilantro::FlatConvexHull3f flat_hull(datacloud.points, true, false);
  datacloud.points = flat_hull.reconstruct<2>(flat_hull.project<2>(datacloud.points));
  const auto& face_v_ind = flat_hull.getFacetVertexIndices();
  const int SegmentNumber = face_v_ind.size();
  cilantro::VectorSet3f src_points(3, SegmentNumber);
  cilantro::VectorSet3f dst_points(3, SegmentNumber);
  if(!face_v_ind.size()){
    return FacetInfoObj;
  } else FacetInfoObj.setFacetValidFlag(true);
  // src_points/dst_points will always ouputs convex hull edges but these edges may not be in the counterclockwise order.
  Eigen::Vector3f FirstA = flat_hull.getVertices3D().col(face_v_ind[0][0]);
  Eigen::Vector3f FirstB = flat_hull.getVertices3D().col(face_v_ind[0][1]);
  int ReverseOrder = 0;
  Vector3 NormalVec(0.0, 0.0, 1.0);     // Initialize it to be a unit vector
  int NormalVecFlag = 0;
  for (size_t i = 1; i < SegmentNumber; i++)
  {
    Eigen::Vector3f SecondA = flat_hull.getVertices3D().col(face_v_ind[i][0]);
    Eigen::Vector3f SAmFB = SecondA - FirstB;
    if(SAmFB.norm()<eps_dist){
      Eigen::Vector3f SecondB = flat_hull.getVertices3D().col(face_v_ind[i][1]);
      Eigen::Vector3f FirstEdge = FirstB - FirstA;
      Eigen::Vector3f SecondEdge = SecondB - SecondA;

      Eigen::Vector3f EdgeNorm = FirstEdge.cross(SecondEdge);
      float NormProj = EdgeNorm(2);
      if(NormProj<0){
        ReverseOrder = 1;
        Vector3 FaceNorm(-1.0 * EdgeNorm(0), -1.0 * EdgeNorm(1), -1.0 * EdgeNorm(2));
        NormalVec.setNormalized(FaceNorm);
        break;
      }
      else {
        switch (NormalVecFlag)
        {
          case 0:
          Vector3 FaceNorm(EdgeNorm(0), EdgeNorm(1), EdgeNorm(2));
          NormalVec.setNormalized(FaceNorm);
          NormalVecFlag = 1;
          break;
        }
      }
    }
  }
  FacetInfoObj.setFacetNorm(NormalVec);
  FacetInfoObj.FacetEdges.reserve(SegmentNumber);
  FacetInfoObj.EdgeNorms.reserve(SegmentNumber);
  for (int i = 0; i < SegmentNumber; i++){
    Eigen::Vector3f FirstA = flat_hull.getVertices3D().col(face_v_ind[i][0]);
    Eigen::Vector3f FirstB = flat_hull.getVertices3D().col(face_v_ind[i][1]);
    if(!ReverseOrder){
      Vector3 EdgeFront(FirstA(0), FirstA(1), FirstA(2));
      Vector3 EdgeBack(FirstB(0), FirstB(1), FirstB(2));
      FacetInfoObj.FacetEdges.push_back(make_pair(EdgeFront, EdgeBack));

      Vector3 EdgeFront2Back = EdgeBack - EdgeFront;
      Vector3 EdgeNorm_i = cross(FacetInfoObj.FacetNorm, EdgeFront2Back);
      FacetInfoObj.EdgeNorms.push_back(EdgeNorm_i/EdgeNorm_i.norm());
    } else {
      // src_point/dst_point order have to be switched.
      // There is no need to change src_point/dst_point order
      Vector3 EdgeFront(FirstB(0), FirstB(1), FirstB(2));
      Vector3 EdgeBack(FirstA(0), FirstA(1), FirstA(2));
      FacetInfoObj.FacetEdges.push_back(make_pair(EdgeFront, EdgeBack));

      Vector3 EdgeFront2Back = EdgeBack - EdgeFront;
      Vector3 EdgeNorm_i = cross(FacetInfoObj.FacetNorm, EdgeFront2Back);
      FacetInfoObj.EdgeNorms.push_back(EdgeNorm_i/EdgeNorm_i.norm());
    }
  }
  return FacetInfoObj;
}

static double ZCoordinate(const double & _X, const double & _Y, const std::vector<Vector3> & ContactPoints){
  // This function is used to select out the ZCoodrinate given _X and _Y.
  int CPNumber = ContactPoints.size();
  std::vector<double> CPDiff(CPNumber);
  for (int i = 0; i < CPNumber; i++)
    CPDiff[i] = (ContactPoints[i].x - _X) * (ContactPoints[i].x - _X) +
                (ContactPoints[i].y - _Y) * (ContactPoints[i].y - _Y);
  int CPDiffIndex = std::distance(CPDiff.begin(), std::min_element(CPDiff.begin(), CPDiff.end()));
  return ContactPoints[CPDiffIndex].z;
}

static PIPInfo PIPGeneratorInner(const Vector3 & EdgeA, const Vector3 & EdgeB, const Vector3 & COM, const Vector3 & COMVel){
  // Here the inputs are EdgeA and EdgeB, which should be chosen to be on the positive dirction of certain edge.
  // Now it is the CoM projection and basically the job is to find the orthogonal intersection point
  Vector3 EdgeA2COM = COM - EdgeA;
  Vector3 EdgeA2B = EdgeB - EdgeA;                      // x
  Vector3 x_unit, y_unit, z_unit;
  EdgeA2B.getNormalized(x_unit);
  double ABLength = EdgeA2B.norm();
  double t = EdgeA2COM.dot(x_unit);
  double InterOnSegRatio = t/ABLength;
  Vector3 COMOnEdge = EdgeA + t * x_unit;
  Vector3 COMOnEdge2COM = COM - COMOnEdge;              // y
  COMOnEdge2COM.getNormalized(y_unit);
  z_unit = -1.0*cross(x_unit, y_unit);                  // z

  double L = COMOnEdge2COM.norm();
  double Ldot = COMVel.dot(y_unit);
  double thetadot = COMVel.dot(z_unit)/L;

  Vector3 g(0.0, 0.0, -9.81);
  Vector3 EdgeA2BProj(EdgeB.x - EdgeA.x, EdgeB.y - EdgeA.y, 0.0);
  Vector3 z_prime_unit = cross(EdgeA2BProj, -1.0 * g);
  z_prime_unit = z_prime_unit/z_prime_unit.norm();
  Vector3 y_prime_unit = cross(z_prime_unit, x_unit);
  double theta = atan2(-1.0 * COMOnEdge2COM.dot(z_prime_unit), COMOnEdge2COM.dot(y_prime_unit));
  double g_proj = -g.dot(y_prime_unit);

  double VertY = EdgeB.z - EdgeA.z;
  double VertX = (EdgeB.x - EdgeA.x)*(EdgeB.x - EdgeA.x) + (EdgeB.y - EdgeA.y) * (EdgeB.y - EdgeA.y);
  VertX = sqrt(VertX);
  double g_angle = abs(atan2(VertY, VertX))/3.1415926535 * 180.0;       // Expressed in degrees

  PIPInfo PIPObj(L, Ldot, theta, thetadot, g_proj, g_angle);
  PIPObj.setPrimeUnits(x_unit, y_prime_unit, z_prime_unit);
  PIPObj.setUnits(x_unit, y_unit, z_unit);    // Here x-y-z do not obey right hand rule.
  PIPObj.setEdgeAnB(EdgeA, EdgeB);
  PIPObj.setIntersection(COMOnEdge);
  if((InterOnSegRatio>=0.0)&&(InterOnSegRatio<=1.0))
          PIPObj.onFlag = true;
  else    PIPObj.onFlag = false;
  return PIPObj;
}

std::vector<PIPInfo> PIPGenerator(const std::vector<Vector3> & ContactPoints, const Vector3 & COMPos, const Vector3 & COMVel){
  // This method utilizes the Projected Support Polygon.
  std::vector<PIPInfo> PIPTotal;
  std::vector<Vector3> SPVertices;
  SPVertices.reserve(ContactPoints.size());
  for (int i = 0; i < ContactPoints.size(); i++)
    SPVertices.push_back(Vector3(ContactPoints[i].x, ContactPoints[i].y, 0.0));
  FacetInfo FlatObj= FlatConvexHullGeneration(SPVertices);
  if(!FlatObj.getFacetValidFlag()) return PIPTotal;
  else {
    for (int i = 0; i < FlatObj.FacetEdges.size(); i++){
      double Edge_i_First_z =   ZCoordinate(FlatObj.FacetEdges[i].first.x,  FlatObj.FacetEdges[i].first.y,    ContactPoints);
      double Edge_i_Second_z =  ZCoordinate(FlatObj.FacetEdges[i].second.x, FlatObj.FacetEdges[i].second.y,   ContactPoints);
      Vector3 FirstEdge(  FlatObj.FacetEdges[i].first.x,   FlatObj.FacetEdges[i].first.y,  Edge_i_First_z);
      Vector3 SecondEdge( FlatObj.FacetEdges[i].second.x,  FlatObj.FacetEdges[i].second.y, Edge_i_Second_z);
      PIPInfo PIPObj = PIPGeneratorInner(FirstEdge, SecondEdge, COMPos, COMVel);
      if(PIPObj.g_angle<70) PIPTotal.push_back(PIPObj);
    }
  }
  return PIPTotal;
}
