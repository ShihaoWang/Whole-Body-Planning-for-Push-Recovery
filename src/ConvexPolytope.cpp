#include "RobotInfo.h"
#include "CommonHeader.h"
#include <set>

#include <cilantro/spatial/convex_polytope.hpp>
#include <cilantro/utilities/point_cloud.hpp>
#include <cilantro/visualization/visualizer.hpp>
#include <cilantro/visualization/common_renderables.hpp>
#include <cilantro/spatial/flat_convex_hull_3d.hpp>
static double eps_dist = 0.001;        // 1 mm

FacetInfo FlatContactHullGeneration(const std::vector<Vector3> & ContactPoints, int & FacetFlag){
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
    FacetFlag = 0;
    return FacetInfoObj;
  } else FacetFlag = 1;
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
