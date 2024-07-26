#include <mc_rbdyn/ZMP.h>
#include <mc_rtc/logging.h>
#include <stdexcept>

// Does not look nice but make sure it's not confused with system headers
#include "libqhullcpp/Qhull.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullFacetSet.h"
#include "libqhullcpp/QhullPoints.h"
#include "libqhullcpp/QhullVertexSet.h"

namespace mc_rbdyn
{

bool zmp(Eigen::Vector3d & zmpOut,
         const sva::ForceVecd & netTotalWrench,
         const Eigen::Vector3d & plane_p,
         const Eigen::Vector3d & plane_n,
         double minimalNetNormalForce) noexcept
{
  assert(minimalNetNormalForce > 0);
  const Eigen::Vector3d & force = netTotalWrench.force();
  const Eigen::Vector3d & moment_0 = netTotalWrench.couple();
  Eigen::Vector3d moment_p = moment_0 - plane_p.cross(force);
  double floorn_dot_force = plane_n.dot(force);
  // Prevent potential division by zero
  if(floorn_dot_force < minimalNetNormalForce)
  {
    mc_rtc::log::error("ZMP cannot be computed, projected force too small {}", floorn_dot_force);
    return false;
  }
  zmpOut = plane_p + plane_n.cross(moment_p) / floorn_dot_force;
  return true;
}

Eigen::Vector3d zmp(const sva::ForceVecd & netTotalWrench,
                    const Eigen::Vector3d & plane_p,
                    const Eigen::Vector3d & plane_n,
                    double minimalNetNormalForce)
{
  if(minimalNetNormalForce <= 0)
  {
    mc_rtc::log::error_and_throw("ZMP cannot be computed: the minimalNetNormalForce must be >0 (divide by zero)");
  }

  Eigen::Vector3d zmpOut;
  if(!zmp(zmpOut, netTotalWrench, plane_p, plane_n, minimalNetNormalForce))
  {
    mc_rtc::log::error_and_throw("ZMP cannot be computed");
  }
  return zmpOut;
}

Eigen::Vector3d zmp(const sva::ForceVecd & netWrench, const sva::PTransformd & zmpFrame, double minimalNetNormalForce)
{
  Eigen::Vector3d n = zmpFrame.rotation().row(2);
  Eigen::Vector3d p = zmpFrame.translation();
  return zmp(netWrench, p, n, minimalNetNormalForce);
}

bool zmp(Eigen::Vector3d & zmpOut,
         const sva::ForceVecd & netWrench,
         const sva::PTransformd & zmpFrame,
         double minimalNetNormalForce) noexcept
{
  Eigen::Vector3d n = zmpFrame.rotation().row(2);
  Eigen::Vector3d p = zmpFrame.translation();
  return zmp(zmpOut, netWrench, p, n, minimalNetNormalForce);
}

Eigen::Matrix3d compute_plane_frame(const std::vector<Eigen::Vector3d> & hull,const Eigen::Vector3d & plane_p,const Eigen::Vector3d & plane_n)
{
  Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
  R.col(2) << plane_n / plane_n.norm();

  size_t i = 0;
  while (i < hull.size())
  {
      auto & p1 = hull[i];
      
      if ((p1 - plane_p).norm() > 1e-4)
      {
        R.col(0) << (p1 - plane_p).cross(plane_n);
        R.col(0) /= R.col(0).norm();
        R.col(1) << (p1 - plane_p) / (p1 - plane_p).norm() ;
      }
      i+=1;
  }
  return R;
}

std::vector<Eigen::Vector3d> plane_hull(const std::vector<Eigen::Vector3d> & pts, const Eigen::Vector3d & plane_p, const Eigen::Vector3d & plane_n)
{
  Eigen::Matrix3d R = compute_plane_frame(pts,plane_p,plane_n);
  std::vector<double> pts_2d;
  pts_2d.reserve(pts.size() * 2);
  for( auto & pt_3d : pts)
  {
   const Eigen::Vector2d pt_2d = (R.transpose() * (pt_3d - plane_p)).segment(0,2);
   pts_2d.push_back(pt_2d.x());
   pts_2d.push_back(pt_2d.y());
  }

  // Run qhull
  orgQhull::Qhull qhull;
  qhull.runQhull("", 2, static_cast<int>(pts.size()), pts_2d.data(), "Qt");
  std::vector<Eigen::Vector3d> hull_pts_3d = {};
  std::vector<int> hull_indx = {};
  
  auto facet = qhull.facetList().first();
  auto p = facet.vertices().first();
  hull_pts_3d.push_back( plane_p + R * (Eigen::Vector3d{p.point().coordinates()[0],p.point().coordinates()[1],0}) );
  hull_indx.push_back(p.id());
  p = facet.vertices().second();
  hull_pts_3d.push_back( plane_p + R * (Eigen::Vector3d{p.point().coordinates()[0],p.point().coordinates()[1],0}) );
  hull_indx.push_back(p.id());
  std::vector<int> facets_indx = {facet.id()};
  facet = facet.neighborFacets().first();
  facets_indx.push_back(facet.id());

  while (hull_indx.size() < qhull.vertexCount())
  {
    auto p = facet.vertices().first();
    if (std::find(hull_indx.begin(), hull_indx.end(), p.id()) != hull_indx.end()) {
        p = facet.vertices().second();
    }
    const Eigen::Vector3d pt3d = plane_p + R * (Eigen::Vector3d{p.point().coordinates()[0],p.point().coordinates()[1],0});
    if( (pt3d - hull_pts_3d.back()).norm() > 1e-4 )
    {
      hull_pts_3d.push_back( pt3d);
    }
  
    hull_indx.push_back(p.id());
    
    if (std::find(facets_indx.begin(), facets_indx.end(), facet.neighborFacets().first().id()) != facets_indx.end()) {
        facet = facet.neighborFacets().back();
    }
    else
    {
        facet = facet.neighborFacets().first();
    }
    facets_indx.push_back(facet.id());
  }
  
  return hull_pts_3d;
}

bool pick_2d_extreme_rays(std::pair<Eigen::Vector2d,size_t> & u_low, std::pair<Eigen::Vector2d,size_t> & u_high, std::vector<Eigen::Vector2d> rays)
{
  assert(rays.size() > 1);
  if (rays.size() == 2)
  {
      u_low = std::pair<Eigen::Vector2d,size_t> {rays[0],0};
      u_high = std::pair<Eigen::Vector2d,size_t> {rays[1],1};
      return true;
  }

  const size_t none_indx = rays.size();
  
  u_high = {Eigen::Vector2d::Zero(),none_indx};
  u_low = {rays.back(),rays.size() - 1}; rays.erase(rays.end());
  while ( u_high.second == none_indx)
  {
      const size_t idx_ray =  rays.size() - 1;
      const Eigen::Vector2d ray = rays.back(); rays.erase(rays.end());
      const double c = (Eigen::Vector3d{u_low.first.x(),u_low.first.y(),0}.cross(Eigen::Vector3d{ray.x(),ray.y(),0})).z();
      if (std::abs(c) < 1e-4 && rays.size()==0)
      {
          
        return false;
          
      }
      else if (c < 0)
      {
          u_low.second = idx_ray;
          u_high.second = u_low.second;
          u_low.first = ray;
          u_high.first = u_low.first;
      }
      else
      {
          u_high.second = idx_ray;
          u_high.first = ray;
      }
  }
  for(size_t indx_u = 0 ; indx_u < rays.size() ; indx_u ++)
  {
      const double c1 = (Eigen::Vector3d{u_low.first.x(),u_low.first.y(),0}.cross(Eigen::Vector3d{rays[indx_u].x(),rays[indx_u].y(),0})).z();
      const double c2 = (Eigen::Vector3d{u_high.first.x(),u_high.first.y(),0}.cross(Eigen::Vector3d{rays[indx_u].x(),rays[indx_u].y(),0})).z();
      if (c1 < 0 && c2 < 0)
      {
          u_low.second = indx_u;
          u_low.first = rays[indx_u];
      }
      else if ( c1 > 0 && c2 > 0)
      {
          u_high.second = indx_u;
          u_high.first = rays[indx_u];
      }
      else if (c1 < 0 && c2 > 0)
      {
          return false;
      }
  }

  return true;
}

std::vector<Eigen::Vector3d> convert_cone2d_to_vertices(const std::vector<Eigen::Vector3d> & vertices,
                                                        const std::vector<Eigen::Vector3d> & rays,
                                                        const Eigen::Vector3d & plane_p,
                                                        const Eigen::Vector3d & plane_n)
{
  const double big_dist = 10;
  std::vector<Eigen::Vector3d> conv_vertices = {};

  if (rays.size() == 0)
  {
      return vertices;
  }
  const Eigen::Matrix3d R = compute_plane_frame(vertices,plane_p,plane_n);
  std::pair<Eigen::Vector2d,size_t> low;
  std::pair<Eigen::Vector2d,size_t> high;
  std::vector<Eigen::Vector2d> pt_2d;
  for(auto & p : rays)
  {
    pt_2d.push_back((R.transpose() * (p - plane_p)).segment(0,2));
  }
  if (!pick_2d_extreme_rays(low,high,pt_2d))
  {
    conv_vertices.push_back(plane_p + R * Eigen::Vector3d{-big_dist,-big_dist,0});  
    conv_vertices.push_back(plane_p + R * Eigen::Vector3d{-big_dist,+big_dist,0});
    conv_vertices.push_back(plane_p + R * Eigen::Vector3d{+big_dist,+big_dist,0});  
    conv_vertices.push_back(plane_p + R * Eigen::Vector3d{+big_dist,-big_dist,0}); 
    return conv_vertices;
  }

  const Eigen::Vector3d r0 = rays[low.second].normalized();
  const Eigen::Vector3d r1 = rays[high.second].normalized();

  conv_vertices = vertices;
  for (auto & v : vertices)
  {    
    conv_vertices.push_back(v + r0 * big_dist);
    conv_vertices.push_back(v + r1 * big_dist);

  }

  return plane_hull(conv_vertices,plane_p,plane_n);
}

void compute_polygon_shadow(std::vector<Eigen::Vector3d> & vertices,
                            std::vector<Eigen::Vector3d> & rays,
                            const std::vector<Eigen::Vector3d> & light_poly,
                            const std::vector<Eigen::Vector3d> & ground_poly,
                            const Eigen::Vector3d & plane_p,
                            const Eigen::Vector3d & plane_n)
{
  
    const auto R = compute_plane_frame(light_poly,plane_p,plane_n);

    std::vector<Eigen::Vector2d> mink_diff;
    for (auto & gv : ground_poly )
    {
      for(auto & lv : light_poly)
      {
        mink_diff.push_back( (R.transpose()*( gv - lv)).segment(0,2) );
      }
    }
    
    std::pair<Eigen::Vector2d,size_t> low;
    std::pair<Eigen::Vector2d,size_t> high;

    vertices.clear();
    rays.clear();

    if (!pick_2d_extreme_rays(low,high,mink_diff))
    {
      const double big_dist = 10;
      vertices.push_back(plane_p + R * Eigen::Vector3d{-big_dist,-big_dist,0});  
      vertices.push_back(plane_p + R * Eigen::Vector3d{-big_dist,+big_dist,0});
      vertices.push_back(plane_p + R * Eigen::Vector3d{+big_dist,+big_dist,0});  
      vertices.push_back(plane_p + R * Eigen::Vector3d{+big_dist,-big_dist,0});  
      return ;
    }

    vertices = ground_poly;
    rays = {   R * Eigen::Vector3d{low.first[0],low.first[1],0},  R * Eigen::Vector3d{high.first[0],high.first[1],0}}; 

}

void project_contact_cone(const mc_rbdyn::Robot & r,
                          std::vector<Eigen::Vector3d> & pos_pts,
                          std::vector<Eigen::Vector3d> & neg_pts,
                          const mc_rbdyn::Contact & contact,
                          const Eigen::Vector3d & plane_p,
                          const Eigen::Vector3d & plane_n)
{
  std::vector<sva::PTransformd> contact_pts = contact.r1Points();
  sva::PTransformd X_0_b;

  if(r.hasSurface(contact.r1Surface()->name()))
  {
    X_0_b = contact.X_b_s().inv() * r.surfacePose(contact.r1Surface()->name()) ;
  }
  else
  {
    contact_pts = contact.r2Points();
    X_0_b = contact.X_b_s().inv() * r.surfacePose(contact.r2Surface()->name()) ;
  }

  for (auto & pt : contact_pts)
  {
 
    const Eigen::Vector3d pt_0 = (pt * X_0_b).translation();

    for (auto & fc : contact.force_span())
    {
      auto f = (contact.X_b_s() * X_0_b).rotation() * fc;
      const double n_dot_f = plane_n.transpose() * f;
      if( n_dot_f != 0)
      {
        const Eigen::Vector3d p = plane_p + plane_n.cross( (pt_0 - plane_p).cross(f) )/n_dot_f;
        if(n_dot_f >=0)
        {
          pos_pts.push_back(p);
        }
        else
        {
          neg_pts.push_back(p);
        }     
      }
    }
  }
}

std::vector<std::vector<Eigen::Vector3d>> full_support_region(const mc_rbdyn::Robot & robot,
                                                              const std::vector<mc_rbdyn::Contact> & contacts,
                                                              const Eigen::Vector3d & plane_p,
                                                              const Eigen::Vector3d & plane_n)
{
  const auto n = plane_n.normalized();
  std::vector<Eigen::Vector3d> pos_pts = {};
  std::vector<Eigen::Vector3d> neg_pts = {};
  for (auto & c : contacts)
  {
    project_contact_cone(robot,pos_pts,neg_pts,c,plane_p,n);
  }


  std::vector<Eigen::Vector3d> neg_poly = {};
  std::vector<Eigen::Vector3d> pos_poly = {};

  if(pos_pts.size() > 2){pos_poly = plane_hull(pos_pts,plane_p,n);}
  if(neg_pts.size() > 2){neg_poly = plane_hull(neg_pts,plane_p,n);}

  if(pos_poly.size() != 0 && neg_poly.size() != 0)
  {
    std::vector<Eigen::Vector3d> pos_r= {};
    std::vector<Eigen::Vector3d> pos_v= {};
    compute_polygon_shadow(pos_v,pos_r,neg_poly,pos_poly,plane_p,n);

    std::vector<Eigen::Vector3d> neg_r= {};
    std::vector<Eigen::Vector3d> neg_v= {};
    compute_polygon_shadow(neg_v,neg_r,pos_poly,neg_poly,plane_p,n);
    
    std::vector<std::vector<Eigen::Vector3d>> regions = {};

    if (pos_r.size() == 0 )
    {
      regions.push_back(pos_v);
    }
    else
    {
      regions.push_back(convert_cone2d_to_vertices(pos_v,pos_r,plane_p,n));
    }

    if (neg_r.size() == 0 )
    {
      regions.push_back(neg_v);
    }
    else
    {
      regions.push_back(convert_cone2d_to_vertices(neg_v,neg_r,plane_p,n));
    }

    return regions;

  }
  else
  {
    return {pos_poly,neg_poly};
  }

}

} // namespace mc_rbdyn
