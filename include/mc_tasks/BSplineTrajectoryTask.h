#pragma once
#include <mc_tasks/SplineTrajectoryTask.h>
#include <mc_trajectory/BSpline.h>

namespace mc_tasks
{

/*! \brief Track a bezier curve with a robot surface.
 *
 * SplineTrajectoryTask takes care of much of the logic (task target updates,
 * orientation waypoint handling, etc.), and brings in all the functionalities
 * from TrajectoryTaskGeneric. This BSplineTrajectoryTask only handles specific
 * aspect of the Bezier curve.
 */
struct MC_TASKS_DLLAPI BSplineTrajectoryTask : public SplineTrajectoryTask<BSplineTrajectoryTask>
{
public:
  /**
   * \brief Creates a trajectory that follows a bspline curve
   *
   * \param robots Robots controlled by the task
   * \param robotIndex Which robot is controlled
   * \param surfaceName Surface controlled by the task, should belong to the
   * \ontrolled robot
   * \param duration Duration of motion (eg time it takes to go from the current
   * \urface position to the curve's final point)
   * \param stiffness Stiffness of the underlying TrajectoryTask (position and
   * \rientation)
   * \param posW Task weight (position)
   * \param oriW Task weight (orientation)
   * \param target Final world pose to reach
   * \param posWp Waypoints in position
   * \param oriWp Waypoints in orientation specified as pairs of (time,
   * orientation). The surface orientation will be interpolated in-between
   * waypoints.
   */
  BSplineTrajectoryTask(const mc_rbdyn::Robots & robots,
                        unsigned int robotIndex,
                        const std::string & surfaceName,
                        double duration,
                        double stiffness,
                        double posW,
                        double oriW,
                        const sva::PTransformd & target,
                        const std::vector<Eigen::Vector3d> & posWp,
                        const std::vector<std::pair<double, Eigen::Matrix3d>> & oriWp = {});

  /*! \brief const accessor to the underlying spline (used by SplineTrajectoryTask)
   *
   * \returns The spline
   */
  const mc_trajectory::BSpline & spline() const
  {
    return bspline;
  };
  /*! \brief accessor to the underlying spline (used by SplineTrajectoryTask)
   *
   * \returns The spline
   */
  mc_trajectory::BSpline & spline()
  {
    return bspline;
  };

  /*! \brief Sets the curve target pose
   * \param target Target pose for the curve
   */
  void target(const sva::PTransformd & target);
  /*!
   * \brief Gets the target position of the curve
   *
   * \returns target position
   */
  Eigen::Vector3d target() const;

  /*! \brief Add interactive GUI elements to control the curve waypoints
   */
  void addToGUI(mc_rtc::gui::StateBuilder & gui);

  /** \brief Control points for the bezier curve (position)
   *
   * \param posWp Vector of position control points for the bezier curve. Should
   * not include the starting and target position.
   */
  void posWaypoints(const std::vector<Eigen::Vector3d> & posWp);

protected:
  mc_trajectory::BSpline bspline;
};

} // namespace mc_tasks
