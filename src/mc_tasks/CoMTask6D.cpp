/*
 * Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include <mc_tasks/CoMTask6D.h>

#include <mc_tasks/MetaTaskLoader.h>

#include <mc_rtc/ConfigurationHelpers.h>
#include <mc_rtc/gui/Transform.h>

#include <mc_tvm/CoM6DFunction.h>


namespace mc_tasks
{

static inline mc_rtc::void_ptr_caster<tasks::qp::CoM6DTask> tasks_error{};
static inline mc_rtc::void_ptr_caster<mc_tvm::CoM6DFunction> tvm_error{};

CoMTask6D::CoMTask6D(const mc_rbdyn::Robots & robots, unsigned int robotIndex, double stiffness, double weight)
: TrajectoryTaskGeneric(robots, robotIndex, stiffness, weight), robot_index_(robotIndex)
{
  const mc_rbdyn::Robot & robot = robots.robot(robotIndex);
  switch(backend_)
  {

    case Backend::Tasks:
      finalize<Backend::Tasks, tasks::qp::CoM6DTask>(robots.mbs(), static_cast<int>(robotIndex), robot.mbc().com);
        break;
    case Backend::TVM:
      finalize<Backend::TVM, mc_tvm::CoM6DFunction>(robot);
      break;
    default:
      mc_rtc::log::error_and_throw("[CoMTask6D] Not implemented for solver backend: {}", backend_);
  }
  type_ = "com6d";
  name_ = "com6d_" + robots.robot(robot_index_).name();
}

void CoMTask6D::reset()
{
  TrajectoryTaskGeneric::reset();
  const mc_rbdyn::Robot & robot = robots.robot(rIndex);
  com(robot.mbc().com);
}

void CoMTask6D::load(mc_solver::QPSolver & solver, const mc_rtc::Configuration & config)
{
  TrajectoryBase::load(solver, config);
  const auto & robot_name = solver.robots().robot(robot_index_).name();
  this->com(solver.robots().robot(robot_name).mbc().com);
  if(config.has("com6d")) { this->com(config("com6d")); }
  if(config.has("above"))
  {
    auto surfaces = mc_rtc::fromVectorOrElement<std::string>(config("above"), "surfaces");
    double height = config("above")("height",com().translation().z());

    auto com = this->com();
    Eigen::Vector3d target = Eigen::Vector3d::Zero();
    auto & robot = robotFromConfig(config, solver.robots(), name());
    for(const auto & s : surfaces) { target += robot.surface(s).X_0_s(robot).translation(); }
    target /= static_cast<double>(surfaces.size());
    this->com(sva::PTransformd(com.rotation(),{target.x(), target.y(), target.z() + height}));
  }
  if(config.has("move_com")) { this->move_com(config("move_com")); }
  if(config.has("move_com_world"))
  { 
    const sva::PTransformd X_offset = config("move_com_world",sva::PTransformd::Identity());
    this->move_com(sva::PTransformd(com().rotation()) * X_offset * sva::PTransformd(com().rotation()).inv()); 
  }

}

void CoMTask6D::move_com(const sva::PTransformd & move)
{
  com(move * com());
}

void CoMTask6D::com(const sva::PTransformd & com)
{
  switch(backend_)
  {
    case Backend::Tasks:
      tasks_error(errorT)->com(com);
      break;
    case Backend::TVM:
      tvm_error(errorT)->com(com);
      break;
    default:
      break;
  }
}

const sva::PTransformd & CoMTask6D::com() const
{
  switch(backend_)
  {
    case Backend::Tasks:
      return tasks_error(errorT)->com();
    case Backend::TVM:
      return tvm_error(errorT)->com();
    default:
      mc_rtc::log::error_and_throw("Not implemented");
  }
}

const sva::PTransformd & CoMTask6D::actual() const
{
  switch(backend_)
  {
    case Backend::Tasks:
      return tasks_error(errorT)->actual();
    case Backend::TVM:
      return tvm_error(errorT)->actual();
    default:
      mc_rtc::log::error_and_throw("Not implemented");
  }
}

void CoMTask6D::addToLogger(mc_rtc::Logger & logger)
{
  TrajectoryBase::addToLogger(logger);
  MC_RTC_LOG_HELPER(name_ + "_pos", actual);
  MC_RTC_LOG_GETTER(name_ + "_target", com);
}

void CoMTask6D::addToGUI(mc_rtc::gui::StateBuilder & gui)
{
  TrajectoryBase::addToGUI(gui);
  gui.addElement({"Tasks", name_},
                 mc_rtc::gui::Transform(
                     "com_target", [this]() -> const sva::PTransformd & { return this->com(); },
                     [this](const sva::PTransformd & com) { this->com(com); }),
                 mc_rtc::gui::Transform("com", [this]() -> const sva::PTransformd & { return this->actual(); }));
}

} // namespace mc_tasks

namespace
{

static auto registered = mc_tasks::MetaTaskLoader::register_load_function(
    "com6d",
    [](mc_solver::QPSolver & solver, const mc_rtc::Configuration & config)
    {
      std::string robot_name = solver.robot().name();
      if(config.has("robot")){ config("robot",robot_name);}

      const auto robotIndex = solver.robots().robotIndex(robot_name);
      auto t = std::make_shared<mc_tasks::CoMTask6D>(solver.robots(), robotIndex);
      t->load(solver, config);
      return t;
    });
}
