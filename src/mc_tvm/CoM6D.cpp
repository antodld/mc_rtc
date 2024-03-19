/*
 * Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include "../../include/mc_tvm/CoM6D.h"

#include <mc_tvm/Robot.h>

namespace mc_tvm
{

CoM6D::CoM6D(NewCoMToken, Robot & robot) : robot_(robot)
{
  // clang-format off
  registerUpdates(
                  Update::CoM6D, &CoM6D::updateCoM,
                  Update::Jacobian, &CoM6D::updateJacobian,
                  Update::Velocity, &CoM6D::updateVelocity,
                  Update::NormalAcceleration, &CoM6D::updateNormalAcceleration,
                  Update::Acceleration, &CoM6D::updateAcceleration,
                  Update::JDot, &CoM6D::updateJDot);
  // clang-format off

  addOutputDependency(Output::CoM6D, Update::CoM6D);
  addInputDependency(Update::CoM6D, robot_, Robot::Output::FK);

  addOutputDependency(Output::Jacobian, Update::Jacobian);
  addInputDependency(Update::Jacobian, robot_, Robot::Output::FV);

  addOutputDependency(Output::Velocity, Update::Velocity);
  addInputDependency(Update::Velocity, robot_, Robot::Output::FV);

  addOutputDependency(Output::NormalAcceleration, Update::NormalAcceleration);
  addInputDependency(Update::NormalAcceleration, robot_, Robot::Output::NormalAcceleration);

  addOutputDependency(Output::Acceleration, Update::Acceleration);
  addInputDependency(Update::Acceleration, robot_, Robot::Output::FA);

  addOutputDependency(Output::JDot, Update::JDot);
  addInputDependency(Update::JDot, robot_, Robot::Output::FV);

  addInternalDependency(Update::Velocity, Update::Jacobian);
  addInternalDependency(Update::NormalAcceleration, Update::Jacobian);
}

void CoM6D::updateCoM()
{
  const auto & r = robot().robot();
  if(r.mass() > 0)
  {
    com_ = r.mbc().com;
  }
  else
  {
    com_ = r.posW();
  }
}

void CoM6D::updateVelocity()
{
  const auto & r = robot().robot();
  velocity_ = r.mbc().comVel.vector();
}

void CoM6D::updateNormalAcceleration()
{
  const auto & r = robot().robot();
  normalAcceleration_ = r.mbc().Jcomdot * rbd::paramToVector(r.mb(), r.mbc().alpha)
                        + sva::MotionVecd(Eigen::Vector3d::Zero(), r.mbc().comVel.angular().cross(r.mbc().comVel.linear())).vector();
}

void CoM6D::updateAcceleration()
{
  const auto & r = robot().robot();
  if(r.mass() > 0)
  {
    acceleration_ = r.mbc().comAcc.vector();
  }
  else
  {
    acceleration_ = Eigen::Vector6d::Zero();
  }
}

void CoM6D::updateJacobian()
{
  const auto & r = robot().robot();
  jac_ = r.mbc().Jcom;
}

void CoM6D::updateJDot()
{
  const auto & r = robot().robot();
  jacDot_ = r.mbc().Jcomdot;
}

} // namespace mc_tvm
