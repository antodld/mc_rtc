/*
 * Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include <mc_tvm/CoM6DFunction.h>

#include <mc_rbdyn/Robot.h>
#include <mc_tvm/Robot.h>

namespace mc_tvm
{

CoM6DFunction::CoM6DFunction(const mc_rbdyn::Robot & robot)
: tvm::function::abstract::Function(6), comAlgo_(robot.tvmRobot().com6DAlgo())
{
  reset();
  // clang-format off
  registerUpdates(Update::Value, &CoM6DFunction::updateValue,
                  Update::Velocity, &CoM6DFunction::updateVelocity,
                  Update::Jacobian, &CoM6DFunction::updateJacobian,
                  Update::NormalAcceleration, &CoM6DFunction::updateNormalAcceleration,
                  Update::JDot, &CoM6DFunction::updateJDot);
  // clang-format on
  addOutputDependency<CoM6DFunction>(Output::Value, Update::Value);
  addOutputDependency<CoM6DFunction>(Output::Velocity, Update::Velocity);
  addOutputDependency<CoM6DFunction>(Output::Jacobian, Update::Jacobian);
  addOutputDependency<CoM6DFunction>(Output::NormalAcceleration, Update::NormalAcceleration);
  addOutputDependency<CoM6DFunction>(Output::JDot, Update::JDot);
  auto & tvm_robot = robot.tvmRobot();
  addVariable(tvm_robot.q(), false);
  auto & comAlgo = tvm_robot.com6DAlgo();
  addInputDependency<CoM6DFunction>(Update::Value, comAlgo, mc_tvm::CoM6D::Output::CoM6D);
  addInputDependency<CoM6DFunction>(Update::Velocity, comAlgo, mc_tvm::CoM6D::Output::Velocity);
  addInputDependency<CoM6DFunction>(Update::Jacobian, comAlgo, mc_tvm::CoM6D::Output::Jacobian);
  addInputDependency<CoM6DFunction>(Update::NormalAcceleration, comAlgo, mc_tvm::CoM6D::Output::NormalAcceleration);
  addInputDependency<CoM6DFunction>(Update::JDot, comAlgo, mc_tvm::CoM6D::Output::JDot);
}

void CoM6DFunction::reset()
{
  com_ = comAlgo_.robot().robot().mbc().com;
  refVel_.setZero();
  refAccel_.setZero();
}

void CoM6DFunction::updateValue()
{
  const auto X_target_current = comAlgo_.com() * com_.inv();
  value_ = sva::transformVelocity(X_target_current).vector();
}

void CoM6DFunction::updateVelocity()
{
  velocity_ = (comAlgo_.velocity() - refVel_);
}

void CoM6DFunction::updateJacobian()
{
  splitJacobian(comAlgo_.jacobian(), comAlgo_.robot().q());
}

void CoM6DFunction::updateNormalAcceleration()
{
  normalAcceleration_ = (comAlgo_.normalAcceleration() - refAccel_);
}

void CoM6DFunction::updateJDot()
{
  splitJacobian(comAlgo_.JDot(), comAlgo_.robot().q());
}

} // namespace mc_tvm
