/*
 * Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <mc_tvm/CoM6D.h>

#include <mc_rbdyn/fwd.h>

#include <tvm/function/abstract/Function.h>

namespace mc_tvm
{

/** This class implements a CoM function for a given robot
 *
 * You can provide:
 * - a reference CoM target
 * - a reference CoM velocity target
 * - a reference CoM acceleration target
 */
struct MC_TVM_DLLAPI CoM6DFunction : tvm::function::abstract::Function
{
  SET_UPDATES(CoM6DFunction, Value, Velocity, Jacobian, NormalAcceleration, JDot)

  /** Constructor
   *
   * Creates a CoM function, the objective is the current robot's CoM
   */
  CoM6DFunction(const mc_rbdyn::Robot & robot);

  /** Reset the objective to the current CoM and reference speed/acceleration to zero */
  void reset();

  /** Get the current objective */
  inline const sva::PTransformd & com() const noexcept { return com_; }

  /** Set the objective */
  inline void com(const sva::PTransformd & com) noexcept { com_ = com; }

  /** Get the robot's current CoM */
  inline const sva::PTransformd & actual() const noexcept { return comAlgo_.com(); }

  /** Get the current objective */
  inline const Eigen::Vector6d & refVel() const noexcept { return refVel_; }

  /** Set the objective */
  inline void refVel(const Eigen::Vector6d & refVel) noexcept { refVel_ = refVel; }

  /** Get the current objective */
  inline const Eigen::Vector6d & refAccel() const noexcept { return refAccel_; }

  /** Set the objective */
  inline void refAccel(const Eigen::Vector6d & refAccel) noexcept { refAccel_ = refAccel; }

private:
  const mc_tvm::CoM6D & comAlgo_;
  //com objective
  sva::PTransformd com_;
  Eigen::Vector6d refVel_;
  Eigen::Vector6d refAccel_;

  void updateValue();
  void updateVelocity();
  void updateJacobian();
  void updateNormalAcceleration();
  void updateJDot();
};

} // namespace mc_tvm
