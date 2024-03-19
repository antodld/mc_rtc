/*
 * Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <mc_tvm/api.h>
#include <mc_tvm/fwd.h>

#include <mc_rtc/shared.h>

#include <RBDyn/CoM.h>

#include <tvm/graph/abstract/Node.h>

namespace mc_tvm
{

/** Center of mass (CoM) of a Robot and related quantities
 *
 * Provides the frame position, jacobian, velocity and normal acceleration.
 * These signals are correctly initialized on the object's creation.
 *
 * Outputs:
 * - Position: position of the CoM in world coordinates
 * - Jacobian: jacobian of the CoM in world coordinates
 * - Velocity: velocity of the CoM in world coordinates
 * - NormalAcceleration: normal acceleration of the CoM in world coordinates
 * - Acceleration: acceleration of the CoM in world coordinates
 * - JDot: derivative of the jacobian of the CoM in world coordinates
 *
 */
struct MC_TVM_DLLAPI CoM6D : public tvm::graph::abstract::Node<CoM6D>
{
  SET_OUTPUTS(CoM6D, CoM6D, Jacobian, Velocity, NormalAcceleration, Acceleration, JDot)
  SET_UPDATES(CoM6D, CoM6D, Jacobian, Velocity, NormalAcceleration, Acceleration, JDot)

  friend struct Robot;

private:
  struct NewCoMToken
  {
  };

public:
  /** Constructor
   *
   * Creates the CoM algorithm for a robot
   *
   * \param robot Robot to which the frame is attached
   *
   */
  CoM6D(NewCoMToken, Robot & robot);

  inline const sva::PTransformd & com() const noexcept { return com_; }

  inline const Eigen::Vector6d & velocity() const noexcept { return velocity_; }

  inline const Eigen::Vector6d & normalAcceleration() const noexcept { return normalAcceleration_; }

  inline const Eigen::Vector6d & acceleration() const noexcept { return acceleration_; }

  inline const Eigen::MatrixXd & jacobian() const noexcept { return jac_; }

  inline const Eigen::MatrixXd & JDot() const noexcept { return jacDot_; }

  inline const Robot & robot() const noexcept { return robot_; }

  inline Robot & robot() noexcept { return robot_; }


private:
  Robot & robot_;

  Eigen::MatrixXd jac_;
  Eigen::MatrixXd jacDot_;

  sva::PTransformd com_;
  void updateCoM();

  Eigen::Vector6d velocity_;
  void updateVelocity();

  Eigen::Vector6d normalAcceleration_;
  void updateNormalAcceleration();

  Eigen::Vector6d acceleration_;
  void updateAcceleration();

  void updateJacobian();

  void updateJDot();
};

} // namespace mc_tvm
