/*
 * Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <mc_tasks/TrajectoryTaskGeneric.h>

namespace mc_tasks
{

/*! \brief Control a robot's CoM */
struct MC_TASKS_DLLAPI CoMTask6D : public TrajectoryTaskGeneric
{
public:
  /*! \brief Constructor
   *
   * \param robots Robots involved in the task
   *
   * \param robotIndex Select the robot which CoM should be controlled
   *
   * \param stiffness Task stiffness
   *
   * \param weight Task weight
   *
   */
  CoMTask6D(const mc_rbdyn::Robots & robots, unsigned int robotIndex, double stiffness = 5.0, double weight = 100.);

  void reset() override;

  /*! \brief Change the CoM target by a given amount
   *
   * \param com Modification applied to the CoM
   *
   */
  void move_com(const sva::PTransformd & com);

  /*! \brief Set the CoM target to a given position
   *
   * \param com New CoM target
   *
   */
  void com(const sva::PTransformd & com);

  /*! \brief Return the current CoM target
   *
   * \returns The current CoM target
   */
  const sva::PTransformd & com() const;

  /** \brief Actual CoM position (computed at the previous iteration)
   *
   * \return Return the current CoM position
   */
  const sva::PTransformd & actual() const;

  /*! \brief Load from configuration */
  void load(mc_solver::QPSolver &, const mc_rtc::Configuration & config) override;

  void flight(const bool s);
  bool flight();

protected:
  void addToGUI(mc_rtc::gui::StateBuilder &) override;
  void addToLogger(mc_rtc::Logger & logger) override;
  void update(mc_solver::QPSolver &) override;

private:
  unsigned int robot_index_;
  bool flight_ = false;
};

} // namespace mc_tasks
