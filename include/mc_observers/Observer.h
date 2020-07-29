/*
 * Copyright 2015-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

/*! Interface used to load observers */

#include <mc_observers/api.h>

#include <mc_rtc/gui/StateBuilder.h>
#include <mc_rtc/log/Logger.h>

namespace mc_control
{
struct MCController;
} // namespace mc_control

namespace mc_rbdyn
{

struct Robots;
}

namespace mc_observers
{

/**
 * @brief State observation API
 *
 * All new observers must inherit from this Observer class, and implement the
 * required virtual functions (at least reset, run, updateRobots)
 */
struct MC_OBSERVERS_DLLAPI Observer
{
  virtual ~Observer() {}

  /**
   * @brief Configure observer
   *
   * @param ctl Controller running the observer
   * @param config Configuration for this observer
   */
  virtual void configure(const mc_control::MCController & /*ctl*/, const mc_rtc::Configuration & /*config*/) {}

  /*! \brief Reset estimator.
   *
   * \param ctl Controller reference (const). Used to access information about
   * the robot and the controller (anchor frame, contacts, desired contact
   * force, etc).
   */
  virtual void reset(const mc_control::MCController & ctl) = 0;

  /*! \brief Compute observer state
   *
   * \param ctl Controller reference (const). Used to access information about
   * the robot and the controller (anchor frame, contacts, desired contact
   * force, etc).
   */
  virtual bool run(const mc_control::MCController & ctl) = 0;

  /*! \brief Update the real robot state from the observed state
   *
   * \param ctl Controller reference (const). Used to access information about
   * the robot and the controller (anchor frame, contacts, desired contact
   * force, etc).
   *
   * \param robot Robot state to write to. Each controller is expected to update
   * the real robot instance with its estimates. The pipeline will only call the
   * updateRobots() function if requested by the user.
   */
  virtual void updateRobots(mc_control::MCController & ctl) = 0;

  /**
   * @brief Set the observer's name
   *
   * @param name Name of the observer
   */
  void name(const std::string & name)
  {
    name_ = name;
  }

  /**
   * @brief Returns the observer's name
   */
  const std::string & name() const;

  /*! \brief Add observer to the logger.
   *
   * Default implementation does nothing, each observer implementation is
   * responsible for logging its own data by overriding this function
   */
  virtual void addToLogger(mc_control::MCController &, const std::string & /* category */ = "") {}
  /*! \brief Remove observer from logger
   *
   * Default implementation does nothing, each observer implementation is
   * responsible for removing all logs entry that it added.
   */
  virtual void removeFromLogger(mc_control::MCController &, const std::string & /* category */ = "") {}
  /*! \brief Add observer information the GUI.
   *
   * Default implementation does nothing, each observer implementation is
   * responsible for adding its own elements to the GUI. Default observers will
   * be shown under the tab "Observers->observer name".
   */
  virtual void addToGUI(mc_rtc::gui::StateBuilder &, std::vector<std::string> /* category */ = {}) {}
  /*! \brief Remove observer from gui
   *
   * Default implementation removes the category {category, observer name}
   */
  virtual void removeFromGUI(mc_rtc::gui::StateBuilder &, std::vector<std::string> /* category */ = {});

  /*! \brief Short description of the observer
   *
   * Used to display a short summary of the observers pipeline to the user, the description should be as short as
   * possible while showing the important information about the observer.
   *
   * The output might end up looking like:
   * Observers: Encoders (Position+Velocity) -> KinematicInertial (cutoff=0.01) -> BodySensor
   *
   * \returns Short description of the observer. Default implementation returns
   * its name.
   */
  virtual const std::string & desc() const;

  virtual const std::string type() const = 0;

protected:
  std::string name_;

  /* Short descriptive description of the observer used for CLI logging */
  std::string desc_;
};

using ObserverPtr = std::shared_ptr<mc_observers::Observer>;

} // namespace mc_observers

#ifdef WIN32
#  define OBSERVER_MODULE_API __declspec(dllexport)
#else
#  if __GNUC__ >= 4
#    define OBSERVER_MODULE_API __attribute__((visibility("default")))
#  else
#    define OBSERVER_MODULE_API
#  endif
#endif
