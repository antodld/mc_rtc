#pragma once

#include <mc_control/fsm/State.h>

#include <mc_control/CompletionCriteria.h>
#include <mc_tasks/MetaTask.h>

namespace mc_control
{

namespace fsm
{

/** Implements a generic tasks state
 *
 * This tasks creates an arbitrary number of MetaTasks and executes them
 * until completion criteria are fullfiled for each tasks.
 *
 * The tasks in the state are configured through the "tasks" entry which a
 * JSON object.
 *
 * The keys are the name of the tasks and values are MetaTask objects as
 * expected by MetaTaskLoader plus an optional "completion" criteria that
 * represents one or more completion criteria.
 *
 * The names of the tasks are only relevant for the state. You can override
 * the actual task's name using the name entry in the Task's configuration.
 *
 * If the "completion" entry is absent and if there is no existing
 * completion entry for the related task then this task is added but not
 * considered as part of the completion criteria (e.g. you can add a
 * CoMTask and an EndEffectorTask but only care for the completion of the
 * later).
 *
 * When the "tasks" entry is read multiple times, the following ensues:
 * - if a new tasks appears then it is added
 * - if an existing task is not repeated, nothing happens for this task
 * - if a task already exists, existing configuration entries are
 *   overwriten by the new entry, non-existing configuration entries are
 *   simply added to the existing configuration
 * Example:
 *
 * \code{.json}
 * # First pass, we simplify task entries for the sake of the example
 * {
 *   "tasks":
 *   {
 *      "t1":
 *      {
 *        "objectiveA": 0.5,
 *        "objectiveB": 1.0,
 *        "completion": { "timeout": 5.0 }
 *      }
 *   }
 * }
 * # After this pass, one task is considered
 * # Second pass
 * {
 *   "tasks":
 *   {
 *     "t1":
 *     {
 *       "objectiveA": 1.0,
 *       "completion": { "eval": 1e-6 }
 *     },
 *     "t2":
 *     {
 *       "objective": 0.5
 *     }
 *   }
 * }
 * # We now have two tasks, and:
 * # - t1's objectiveA is changed to 1.0, objectiveB is the same
 * # - t1 completion criteria is replaced
 * # Third pass
 * {
 *   "tasks":
 *   {
 *     "t1":
 *     {
 *       "completion": {}
 *     },
 *     "t2":
 *     {
 *       "completion": { "eval": 1e-6 }
 *     }
 *   }
 * }
 * # We still have two tasks, objectives are unchanged but:
 * - t1 has no more completion criteria
 * - t2 has a completion criteria
 * \endcode
 *
 * Using the AddContacts and RemoveContacts entries, one can specify that
 * should be added/removed at the start of the state.
 *
 * If RemovePostureTask is set to true, the posture task is removed at the
 * start of the state.
 *
 */

struct MC_CONTROL_DLLAPI MetaTasksState : State
{
  void configure(const mc_rtc::Configuration & config) override;

  void start(Controller&) override;

  bool run(Controller&) override;

  void teardown(Controller&) override;
protected:
  std::map<std::string, mc_rtc::Configuration> tasks_configs_;
  std::vector<mc_tasks::MetaTaskPtr> tasks_;
  std::vector<std::pair<size_t, mc_control::CompletionCriteria>> criterias_;
  bool finished_first_ = false;
};

} // namespace fsm

} // namespace mc_control
