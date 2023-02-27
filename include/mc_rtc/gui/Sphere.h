#pragma once

#include <mc_rtc/gui/Visual.h>
#include <mc_rtc/visual_utils.h>

namespace mc_rtc::gui
{

/** Creates a Sphere
 *
 * \tparam GetPos Callback to get the sphere position
 *
 * \tparam GetRadius A double (fixed radius) or callback to get the radius
 *
 * \tparam GetColor A color (fixed color) or callback to get the color
 *
 * If \tparam GetRadius and/or \tparam GetColor are callbacks they are invoked immediately
 */
template<typename GetPos, typename GetRadius = double, typename GetColor = const mc_rtc::gui::Color &>
auto Sphere(const std::string & name,
            GetRadius radius_fn,
            GetPos get_pos_fn,
            GetColor color_fn = mc_rtc::gui::Color::Red)
{
  auto sphere =
      mc_rtc::makeVisualSphere(details::GetValueOrCallbackValue(radius_fn), details::GetValueOrCallbackValue(color_fn));
  std::function<rbd::parsers::Visual &()> get_visual_fn = [sphere]() mutable -> rbd::parsers::Visual & {
    return sphere;
  };
  if constexpr(std::is_invocable_v<GetRadius>)
  {
    get_visual_fn = [get_visual_fn, radius_fn]() -> rbd::parsers::Visual & {
      auto & sphere = get_visual_fn();
      mc_rtc::getVisualSphere(sphere).radius = radius_fn();
      return sphere;
    };
  }
  if constexpr(std::is_invocable_v<GetColor>)
  {
    get_visual_fn = [get_visual_fn, color_fn]() -> rbd::parsers::Visual & {
      auto & visual = get_visual_fn();
      mc_rtc::details::setVisualColor(visual, color_fn());
      return visual;
    };
  }
  return Visual(name, get_visual_fn, get_pos_fn);
}

} // namespace mc_rtc::gui
