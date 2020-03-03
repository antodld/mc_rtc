#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Copyright 2015-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

from PyQt5 import QtCore, QtGui, QtWidgets

from PyQt5.QtWidgets import QWidget, QVBoxLayout

import copy
import math
import matplotlib
import numpy as np
matplotlib.use('Qt5Agg')
import matplotlib.pyplot

from matplotlib.animation import FuncAnimation
from matplotlib.figure import Figure
from matplotlib.patches import Patch, Polygon, Rectangle

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas,\
                                               NavigationToolbar2QT as NavigationToolbar

from collections import OrderedDict
from math import asin, atan2

from mc_log_types import LineStyle, PlotSide


def rpyFromMat(E):
    """Same as mc_rbdyn::rpyFromMat."""
    roll = atan2(E[1][2], E[2][2]);
    pitch = -asin(E[0][2]);
    yaw = atan2(E[0][1], E[0][0]);
    return [roll, pitch, yaw]


def rpyFromQuat(quat):
    """Same as mc_rbdyn::rpyFromQuat."""
    import eigen
    return rpyFromMat(list(eigen.Quaterniond(*quat).toRotationMatrix()))

class PlotPolygonAxis(object):
  def __init__(self, parent, axis):
    self.figure = parent
    self._axis = axis.twinx()
    self._axis.set_zorder(-1)
    self._axis.set_yticks([])
    self._axis.set_ylim(0, 1)
    self._axis.get_yaxis().set_visible(False)
    self.plotsCount = {}
    self.plots = OrderedDict()
    self.colors = OrderedDict()
  def _plot_string(self, x, y, y_label, style):
    if y_label in self.plots:
      self.plotsCount[y_label] += 1
      return False
    self.plotsCount[y_label] = 1
    self.plots[y_label] = []
    i = 0
    i0 = 0
    label = y[i0]
    while i < len(y):
      if y[i] == label and i + 1 != len(y):
        i += 1
        continue
      if len(label) == 0:
        i += 1
        continue
      if label not in self.colors:
        self.colors[label] = self.figure._next_poly_color()
      if i + 1 < len(y):
        xi = x[i + 1]
      else:
        xi = x[i]
      self.plots[y_label].append(self._axis.add_patch(Rectangle((x[i0], 0), xi - x[i0], 1, label = label, facecolor = self.colors[label])))
      i += 1
      i0 = i
      if i < len(y):
        label = y[i0]
    return True
  def legend(self):
    if not len(self.plots):
      self._axis.clear()
      return
    xL = 1.0
    if len(self.figure._right()):
      xL = 1.025
    self._axis.legend(bbox_to_anchor=(xL, 0., 0.2, 0.1), mode="expand", borderaxespad=0.5)
  def remove_plot(self, y):
    if y not in self.plots:
      return
    self.plotsCount[y] -= 1
    if self.plotsCount[y]:
      return
    for plt in self.plots[y]:
      plt.remove()
    del self.plots[y]
    self.figure.draw()

class PlotYAxis(object):
  def __init__(self, parent, x_axis = None, poly = None):
    self.figure = parent
    if x_axis is None:
      self.is_left = True
      self._axis = parent.fig.add_subplot(111)
      self._axis.autoscale(enable = True, axis = 'both', tight = False)
      self._x_axis = self._axis
    else:
      self.is_left = False
      self._axis = x_axis.twinx()
      self._x_axis = x_axis
    if poly is None:
      self._polyAxis = PlotPolygonAxis(parent, self._axis)
    else:
      self._polyAxis = poly
    self._axis.set_facecolor((1, 1, 1, 0))
    box = self._axis.get_position()
    self._axis.autoscale_view(False,True,True)
    self._axis.format_coord = parent.format_coord
    self._axis.get_yaxis().set_visible(False)
    self.grid = LineStyle(linestyle = '--')
    self.plots = OrderedDict()
    self._label_fontsize = 10
    self._legend_ncol = 3
    self.data = {}

  def __len__(self):
    return len(self.plots)

  def _data(self):
    return self.figure.data

  def _legend_fontsize(self):
    return self.figure._legend_fontsize

  def _labelpad(self):
    return self.figure._labelpad

  def _x_label_fontsize(self):
    return self.figure._x_label_fontsize

  def axis(self):
    return self._axis

  def drawGrid(self):
    if len(self.plots):
      self._axis.grid(color = self.grid.color, linestyle = self.grid.linestyle, linewidth = self.grid.linewidth, visible = self.grid.visible, which = 'both')

  def legend(self):
    if not len(self.plots):
      return
    if self.is_left:
      loc = 3
      top_anchor = 1.02
    else:
      loc = 2
      top_anchor = -0.14
      if len(self.figure.x_label()):
        top_anchor = -0.175
    self._axis.legend(bbox_to_anchor=(0., top_anchor, 1., .102), loc=loc, ncol=self._legend_ncol, mode="expand", borderaxespad=0.5, fontsize=self._legend_fontsize())
    self._polyAxis.legend()

  def legendNCol(self, n = None):
    if n is None:
      return self._legend_ncol
    self._legend_ncol = n
    self.legend()

  def legendRows(self):
    return math.ceil(len(self.plots) / float(self._legend_ncol))

  def legendOffset(self, offset, sign):
    if self.legendRows() > 3:
      offset = offset + sign * 0.035 * (self.legendRows() - 3)
    return offset

  # If this axis is empty but the other is not, fix my own limits, in the
  # opposite condition call the opposite function
  def fixLimits(self, axis):
    if len(axis) != 0 and len(self) == 0:
      point = (axis.axis().dataLim.get_points()[1] + axis.axis().dataLim.get_points()[0])/2
      plt, = self._axis.plot([point[0]], [point[1]], visible = False)
      plt.remove()
      del plt
      self._axis.relim()
    elif len(self) != 0 and len(axis) == 0:
      axis.fixLimits(self)

  def getLimits(self, frame, idx):
    if not len(self):
      return None
    min_ = np.nanmin(self.data.values()[0][idx][:frame])
    max_ = np.nanmax(self.data.values()[0][idx][:frame])
    for i in range(1, len(self.data.values())):
      data = self.data.values()[i][idx][:frame]
      min_ = min(np.nanmin(data), min_)
      max_ = max(np.nanmax(data), max_)
    return min_, max_


  def setLimits(self, xlim = None, ylim = None, frame = None):
    if not len(self):
      return xlim
    dataLim = self._axis.dataLim.get_points()
    def setLimit(lim, idx, set_lim):
      if lim is not None:
        min_, max_ = lim
      elif frame is not None:
        min_, max_ = self.getLimits(frame, idx)
      else:
        range_ = dataLim[1][idx] - dataLim[0][idx]
        min_ = dataLim[0][idx] - range_ * 0.01
        max_ = dataLim[1][idx] + range_ * 0.01
      set_lim([min_, max_])
      return min_, max_
    setLimit(ylim, 1, self._axis.set_ylim)
    return setLimit(xlim, 0, self._x_axis.set_xlim)

  def _label(self, get_label, set_label, l, size):
    if l is None:
      l = get_label()
    set_label(l, fontsize = size, labelpad = self._labelpad())

  def _label_property(self, get_label, set_label, l = None):
    if l is None:
      return get_label()
    set_label(l)

  def _x_label(self, l = None):
    self._label(self.x_label, self._x_axis.set_xlabel, l, self._x_label_fontsize())

  def x_label(self, l = None):
    return self._label_property(self._x_axis.get_xlabel, self._x_label, l)

  def _y_label(self, l = None):
    self._label(self.y_label, self._axis.set_ylabel, l, self._label_fontsize)

  def y_label(self, l = None):
    return self._label_property(self._axis.get_ylabel, self._y_label, l)

  def y_label_fontsize(self, fontsize = None):
    if fontsize is None:
      return self._label_fontsize
    self._label_fontsize = fontsize
    self._y_label()

  def animate(self, frame):
    for y_label in self.plots.keys():
      self.plots[y_label].set_data(self.data[y_label][0][:frame], self.data[y_label][1][:frame])
    return self.plots.values()

  def _plot(self, x, y, y_label, style = None):
    if type(y[0]) is unicode:
      return self._polyAxis._plot_string(x, y, y_label, style)
    if style is None:
      return self._plot(x, y, y_label, LineStyle(color = self.figure._next_color()))
    if y_label in self.plots:
      return False
    self._axis.get_yaxis().set_visible(True)
    self.plots[y_label] = self._axis.plot(x, y, label = y_label, color = style.color, linestyle = style.linestyle, linewidth = style.linewidth)[0]
    self.data[y_label] = (x, y)
    self.legend()
    return True

  def startAnimation(self, i0):
    for y_label in self.plots.keys():
      style = self.style(y_label)
      self.plots[y_label] = self._axis.plot(self.data[y_label][0][i0], self.data[y_label][1][i0], label = y_label, color = style.color, linestyle = style.linestyle, linewidth = style.linewidth)[0]

  def stopAnimation(self):
    for y_label in self.plots.keys():
      style = self.style(y_label)
      self.plots[y_label].remove()
      self.plots[y_label] = self._axis.plot(self.data[y_label][0], self.data[y_label][1], label = y_label, color = style.color, linestyle = style.linestyle, linewidth = style.linewidth)[0]

  def add_plot(self, x, y, y_label, style = None):
    return self._plot(self._data()[x], self._data()[y], y_label, style)

  def add_diff_plot(self, x, y, y_label):
    dt = self._data()[x][1] - self._data()[x][0]
    return self._plot(self._data()[x][1:], np.diff(self._data()[y])//dt, y_label)

  def _add_rpy_plot(self, x_label, y, idx):
    assert (idx >= 0 and idx <= 2),"index must be 0, 1 or 2"
    rpy_label = ['roll', 'pitch', 'yaw']
    y_label = "{}_{}".format(y, rpy_label[idx])
    fmt = ""
    if "{}_qw".format(y) in self._data().keys():
      fmt = "q"
    qw, qx, qy, qz = [ self._data()[k] for k in [ "{}_{}{}".format(y, fmt, ax) for ax in ["w", "x", "y", "z"] ] ]
    data = [ rpyFromQuat([w, x, y, z])[idx] for w, x, y, z in zip(qw, qx, qy, qz) ]
    return self._plot(self._data()[x_label], data, y_label)

  def add_roll_plot(self, x, y):
    return self._add_rpy_plot(x, y, 0)

  def add_pitch_plot(self, x, y):
    return self._add_rpy_plot(x, y, 1)

  def add_yaw_plot(self, x, y):
    return self._add_rpy_plot(x, y, 2)

  def add_rpy_plot(self, x, y):
    r = self.add_roll_plot(x, y)
    p = self.add_pitch_plot(x, y)
    y = self.add_yaw_plot(x, y)
    return r or p or y

  def remove_plot(self, y):
    if y not in self.plots:
      self._polyAxis.remove_plot(y)
      return
    self.plots[y].remove()
    del self.plots[y]
    del self.data[y]
    if len(self.plots):
      self._axis.relim()
      self.legend()
    else:
      self._axis.get_yaxis().set_visible(False)
      self._axis.clear()

  def clear(self):
    self.plots = {}
    self._axis.clear()

  # Get or set the style of a given plot
  def style(self, y, style = None):
    if y not in self.plots:
      raise KeyError("No plot named {}".format(y))
    plt = self.plots[y]
    if style is None:
      return LineStyle(plt.get_color(), plt.get_linestyle(), plt.get_linewidth(), label = plt.get_label())
    plt.set_color(style.color)
    plt.set_linestyle(style.linestyle)
    plt.set_linewidth(style.linewidth)
    if len(style.label):
      plt.set_label(style.label)

class PlotFigure(object):
  def __init__(self):
    self.fig = matplotlib.pyplot.figure(figsize=(5, 4), dpi=100)
    self.axes = {}
    self.axes[PlotSide.LEFT] = PlotYAxis(self)
    self.axes[PlotSide.RIGHT] = PlotYAxis(self, self._left().axis(), self._left()._polyAxis)
    self.animation = None

    self._title_fontsize = 12
    self._x_label_fontsize = 10
    self._labelpad = 10
    self._tick_labelsize = 10
    self._legend_fontsize = 10
    self._top_offset = 0.9
    self._bottom_offset = 0.1

    self.data = None
    self.computed_data = {}

    self.color = 0
    cm = matplotlib.cm.Set1
    self.Ncolor = min(cm.N, 12)
    cm2rgb = (np.array(cm(x)[0:3]) for x in np.linspace(0, 1, self.Ncolor))
    self.colors = ['#%02x%02x%02x' % tuple((255 * rgb).astype(int)) for rgb in cm2rgb]

    self.polyColor = 0
    cm = matplotlib.cm.Pastel1
    cm2rgb = (np.array(cm(x)[0:3] + (0.5,)) for x in np.linspace(0, 1, min(cm.N, 32)))
    self.polyColors = [rgb for rgb in cm2rgb]

    self.x_data = 't'

  # Helper function to call something on all axes, call expects an axis argument
  def _axes(self, call):
    call(self._left())
    call(self._right())

  # Shortcut to the left axis
  def _left(self):
    return self.axes[PlotSide.LEFT]

  # Shortcut to the right axis
  def _right(self):
    return self.axes[PlotSide.RIGHT]

  def _drawGrid(self):
    self._axes(lambda axis: axis.drawGrid())

  def _legend(self):
    self._axes(lambda axis: axis.legend())

  def draw(self, x_limits = None, y1_limits = None, y2_limits = None, frame = None):
    self._left().fixLimits(self._right())
    x_limits = self._left().setLimits(x_limits, y1_limits, frame = frame)
    self._right().setLimits(x_limits, y2_limits, frame = frame)
    self._legend()
    self._drawGrid()
    top_offset = self._left().legendOffset(self._top_offset, -1)
    bottom_offset = self._right().legendOffset(self._bottom_offset, 1)
    left_offset = 0.125
    right_offset = 0.9
    if len(self._left()._polyAxis.plots):
      left_offset, right_offset = 0.05, right_offset - (left_offset - 0.05)
    self.fig.subplots_adjust(left = left_offset, right = right_offset, top = top_offset, bottom = bottom_offset)

  def animate(self, frame, x_limits = None, y1_limits = None, y2_limits = None):
    ret = self._left().animate(frame)
    ret.extend(self._right().animate(frame))
    PlotFigure.draw(self, x_limits = x_limits, y1_limits = y1_limits, y2_limits = y2_limits, frame = frame)
    return ret

  def stopAnimation(self):
    self._axes(lambda a: a.stopAnimation())

  def setData(self, data):
    self.data = data

  def setColors(self, colors):
    self.colors = colors
    self.Ncolor = len(self.colors)

  def show(self):
    self.fig.show()

  def top_offset(self, off = None):
    if off is None:
      return self._top_offset
    else:
      self._top_offset = off

  def bottom_offset(self, off = None):
    if off is None:
      return self._bottom_offset
    else:
      self._bottom_offset = off

  def title(self, title = None):
    if title is None:
      if self.fig._suptitle is None:
        return ""
      return self.fig._suptitle.get_text()
    self.fig.suptitle(title)

  def title_fontsize(self, fontsize = None):
    if fontsize is None:
      return self._title_fontsize
    self._title_fontsize = fontsize
    self.fig.suptitle(self.title(), fontsize = self._title_fontsize)

  def tick_fontsize(self, size = None):
    if size is None:
      return self._tick_labelsize
    self._tick_labelsize = size
    self._axes(lambda a: a.axis().tick_params(labelsize = self._tick_labelsize))

  def labelpad(self, pad = None):
    if pad is None:
      return self._labelpad
    self._labelpad = pad
    self._left()._x_label()
    self._axes(lambda axis: axis._y_label())

  def _x_label(self, label = None):
    self._left()._x_label(label)

  def _y1_label(self, label = None):
    self._left()._y_label(label)

  def _y2_label(self, label = None):
    self._right()._y_label(label)

  def x_label(self, label = None):
    return self._left().x_label(label)

  def x_label_fontsize(self, fontsize = None):
    if fontsize is None:
      return self._x_label_fontsize
    self._x_label_fontsize = fontsize
    self._x_label()

  def y1_label(self, label = None):
    return self._left().y_label(label)

  def y1_label_fontsize(self, fontsize = None):
    return self._left().y_label_fontsize(fontsize)

  def y2_label(self, label = None):
    return self._right().y_label(label)

  def y2_label_fontsize(self, fontsize = None):
    return self._right().y_label_fontsize(fontsize)

  def _next_poly_color(self):
    self.polyColor += 1
    return self.polyColors[ (self.polyColor - 1) % len(self.polyColors) ]

  def _next_color(self):
    self.color += 1
    return self.colors[ (self.color - 1) % self.Ncolor ]

  def legend_fontsize(self, size = None):
    if size is None:
      return self._legend_fontsize
    self._legend_fontsize = size
    self._legend()

  def y1_legend_ncol(self, n = None):
    return self._left().legendNCol(n)

  def y2_legend_ncol(self, n = None):
    return self._right().legendNCol(n)

  def add_plot_left(self, x, y, y_label, style = None):
    return self._left().add_plot(x, y, y_label, style)

  def add_plot_right(self, x, y, y_label, style = None):
    return self._right().add_plot(x, y, y_label, style)

  def add_diff_plot_left(self, x, y, y_label):
    return self._left().add_diff_plot(x, y, y_label)

  def add_diff_plot_right(self, x, y, y_label):
    return self._right().add_diff_plot(x, y, y_label)

  def add_roll_plot_left(self, x, y):
    return self._left().add_roll_plot(x, y)

  def add_pitch_plot_left(self, x, y):
    return self._left().add_pitch_plot(x, y)

  def add_yaw_plot_left(self, x, y):
    return self._left().add_yaw_plot(x, y)

  def add_roll_plot_right(self, x, y):
    return self._right().add_roll_plot(x, y)

  def add_pitch_plot_right(self, x, y):
    return self._right().add_pitch_plot(x, y)

  def add_yaw_plot_right(self, x, y):
    return self._right().add_yaw_plot(x, y)

  def add_rpy_plot_left(self, x, y):
    return self._left().add_rpy_plot(x, y)

  def add_rpy_plot_right(self, x, y):
    return self._right().add_rpy_plot(x, y)

  def _remove_plot(self, SIDE, y_label):
    self.axes[SIDE].remove_plot(y_label)
    if len(self._left()) == 0 and len(self._right()) == 0:
      self.color = 0

  def remove_plot_left(self, y_label):
    self._remove_plot(PlotSide.LEFT, y_label)

  def remove_plot_right(self, y_label):
    self._remove_plot(PlotSide.RIGHT, y_label)

  def format_coord(self, x, y):
    display_coord = self.axes[PlotSide.RIGHT].axis().transData.transform((x,y))
    inv = self.axes[PlotSide.LEFT].axis().transData.inverted()
    ax_coord = inv.transform(display_coord)
    if len(self._left()) and len(self._right()):
      return "x: {:.3f}    y1: {:.3f}    y2: {:.3f}".format(x, ax_coord[1], y)
    elif len(self._left()):
      return "x: {:.3f}    y1: {:.3f}".format(x, ax_coord[1])
    elif len(self._right()):
      return "x: {:.3f}    y2: {:.3f}".format(x, y)
    else:
      return "x: {:.3f}".format(x)

  def clear_all(self):
    self.color = 0
    self._axes(lambda a: a.clear())

  def style_left(self, y, styleIn = None):
    return self._left().style(y, styleIn)

  def style_right(self, y, styleIn = None):
    return self._right().style(y, styleIn)

class SimpleAxesDialog(QtWidgets.QDialog):
  def __init__(self, parent):
    QtWidgets.QDialog.__init__(self, parent)
    self.setWindowTitle('Edit axes limits')
    self.setModal(True)
    self.layout = QtWidgets.QGridLayout(self)
    self.layout.addWidget(QtWidgets.QLabel("Min"), 0, 1)
    self.layout.addWidget(QtWidgets.QLabel("Max"), 0, 2)

    self.layout.addWidget(QtWidgets.QLabel("X"), 1, 0)
    x_limits = parent.x_limits
    if x_limits is None:
      x_limits = parent._left().axis().get_xlim()
    self.x_min = QtWidgets.QLineEdit(str(x_limits[0]))
    self.x_min.setValidator(QtGui.QDoubleValidator())
    self.layout.addWidget(self.x_min, 1, 1)
    self.x_max = QtWidgets.QLineEdit(str(x_limits[1]))
    self.x_max.setValidator(QtGui.QDoubleValidator())
    self.x_init = [float(self.x_min.text()), float(self.x_max.text())]
    self.layout.addWidget(self.x_max, 1, 2)

    self.layout.addWidget(QtWidgets.QLabel("Y1"), 2, 0)
    y1_limits = parent.y1_limits
    if y1_limits is None:
      y1_limits = parent._left().axis().get_ylim()
    self.y1_min = QtWidgets.QLineEdit(str(y1_limits[0]))
    self.y1_min.setValidator(QtGui.QDoubleValidator())
    self.layout.addWidget(self.y1_min, 2, 1)
    self.y1_max = QtWidgets.QLineEdit(str(y1_limits[1]))
    self.y1_max.setValidator(QtGui.QDoubleValidator())
    self.y1_init = [float(self.y1_min.text()), float(self.y1_max.text())]
    self.layout.addWidget(self.y1_max, 2, 2)

    self.layout.addWidget(QtWidgets.QLabel("Y2"), 3, 0)
    y2_limits = parent.y2_limits
    if y2_limits is None:
      y2_limits = parent._right().axis().get_ylim()
    self.y2_min = QtWidgets.QLineEdit(str(y2_limits[0]))
    self.y2_min.setValidator(QtGui.QDoubleValidator())
    self.layout.addWidget(self.y2_min, 3, 1)
    self.y2_max = QtWidgets.QLineEdit(str(y2_limits[1]))
    self.y2_max.setValidator(QtGui.QDoubleValidator())
    self.y2_init = [float(self.y2_min.text()), float(self.y2_max.text())]
    self.layout.addWidget(self.y2_max, 3, 2)

    confirmLayout = QtWidgets.QHBoxLayout()
    okButton = QtWidgets.QPushButton("Ok", self)
    confirmLayout.addWidget(okButton)
    okButton.clicked.connect(self.accept)
    applyButton = QtWidgets.QPushButton("Apply", self)
    confirmLayout.addWidget(applyButton)
    applyButton.clicked.connect(self.apply)
    cancelButton = QtWidgets.QPushButton("Cancel", self)
    confirmLayout.addWidget(cancelButton)
    cancelButton.clicked.connect(self.reject)
    self.layout.addLayout(confirmLayout, 4, 0, 1, 3)

  def apply(self):
    changed = False
    x_limits = [float(self.x_min.text()), float(self.x_max.text())]
    if x_limits != self.x_init:
      changed = True
      self.parent().x_locked.setChecked(True)
      self.parent().x_limits = x_limits
    y1_limits = [float(self.y1_min.text()), float(self.y1_max.text())]
    if y1_limits != self.y1_init:
      changed = True
      self.parent().y1_locked.setChecked(True)
      self.parent().y1_limits = y1_limits
    y2_limits = [float(self.y2_min.text()), float(self.y2_max.text())]
    if y2_limits != self.y2_init:
      changed = True
      self.parent().y2_locked.setChecked(True)
      self.parent().y2_limits = y2_limits
    if changed:
      self.parent().draw()

  def accept(self):
    QtWidgets.QDialog.accept(self)
    self.apply()

class PlotCanvasWithToolbar(PlotFigure, QWidget):
  def __init__(self, parent = None):
    PlotFigure.__init__(self)
    QWidget.__init__(self, parent)

    self.canvas = FigureCanvas(self.fig)
    self.canvas.mpl_connect('draw_event', self.on_draw)
    self.toolbar = NavigationToolbar(self.canvas, self)

    vbox = QVBoxLayout(self)
    vbox.addWidget(self.canvas)
    vbox.addWidget(self.toolbar)

  def setupLockButtons(self, layout):
    self.x_locked = QtWidgets.QPushButton(u"🔒X", self)
    self.x_locked.setCheckable(True)
    layout.addWidget(self.x_locked)
    self.x_locked.toggled.connect(self.x_locked_changed)
    self.x_limits = None

    self.y1_locked = QtWidgets.QPushButton(u"🔒 Y1", self)
    self.y1_locked.setCheckable(True)
    layout.addWidget(self.y1_locked)
    self.y1_locked.toggled.connect(self.y1_locked_changed)
    self.y1_limits = None

    self.y2_locked = QtWidgets.QPushButton(u"🔒 Y2", self)
    self.y2_locked.setCheckable(True)
    layout.addWidget(self.y2_locked)
    self.y2_locked.toggled.connect(self.y2_locked_changed)
    self.y2_limits = None

  def setupAnimationButtons(self, layout):
    animationLayout = QtWidgets.QHBoxLayout()
    self.animation = None
    self.animationButton = QtWidgets.QPushButton("Start animation")
    self.animationButton.setCheckable(True)
    self.animationButton.toggled.connect(self.startStopAnimation)
    animationLayout.addWidget(self.animationButton)
    self.lockAxesButton = QtWidgets.QPushButton("Lock axes")
    self.lockAxesButton.released.connect(self.lockAxes)
    animationLayout.addWidget(self.lockAxesButton)
    self.saveAnimationButton = QtWidgets.QPushButton("Save animation")
    self.saveAnimationButton.released.connect(self.saveAnimation)
    animationLayout.addWidget(self.saveAnimationButton)
    layout.addLayout(animationLayout)

  def startStopAnimation(self):
    if self.animationButton.isChecked():
      if self.startAnimation():
        self.animationButton.setText("Stop animation")
      else:
        self.animationButton.setChecked(False)
    else:
      self.stopAnimation()
      self.animationButton.setText("Start animation")

  def restartAnimation(self):
    if self.animationButton.isChecked():
      self.stopAnimation()
      self.startAnimation()

  def getFrameRange(self):
    if self.data is None or len(self.data) == 0:
      return 0, 0
    x_data = self.data[self.x_data]
    i0 = 0
    while i0 < len(x_data) and np.isnan(x_data[i0]):
      i0 += 1
    iN = i0
    while iN + 1 < len(x_data) and not np.isnan(x_data[iN + 1]):
      iN += 1
    assert(iN > i0 and i0 < len(x_data)),"Strange time range"
    return i0, iN

  def lockAxes(self):
    i0, iN = self.getFrameRange()
    if i0 == iN:
      return
    self.x_limits = self._left().getLimits(iN, 0) or self._right().getLimits(iN, 0)
    if self.x_limits is None:
      return
    self.x_locked.setChecked(True)
    self.y1_limits = self._left().getLimits(iN, 1)
    if self.y1_limits is not None:
      self.y1_locked.setChecked(True)
    self.y2_limits = self._right().getLimits(iN, 1)
    if self.y2_limits is not None:
      self.y2_locked.setChecked(True)

  def startAnimation(self):
    interval = 50 # ms
    i0, iN = self.getFrameRange()
    if i0 == iN:
      return False
    x_data = self.data[self.x_data]
    dt = (x_data[i0 + 1] - x_data[i0]) * 1000 # dt in ms
    step = int(math.ceil(interval/dt))
    self._axes(lambda axis: axis._axis.clear())
    self.animation = FuncAnimation(self.fig, self.animate, frames = range(i0 + 1, iN, step), interval = interval)
    self._axes(lambda a: a.startAnimation(i0))
    self.draw()
    return True

  def animate(self, frame):
    return PlotFigure.animate(self, frame, self.x_limits, self.y1_limits, self.y2_limits)

  def stopAnimation(self):
    self.animation.event_source.stop()
    PlotFigure.stopAnimation(self)
    self.draw()

  def saveAnimation(self):
    fpath = QtWidgets.QFileDialog.getSaveFileName(self, "Output file", filter = "Video (*.mp4)")[0]
    if not len(fpath):
      return
    if self.animationButton.isChecked():
      self.animation.save(fpath)
    else:
      self.startAnimation()
      self.animation.save(fpath)
      self.stopAnimation()

  def axesDialog(self):
    SimpleAxesDialog(self).exec_()

  def on_draw(self, event):
    if self.x_limits is not None:
      self.x_limits = self._left().axis().get_xlim()
    if self.y1_limits is not None:
      self.y1_limits = self._left().axis().get_ylim()
    if self.y2_limits is not None:
      self.y2_limits = self._right().axis().get_ylim()

  def draw(self):
    PlotFigure.draw(self, self.x_limits, self.y1_limits, self.y2_limits)
    self.canvas.draw()

  def _y_lock_changed(self, name, cbox, get_lim):
    if cbox.isChecked():
      cbox.setText(u"🔓 {}".format(name))
      return get_lim()
    else:
      cbox.setText(u"🔒{}".format(name))
      return None

  def x_locked_changed(self, status):
    self.x_limits = self._y_lock_changed("X", self.x_locked, self._left().axis().get_xlim)
    if self.x_limits is None:
      self.draw()

  def y1_locked_changed(self, status):
    self.y1_limits = self._y_lock_changed("Y1", self.y1_locked, self._left().axis().get_ylim)
    if self.y1_limits is None:
      self.draw()

  def y2_locked_changed(self, status):
    self.y2_limits = self._y_lock_changed("Y2", self.y2_locked, self._right().axis().get_ylim)
    if self.y2_limits is None:
      self.draw()
