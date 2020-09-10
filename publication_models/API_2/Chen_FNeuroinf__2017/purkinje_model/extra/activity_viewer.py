from __future__ import print_function
from __future__ import unicode_literals

import numpy as np
from pylab import *
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui, QtOpenGL
import pyqtgraph.opengl as gl
import random
import sys
import os
import random
from numpy import outer
from matplotlib.backends import qt_compat
from matplotlib import cm
import numpy as np
from PyQt5 import QtGui, QtCore
from numpy import arange, sin, pi
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

def getJetColor(v, min_v, max_v):
    if v < min_v:
        v = min_v
    if v > max_v:
        v = max_v
    gray_scale = (v - min_v) / (max_v - min_v)
    return cm.jet(gray_scale)

def readData(data_file):
    dataset = {}
    file = open(data_file, 'r')
    for line in file:
        line_secs = line.split()
        # entry info
        if line_secs[0] == '#Entries:':
            dataset["Entries"] = line_secs[1:]
            dataset["Data"] = []
        else:
            data_line = [float(value) for value in line_secs]
            dataset["Data"].append(data_line)
    return dataset
    
def SI2NEURON(dataset):
    new_dataset = {}
    new_dataset["Entries"] = dataset["Entries"]
    new_dataset["Data"] = []
    for data_line in dataset["Data"]:
        new_data_line = []
        # WARNING: Using a variable name that is reserved (['time']).
        time = data_line[0] * 1000
        # WARNING: Using a variable name that is reserved (['time']).
        new_data_line.append(time)
        for v in data_line[1:]:
            new_data_line.append(v*1e6)
        new_dataset["Data"].append(new_data_line)
    return new_dataset

class ColorFigCanvas(FigureCanvas):
    def __init__(self, parent=None, min_v = 0.0, max_v = 1.0, width = 80.0, height= 600.0):
        fig = Figure()
        dpi = fig.get_dpi()
        fig.set_size_inches(width/float(dpi),height/float(dpi))
        self.axes = fig.add_subplot(111, frame_on = True)
        # We want the axes cleared every time plot() is called
        # self.axes.hold(False)
        
        self.compute_initial_figure(min_v, max_v)
        
        #
        # WARNING: Using a variable name that is reserved (['__init__']).
        FigureCanvas.__init__(self, fig)
        self.setParent(parent)
        #FigureCanvas.setMinimumSize(self, 100, 2)
        self.resize(width, height)
        FigureCanvas.setSizePolicy(self,
                                   QtGui.QSizePolicy.Minimum,
                                   QtGui.QSizePolicy.Minimum)
        FigureCanvas.updateGeometry(self)
    
    def compute_initial_figure(self, min_v, max_v):
        a=np.outer(arange(min_v, max_v, (max_v - min_v) / 100.0), np.ones(1))
        self.axes.imshow(a,cmap=cm.jet,origin="lower", extent=[ 0.0,1.0, min_v,max_v], aspect = 80.0 / (max_v - min_v))
        self.axes.axes.get_xaxis().set_visible(False)


class TriActivitySeriesDisplay(QtGui.QMainWindow):
    """
    Partition Display
    Parameters:
        * mesh                STEPS Tetmesh
        * tri_partitions      Partition of triangles
        * title               Display title
        * x                   X cordinate of the display
        * y                   Y cordinate of the display
        * w                   Width of the display
        * h                   Height of the display
        * scale               Scaling between STEPS mesurement and display pixel
        * series_color_map           Color map for each partition, it is a Python dict with key-value pairs as color_map[partition_id] = [red, green, blue, alpha], where partition_id is the data stored in tet_partitions. The color is defined by four parameters [red, green, blue, alpha], each with range from 0.0 to 1.0. If the partition name is not in the color map a random color will be generated. A partition assigned to None will always be colored as nontransparent red [1.0, 0.0, 0.0, 1.0].
    """
    
    def __init__(self, mesh, tri_partitions, activity_data, min_v = None, max_v = None, fps = 60, title = "TriActivitySeriesDisplay", x = 100, y = 100, w = 800, h = 600, data_unit = "", time_unit = "", scale = 1e6):
        """
        Constructor.
        """
        # WARNING: Using a variable name that is reserved (['__init__']).
        QtGui.QMainWindow.__init__(self)
        self.setGeometry(x, y, w, h)
        
        # filter data
        self.min_v = min_v
        self.max_v = max_v
        
        if min_v == None:
            self.min_v = float('Inf')
        if max_v == None:
            self.max_v = -float('Inf')
        
        self.series_data = activity_data["Data"]
        self.roi_entries = activity_data["Entries"][1:]
        self.n_entries = len(self.roi_entries)
        
        for tpn_data in self.series_data:
            for v in tpn_data[1:]:
                if self.min_v > v:
                    self.min_v = v
                if self.max_v < v:
                    self.max_v = v
        if min_v != None:
            self.min_v = min_v
        if max_v != None:
            self.max_v = max_v

        self.curr_tpn = 0
        self.ntpns = len(self.series_data)
        
        self.setGeometry(x, y, w, h)
        center_widget = QtGui.QSplitter(self)
        main_layout = QtGui.QHBoxLayout()
        
        self.view_widget = gl.GLViewWidget(center_widget)
        self.view_widget.resize(w - 50.0, h)
        self.view_widget.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        
        main_layout.addWidget(self.view_widget)
        info_panel = QtGui.QWidget(center_widget)
        info_panel.resize(50.0, h)
        info_panel.setSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
                                   
        info_panel_layout = QtGui.QVBoxLayout()
        
        label = QtGui.QLabel("Unit: %s" % (data_unit), info_panel)
        label.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignVCenter)
        info_panel_layout.addWidget(label)
        
        info_panel_layout.addWidget(ColorFigCanvas(info_panel, self.min_v, self.max_v, height = (h - 100)))
        
        self.time_unit = time_unit
        self.time_text = QtGui.QLabel("%f %s" % (self.series_data[self.curr_tpn][0], time_unit), info_panel)
        self.time_text.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignVCenter)
        info_panel_layout.addWidget(self.time_text)
        
        self.play_btn = QtGui.QPushButton("Play", info_panel)
        info_panel_layout.addWidget(self.play_btn)
        self.play_btn.clicked.connect(self.playBtnClicked)
        self.inplay = False
        
        self.reset_btn = QtGui.QPushButton("Reset", info_panel)
        info_panel_layout.addWidget(self.reset_btn)
        self.reset_btn.clicked.connect(self.resetBtnClicked)
        
        label = QtGui.QLabel("FPS:")
        label.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignVCenter)
        info_panel_layout.addWidget(label)
        
        self.fps = fps
        self.fps_edit = QtGui.QLineEdit(info_panel)
        validator = QtGui.QDoubleValidator()
        self.fps_edit.setValidator(validator)
        self.fps_edit.setText("%i" % (self.fps))
        self.fps_edit.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignVCenter)
        self.fps_edit.textChanged.connect(self.fpsChanged)
        self.fps_edit.textChanged.emit(self.fps_edit.text())
        info_panel_layout.addWidget(self.fps_edit)
        
        info_panel_layout.addStretch(1)
    
        info_panel.setLayout(info_panel_layout)
        main_layout.addWidget(info_panel)
        
        center_widget.setLayout(main_layout)
        
        self.setCentralWidget(center_widget)

        self.mesh = mesh
        
        self.tri_part_table = {}
        for tri in tri_partitions:
            if tri_partitions[tri] not in self.tri_part_table.keys():
                self.tri_part_table[tri_partitions[tri]] = []
            self.tri_part_table[tri_partitions[tri]].append(tri)
        
        self.scale = scale
        self.setWindowTitle(title)
        
        self.bound_min = [v * self.scale for v in mesh.bbox.min]
        self.bound_max = [v * self.scale for v in mesh.bbox.max]
        
        self.center = [0.0, 0.0, 0.0]
        
        axis = [0, 1, 2]
        new_center = [(self.bound_min[i] + self.bound_max[i]) / 2.0 for i in axis]
        
        pan_dist = [new_center[i] - self.center[i] for i in axis]
        self.view_widget.pan(pan_dist[0], pan_dist[1], pan_dist[2])
        
        dist_x = self.bound_max[0] - self.bound_min[0]
        dist_y = self.bound_max[1] - self.bound_min[1]
        dist_z = self.bound_max[2] - self.bound_min[2]
        
                
        if dist_y >= dist_x and dist_y >= dist_z:
            self.main_axis = 1
            self.view_widget.setCameraPosition(distance=dist_y * 2)
        
        elif dist_z >= dist_x and dist_z >= dist_y:
            self.view_widget.setCameraPosition(distance=dist_z * 2)
            self.main_axis = 2

        else:
            self.view_widget.setCameraPosition(distance=dist_x * 2)
            self.main_axis = 0

        self.draw_parts = []
        self.center = new_center

        for r in range(self.n_entries):
            # WARNING: Using a variable name that is reserved (['r']).
            color = getJetColor(self.series_data[self.curr_tpn][r+1], self.min_v, self.max_v)
            # WARNING: Using a variable name that is reserved (['r']).
            part = TriPartitionMesh(self, mesh.stepsMesh, self.tri_part_table[self.roi_entries[r]], color = color)
            self.draw_parts.append(part)
            self.view_widget.addItem(part)
        self.show()

        self.timer = QtCore.QTimer()
        
        if self.ntpns == 1:
            self.play_btn.setEnabled(False)
            self.reset_btn.setEnabled(False)
        
    def playBtnClicked(self):
        if self.inplay:
            self.inplay = False
            self.play_btn.setText("Play")
        else:
            self.inplay = True
            self.play_btn.setText("Stop")
            self.nextTpn()

    def resetBtnClicked(self):
        self.play_btn.setEnabled(True)
        self.curr_tpn = 0
        self.time_text.setText("%f %s" % (self.series_data[self.curr_tpn][0], self.time_unit))
        for r in range(self.n_entries):
            # WARNING: Using a variable name that is reserved (['r']).
            color = getJetColor(self.series_data[self.curr_tpn][r+1], self.min_v, self.max_v)
            # WARNING: Using a variable name that is reserved (['r']).
            self.draw_parts[r].setColor(color)
    
    def nextTpn(self):
        if not self.inplay:
            return
        self.curr_tpn += 1
        if self.curr_tpn % self.fps == 0:
            self.time_text.setText("%f %s" % (self.series_data[self.curr_tpn][0], self.time_unit))
        for r in range(self.n_entries):
            # WARNING: Using a variable name that is reserved (['r']).
            color = getJetColor(self.series_data[self.curr_tpn][r+1], self.min_v, self.max_v)
            #print(self.roi_entries[r], " ", self.series_data[self.curr_tpn][r+1], " ", color)
            # WARNING: Using a variable name that is reserved (['r']).
            self.draw_parts[r].setColor(color)
        if self.curr_tpn < self.ntpns - 1:
            self.timer.singleShot(1000 / self.fps, self.nextTpn)
        else:
            self.play_btn.setEnabled(False)
            self.play_btn.setText("Play")
            self.inplay = False

    def fpsChanged(self, *args, **kwargs):
        sender = self.sender()
        validator = sender.validator()
        state = validator.validate(sender.text(), 0)[0]
        if state == QtGui.QValidator.Acceptable:
            color = '#c4df9b' # green
            self.fps = float(sender.text())
        elif state == QtGui.QValidator.Intermediate:
            color = '#fff79a' # yellow
        else:
            color = '#f6989d' # red
        sender.setStyleSheet('QLineEdit { background-color: %s }' % color)

class TriPartitionMesh(gl.GLMeshItem):
    """
    Static mesh component for a compartment
    Parameters:
        * display                 Parent display
        * steps_mesh              STEPS mesh
        * tri_list                Triangle list of a section of the mesh
        * color                   Color of the component
    """
    def __init__(self, display, mesh, tri_list, color):
        """
        Constructor
        """
        self.display = display
        self.mesh = mesh
        
        surface_tris = np.array(tri_list, dtype = np.uint32)
        
        v_set_size = mesh.getTriVerticesSetSizeNP(surface_tris)
        tris_data = np.zeros(surface_tris.size * 3, dtype = np.uint32)
        v_set = np.zeros(v_set_size, dtype = np.uint32)
        verts_data = np.zeros(v_set_size * 3)
        mesh.getTriVerticesMappingSetNP(surface_tris, tris_data, v_set)
        mesh.getBatchVerticesNP(v_set, verts_data)
        verts_data *= display.scale
        tris_data.shape = -1, 3
        verts_data.shape = -1, 3
        mesh_data = gl.MeshData(vertexes=verts_data, faces = tris_data)
        # WARNING: Using a variable name that is reserved (['__init__']).
        gl.GLMeshItem.__init__(self, meshdata=mesh_data, smooth=False, computeNormals =True, shader='balloon', glOptions='additive')
        self.setColor(color)



