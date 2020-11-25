import os
from qtpy.QtWidgets import QMainWindow, QVBoxLayout, QTableWidgetItem, QTableWidgetSelectionRange
from qtpy.QtGui import QColor, QPen
from jupyter_notebooks.code import load_ui
import pyqtgraph as pg
import numpy as np
from collections import OrderedDict

MARKER_HEIGHT, MARKER_WIDTH = 20, 20


class Interface(QMainWindow):

    o_roi = None
    o_api = None
    live_image = None

    pixel_marker = {'x': 0,
                    'y': 0}
    cross_of_pixel_to_fit = None


    # list_roi = OrderedDict()
    # default_roi = {'x0': 0, 'y0': 0, 'x1': 50, 'y1': 50, 'id': None}
    # live_image = None
    # roi_width = 0.01
    # integrated_image_size = {'width': -1, 'height': -1}
    # splitter_2_has_been_resized = False

    def __init__(self, parent=None, o_roi=None, o_api=None):
        super(Interface, self).__init__(parent)

        self.o_roi = o_roi
        self.o_api = o_api

        # self.sample_projections = main_api.sample_projections  # x, y, lambda
        # self.sample_projections_lambda_x_y = self.sample_projections.transpose(2, 0, 1)  # lambda, x, y
        # self.tof_array = main_api.tof_array
        # self.lambda_array = main_api.lambda_array

        ui_full_path = os.path.join(os.path.dirname(__file__), os.path.join('ui', 'ui_fit.ui'))
        self.ui = load_ui(ui_full_path, baseinstance=self)
        self.setWindowTitle("Fit Interface")

        self.init_widgets()
        self.initialize_pixel_marker()
        self.display_image()
        self.display_profile()
        self.display_roi()
        self.display_cross_of_pixel_to_fit()
        self.display_box_around_pixel_to_fit()

    def initialize_pixel_marker(self):
        for _index_roi in self.o_roi.list_roi:
            _roi = self.o_roi.list_roi[_index_roi]
            _x0 = _roi['x0']
            _y0 = _roi['y0']
            _x1 = _roi['x1']
            _y1 = _roi['y1']

            x = np.mean([_x0, _x1])
            y = np.mean([_y0, _y1])

            self.pixel_marker = {'x': x,
                                 'y': y}

    def init_widgets(self):
        # pyqtgraphs
        self.ui.image_view = pg.ImageView(view=pg.PlotItem())
        self.ui.image_view.ui.roiBtn.hide()
        self.ui.image_view.ui.menuBtn.hide()
        verti_layout1 = QVBoxLayout()
        verti_layout1.addWidget(self.ui.image_view)
        self.ui.widget_image.setLayout(verti_layout1)

        self.ui.plot_view = pg.PlotWidget()
        verti_layout2 = QVBoxLayout()
        verti_layout2.addWidget(self.ui.plot_view)
        self.ui.widget_plot.setLayout(verti_layout2)

        # splitters
        self.ui.splitter.setSizes([550, 50])
        self.ui.splitter_2.setSizes([200, 2])

    def pixel_marker_changed(self):
        region = self.pixel_marker_item.getArraySlice(self.live_image,
                                                      self.ui.image_view.imageItem)

        x0 = region[0][0].start
        y0 = region[0][1].start
        self.pixel_marker['x'] = x0
        self.pixel_marker['y'] = y0

        self.display_cross_of_pixel_to_fit()
        

    def display_box_around_pixel_to_fit(self):
        x, y = self.pixel_marker['x'], self.pixel_marker['y']

        self.pixel_marker_item = pg.ROI([x, y],
                                   [MARKER_WIDTH, MARKER_HEIGHT],
                                   scaleSnap=True)
        self.ui.image_view.addItem(self.pixel_marker_item)
        self.pixel_marker_item.sigRegionChanged.connect(self.pixel_marker_changed)

    def display_cross_of_pixel_to_fit(self):

        if self.cross_of_pixel_to_fit:
            self.ui.image_view.removeItem(self.cross_of_pixel_to_fit)

        x, y = self.pixel_marker['x'], self.pixel_marker['y']

        pos = []
        adj = []

        # vertical guide
        pos.append([x + MARKER_WIDTH / 2, y - MARKER_HEIGHT / 2])
        pos.append([x + MARKER_WIDTH / 2, y + MARKER_HEIGHT + MARKER_HEIGHT / 2])
        adj.append([0, 1])

        # horizontal guide
        pos.append([x - MARKER_WIDTH / 2, y + MARKER_HEIGHT / 2])
        pos.append([x + MARKER_WIDTH + MARKER_WIDTH / 2, y + MARKER_HEIGHT / 2])
        adj.append([2, 3])

        pos = np.array(pos)
        adj = np.array(adj)

        line_color = (255, 0, 0, 255, 1)
        lines = np.array([line_color for _ in np.arange(len(pos))],
                         dtype=[('red', np.ubyte), ('green', np.ubyte),
                                ('blue', np.ubyte), ('alpha', np.ubyte),
                                ('width', float)])
        self.cross_of_pixel_to_fit = pg.GraphItem()
        self.ui.image_view.addItem(self.cross_of_pixel_to_fit)
        self.cross_of_pixel_to_fit.setData(pos=pos,
                                  adj=adj,
                                  pen=lines,
                                  symbol=None,
                                  pxMode=False)

    def display_roi(self):
        for _index_roi in self.o_roi.list_roi:
            _roi = self.o_roi.list_roi[_index_roi]
            _id = _roi['id']

            pos = _id.pos()
            size = _id.size()

            new_id = pg.ROI(pos, size, movable=False)
            self.ui.image_view.addItem(new_id)

    def display_image(self):
        self.live_image = self.o_roi.live_image
        live_image = np.transpose(self.live_image)
        self.ui.image_view.setImage(live_image)

    def display_profile(self):
        normalize_projections_lambda_x_y = self.o_api.normalize_projections.transpose(2, 0, 1)
        list_roi = self.o_roi.list_roi

        profile_y_axis = []
        total_number_of_pixels_in_rois = 0
        for _index_roi in list_roi.keys():
            _roi = list_roi[_index_roi]
            _x0 = _roi['x0']
            _y0 = _roi['y0']
            _x1 = _roi['x1']
            _y1 = _roi['y1']
            total_number_of_pixels_in_rois += (_y1 - _y0 + 1) * (_x1 - _x0 + 1)

        for _projection in normalize_projections_lambda_x_y:

            total_counts_of_roi = 0
            for _index_roi in list_roi.keys():
                _roi = list_roi[_index_roi]
                _x0 = _roi['x0']
                _y0 = _roi['y0']
                _x1 = _roi['x1']
                _y1 = _roi['y1']
                total_counts_of_roi += np.sum(_projection[_y0: _y1 + 1, _x0: _x1 + 1])

            mean_counts_of_roi = total_counts_of_roi / total_number_of_pixels_in_rois
            profile_y_axis.append(mean_counts_of_roi)

        self.profile_y_axis = profile_y_axis

        x_axis = self.o_roi.lambda_array
        self.ui.plot_view.plot(x_axis, profile_y_axis)

    def apply_clicked(self):
        self.close()
