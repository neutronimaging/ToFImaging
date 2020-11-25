import os
from qtpy.QtWidgets import QMainWindow, QVBoxLayout, QTableWidgetItem, QTableWidgetSelectionRange
from qtpy.QtGui import QColor, QPen
from jupyter_notebooks.code import load_ui
import pyqtgraph as pg
import numpy as np
from collections import OrderedDict


class Interface(QMainWindow):

    o_roi = None
    o_api = None
    live_image = None

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
        self.display_image()
        self.display_profile()
        self.display_roi()

    def init_widgets(self):
        # pyqtgraphs
        self.ui.image_view = pg.ImageView()
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
