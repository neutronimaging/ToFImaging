import os
from qtpy.QtWidgets import QMainWindow, QVBoxLayout, QTableWidgetItem, QTableWidgetSelectionRange
from qtpy.QtGui import QColor, QPen
from jupyter_notebooks.code import load_ui
import pyqtgraph as pg
import numpy as np
from collections import OrderedDict


class Interface(QMainWindow):

    # list_roi = OrderedDict()
    # default_roi = {'x0': 0, 'y0': 0, 'x1': 50, 'y1': 50, 'id': None}
    # live_image = None
    # roi_width = 0.01
    # integrated_image_size = {'width': -1, 'height': -1}
    # splitter_2_has_been_resized = False

    def __init__(self, parent=None):
        super(Interface, self).__init__(parent)

        # self.sample_projections = main_api.sample_projections  # x, y, lambda
        # self.sample_projections_lambda_x_y = self.sample_projections.transpose(2, 0, 1)  # lambda, x, y
        # self.tof_array = main_api.tof_array
        # self.lambda_array = main_api.lambda_array

        ui_full_path = os.path.join(os.path.dirname(__file__), os.path.join('ui', 'ui_fit.ui'))
        self.ui = load_ui(ui_full_path, baseinstance=self)
        self.setWindowTitle("Fit Interface")

        # self.init_widgets()
        # self.display_image()

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

        # table
        nbr_columns = self.ui.table_roi.columnCount()
        for _col in range(nbr_columns):
            self.ui.table_roi.setColumnWidth(_col, 70)

        # splitters
        self.ui.splitter_2.setSizes([400, 0])

        # x_axis buttons
        self.ui.tof_radioButton.setText(u"TOF (\u03BCs)")
        self.ui.lambda_radioButton.setText(u"\u03BB (\u212B)")

    def display_image(self):
        sample_projections = self.sample_projections
        self.live_image = np.mean(sample_projections, axis=2)
        live_image = np.transpose(self.live_image)
        self.ui.image_view.setImage(live_image)
        self.integrated_image_size['height'], self.integrated_image_size['width'] = np.shape(live_image)

    def apply_clicked(self):
        self.close()
