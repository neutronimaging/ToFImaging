import os
from qtpy.QtWidgets import QMainWindow, QVBoxLayout
from jupyter_notebooks.code import load_ui
import pyqtgraph as pg
import numpy as np


class Interface(QMainWindow):

    def __init__(self, parent=None, working_dir="", sample_projections=None, spectra_file=None):
        super(Interface, self).__init__(parent)

        self.sample_projections = sample_projections

        ui_full_path = os.path.join(os.path.dirname(__file__), os.path.join('ui', 'ui_roi_selection.ui'))
        self.ui = load_ui(ui_full_path, baseinstance=self)
        self.setWindowTitle("Region of Interest Tool")

        self.init_widgets()
        self.display_image()

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

    def display_image(self):
        sample_projections = self.sample_projections
        live_image = np.transpose(np.mean(sample_projections, axis=2))
        self.ui.image_view.setImage(live_image)

    def update_table_roi(self):
        pass

    def add_roi_button_clicked(self):
        pass

    def remove_roi_button_clicked(self):
        pass

    def cancel_clicked(self):
        self.close()

    def apply_clicked(self):
        print('create mask!')
