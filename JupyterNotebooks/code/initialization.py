from qtpy.QtWidgets import QVBoxLayout
import pyqtgraph as pg
from pyqtgraph.dockarea import *
from qtpy.QtGui import QIcon
import numpy as np
import os

from JupyterNotebooks.code.config_handler import ConfigHandler


class Initialization:

    def __init__(self, parent=None):
        self.parent = parent

        # load config
        o_config = ConfigHandler(parent=self.parent)
        o_config.load()

    def widgets(self):

        # disable all tabs except first one (prepare data)
        self.parent.ui.toolBox.setItemEnabled(1, False)
        self.parent.ui.toolBox.setItemEnabled(2, False)

        # labels
        self.parent.ui.kernel_size_custom_lambda_label.setText(u"\u03BB:")
        self.parent.kernel_dimension_changed()

        # hide normalization if not needed
        if self.parent.debugging_mode:
            normalization_flag_value = True
        else:
            normalization_flag_value = self.parent.o_api.normalization_flag_ui.value
        self.parent.ui.prepare_data_normalization_groupBox.setVisible(normalization_flag_value)

        # prepare data tab
        self.parent.ui.raw_image_view = pg.ImageView(view=pg.PlotItem(), name='raw image')
        self.parent.ui.raw_image_view.ui.roiBtn.hide()
        self.parent.ui.raw_image_view.ui.menuBtn.hide()
        prepare_data_vertical_layout1 = QVBoxLayout()
        prepare_data_vertical_layout1.addWidget(self.parent.ui.raw_image_view)
        self.parent.ui.prepare_data_raw_widget.setLayout(prepare_data_vertical_layout1)

        self.parent.ui.process_image_view = pg.ImageView(view=pg.PlotItem(), name='process image')
        self.parent.ui.process_image_view.ui.roiBtn.hide()
        self.parent.ui.process_image_view.ui.menuBtn.hide()
        prepare_data_vertical_layout2 = QVBoxLayout()
        prepare_data_vertical_layout2.addWidget(self.parent.ui.process_image_view)
        self.parent.ui.prepare_data_process_widget.setLayout(prepare_data_vertical_layout2)

        # self.parent.ui.process_image_view.view.getViewBox().setXLink('raw image')
        # self.parent.ui.process_image_view.view.getViewBox().setYLink('raw image')

        # fit tab
        self.parent.ui.image_view = pg.ImageView(view=pg.PlotItem())
        self.parent.ui.image_view.ui.roiBtn.hide()
        self.parent.ui.image_view.ui.menuBtn.hide()
        vertical_layout1 = QVBoxLayout()
        vertical_layout1.addWidget(self.parent.ui.image_view)
        self.parent.ui.widget_image.setLayout(vertical_layout1)

        self.parent.ui.plot_view = pg.PlotWidget()
        vertical_layout2 = QVBoxLayout()
        vertical_layout2.addWidget(self.parent.ui.plot_view)
        self.parent.ui.widget_plot.setLayout(vertical_layout2)

        # splitters
        self.parent.ui.splitter.setSizes([550, 50])
        self.parent.ui.splitter_2.setSizes([200, 2])

        self.parent.ui.prepare_data_main_splitter.setSizes([100, 500])
        self.parent.ui.prepare_data_preview_splitter.setSizes([200, 200])

        # sliders
        if self.parent.debugging_mode:
            half_number_of_files = 100
        else:
            half_number_of_files = np.int(len(self.parent.o_roi.lambda_array) / 2)
        self.parent.ui.right_number_of_files_to_exclude_slider.setMaximum(half_number_of_files)
        self.parent.ui.left_number_of_files_to_exclude_slider.setMaximum(half_number_of_files)
        self.parent.ui.right_number_of_files_to_exclude_slider.setValue(
                self.parent.nbr_files_to_exclude_from_plot['left'])
        self.parent.ui.left_number_of_files_to_exclude_slider.setValue(
                self.parent.nbr_files_to_exclude_from_plot['right'])

        # icons
        _file_path = os.path.dirname(__file__)
        settings_icon_1 = os.path.abspath(os.path.join(_file_path, 'static/settings_icon.png'))
        self.parent.ui.step3_fit_pixel_settings_pushButton.setIcon(QIcon(settings_icon_1))
        settings_icon_2 = os.path.abspath(os.path.join(_file_path, 'static/settings_icon.png'))
        self.parent.ui.step4_fit_roi_settings_pushButton.setIcon(QIcon(settings_icon_2))

        # result tab

        # gaussian fitting
        area = DockArea()
        result_vertical_layout = QVBoxLayout()
        result_vertical_layout.addWidget(area)
        self.parent.ui.result_widget.setLayout(result_vertical_layout)

        dock_input_image, self.parent.ui.dock_input_image_view = \
            Initialization._set_up_dock(dock_name="Input Image")
        dock_edge_position, self.parent.ui.dock_edge_position_view = \
            Initialization._set_up_dock(dock_name="Edge Position")
        dock_edge_height, self.parent.ui.dock_edge_height_view = \
            Initialization._set_up_dock(dock_name="Edge Height")
        dock_edge_width, self.parent.ui.dock_edge_width_view = \
            Initialization._set_up_dock(dock_name="Edge Width")
        dock_edge_slope, self.parent.ui.dock_edge_slope_view = \
            Initialization._set_up_dock(dock_name="Edge Slope")
        dock_median_image, self.parent.ui.dock_median_image_view = \
            Initialization._set_up_dock(dock_name="Image Median")

        self.parent.ui.dock_input_image_view.view.getViewBox().setXLink("Edge Position View")
        self.parent.ui.dock_input_image_view.view.getViewBox().setYLink("Edge Position View")

        self.parent.ui.dock_edge_position_view.view.getViewBox().setXLink("Edge Height View")
        self.parent.ui.dock_edge_position_view.view.getViewBox().setYLink("Edge Height View")

        self.parent.ui.dock_edge_height_view.view.getViewBox().setXLink("Edge Width View")
        self.parent.ui.dock_edge_height_view.view.getViewBox().setYLink("Edge Width View")

        self.parent.ui.dock_edge_width_view.view.getViewBox().setXLink("Edge Slope View")
        self.parent.ui.dock_edge_width_view.view.getViewBox().setYLink("Edge Slope View")

        self.parent.ui.dock_edge_slope_view.view.getViewBox().setXLink("Image Median View")
        self.parent.ui.dock_edge_slope_view.view.getViewBox().setYLink("Image Median View")

        area.addDock(dock_input_image, 'left')
        area.addDock(dock_edge_position, 'right', dock_input_image)
        area.addDock(dock_edge_height, 'right', dock_edge_position)
        area.addDock(dock_edge_width, 'bottom', dock_input_image)
        area.addDock(dock_edge_slope, 'bottom', dock_edge_position)
        area.addDock(dock_median_image, 'bottom', dock_edge_height)

        # advanced fitting
        area1 = DockArea()
        result_vertical_layout1 = QVBoxLayout()
        result_vertical_layout1.addWidget(area1)
        self.parent.ui.advanced_result_widget.setLayout(result_vertical_layout1)

        dock_input_image1, self.parent.ui.advanced_dock_input_image_view = \
            Initialization._set_up_dock(dock_name="Input Image")
        dock_edge_position1, self.parent.ui.advanced_dock_edge_position_view = \
            Initialization._set_up_dock(dock_name="Edge Position")
        dock_edge_height1, self.parent.ui.advanced_dock_edge_height_view = \
            Initialization._set_up_dock(dock_name="Edge Height")
        dock_edge_width1, self.parent.ui.advanced_dock_edge_width_view = \
            Initialization._set_up_dock(dock_name="Edge Width")

        self.parent.ui.advanced_dock_input_image_view.view.getViewBox().setXLink("Edge Position View")
        self.parent.ui.advanced_dock_input_image_view.view.getViewBox().setYLink("Edge Position View")

        self.parent.ui.advanced_dock_edge_position_view.view.getViewBox().setXLink("Edge Height View")
        self.parent.ui.advanced_dock_edge_position_view.view.getViewBox().setYLink("Edge Height View")

        self.parent.ui.advanced_dock_edge_height_view.view.getViewBox().setXLink("Edge Width View")
        self.parent.ui.advanced_dock_edge_height_view.view.getViewBox().setYLink("Edge Width View")

        area1.addDock(dock_input_image1, 'left')
        area1.addDock(dock_edge_position1, 'right', dock_input_image1)
        area1.addDock(dock_edge_height1, 'right', dock_edge_position1)
        area1.addDock(dock_edge_width1, 'bottom', dock_input_image1)

    @staticmethod
    def _set_up_dock(size=(300, 300), dock_name=''):
        dock_input_image = Dock(dock_name, size=size)
        dock_input_image_view = pg.ImageView(view=pg.PlotItem(), name=dock_name + " View")
        dock_input_image_view.ui.roiBtn.hide()
        dock_input_image_view.ui.menuBtn.hide()
        dock_input_image.addWidget(dock_input_image_view)
        return dock_input_image, dock_input_image_view
