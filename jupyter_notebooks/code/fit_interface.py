import os
from qtpy.QtWidgets import QMainWindow
from jupyter_notebooks.code import load_ui
import numpy as np
import copy
import inflect
import logging
import versioneer
import warnings

from jupyter_notebooks.code.decorators import wait_cursor
from jupyter_notebooks.code.fit_handler import FitHandler
from jupyter_notebooks.code.utilities.get import Get
from jupyter_notebooks.code.initialization import Initialization
from jupyter_notebooks.code.export import Export
from jupyter_notebooks.code.import_data import Import
from jupyter_notebooks.code.event_handler import EventHandler
from jupyter_notebooks.code.prepare_data import PrepareData
from jupyter_notebooks.code.step3_settings_handler import Step3SettingsHandler
from jupyter_notebooks.code.step4_settings_handler import Step4SettingsHandler
from jupyter_notebooks.code.display import Display
from jupyter_notebooks.code.log_launcher import LogLauncher
from jupyter_notebooks.code.gui_handler import GuiHandler

# warnings.filterwarnings('ignore')

MARKER_HEIGHT, MARKER_WIDTH = 20, 20
COLOR_LAMBDA_RANGE = [250, 128, 247]
COLOR_ROUGH_LAMBDA = [50, 255, 50]


class Interface(QMainWindow):

    config = None
    debugging_mode = False

    o_roi = None
    o_api = None
    live_image = None

    histogram_level = {"process_image_view": [],
                       "image_view": []}

    list_roi = None
    ob_list_roi = None

    # pyqtgraph items in fitting tab
    roi_id = None
    bragg_peak_range_ui = None
    cross_of_pixel_to_fit = None
    pixel_marker_item = None
    rough_peak_ui = None

    moving_average_config = None

    is_with_normalization = False

    sample_projections = None
    ob_projections = None

    default_kernel_size = {'y': 3, 'x': 3, 'lambda': 3}
    kernel_size = copy.deepcopy(default_kernel_size)
    default_kernel_size_label = {'3d': u"y:3  x:3  \u03BB:3",
                                 '2d': "y:3  x:3"}

    pixel_marker = {'x': 0,
                    'y': 0}
    cross_of_pixel_to_fit = None

    default_bragg_peak_range = None

    nbr_files_to_exclude_from_plot = {'left': 40,
                                      'right': 40}
    rough_peak_index_position = 0
    rough_peak_ui = None

    pixel_fit_result = None
    full_fit_result = None

    # list_roi = OrderedDict()
    # default_roi = {'x0': 0, 'y0': 0, 'x1': 50, 'y1': 50, 'id': None}
    # live_image = None
    # roi_width = 0.01
    # integrated_image_size = {'width': -1, 'height': -1}
    # splitter_2_has_been_resized = False

    step3_settings_ui = None
    step4_settings_ui = None
    log_id = None

    step3_config = {'is_automatic_masking': False,
                    'threshold_low': 0.05,
                    'threshold_high': 0.95,
                    'is_perform_savitzky_golay_filtering': False,
                    'window_size': 5,
                    'order_number': 1,
                    'is_interpolation_factor': False,
                    'interpolation_factor_value': 0,
                    'is_cross_section_mode': False,
                    }

    step4_config = {'is_cross_section_mode': False,
                    'is_perform_savitzky_golay_filtering': False,
                    'window_size': 5,
                    'order_number': 1,
                    'boundary_conditions_position': None,
                    'boundary_conditions_width': None,
                    'estimated_bragg_edge_position_value': np.NaN,
                    'estimated_bragg_edge_width_value': np.NaN,
                    'estimated_bragg_edge_height_value': np.NaN,
                    'estimated_bragg_edge_position_range': np.NaN,
                    'estimated_bragg_edge_width_range': np.NaN,
                    'is_automatic_masking': False,
                    'threshold_low': 0.05,
                    'threshold_high': 0.95,
                    'is_interpolation_factor': False,
                    'interpolation_factor_value': 0,
                    }

    def __init__(self, parent=None, o_roi=None, o_api=None):
        super(Interface, self).__init__(parent)

        if o_roi is None:
            self.debugging_mode = True

        self.o_roi = o_roi
        self.o_api = o_api

        if not self.debugging_mode:
            self.live_image = self.o_roi.live_image
            self.sample_projections = self.o_api.sample_projections
            self.untouched_sample_projections = copy.deepcopy(self.o_api.sample_projections)
            self.ob_projections = self.o_api.ob_projections
            self.list_roi = self.o_roi.list_roi

        self.is_mcp_correction_requested = o_api.is_mcp_corrected_ui.value

        ui_full_path = os.path.join(os.path.dirname(__file__), os.path.join('ui', 'ui_fit.ui'))
        self.ui = load_ui(ui_full_path, baseinstance=self)
        self.setWindowTitle("Fit Interface")

        o_init = Initialization(parent=self)
        o_init.widgets()
        if not self.debugging_mode:
            o_display = Display(parent=self)
            o_display.prepare_data_raw_image()

        # configuration of config
        o_get = Get(parent=self)
        log_file_name = o_get.log_file_name()
        logging.basicConfig(filename=log_file_name,
                            filemode='a',
                            format='[%(levelname)s] - %(asctime)s - %(message)s',
                            level=logging.INFO)
        logging.info("*** Starting a new session ***")
        # logging.info(f" Version: {versioneer.get_version()}")

    def log_button_clicked(self):
        LogLauncher(parent=self)

    def moving_average_size_radioButton_clicked(self):
        o_gui = GuiHandler(parent=self)
        o_gui.moving_average_size_radioButton_clicked()

    @wait_cursor
    def prepare_data_button_clicked(self):
        o_event = PrepareData(parent=self)
        o_event.prepare_data()

    def activate_moving_average_clicked(self):
        status = self.ui.activate_moving_average_checkBox.isChecked()
        self.ui.moving_average_groupBox.setEnabled(status)

    def import_prepared_data_clicked(self):
        o_import = Import(parent=self)
        o_import.run()

    def export_prepared_data_clicked(self):
        o_export = Export(parent=self)
        o_export.run()

    def pixel_marker_changed(self):
        region = self.pixel_marker_item.getArraySlice(self.live_image,
                                                      self.ui.image_view.imageItem)

        # x0 = region[0][0].start
        # y0 = region[0][1].start
        # x1 = region[0][0].stop
        # y1 = region[0][1].stop

        y0 = region[0][0].start
        x0 = region[0][1].start
        y1 = region[0][0].stop
        x1 = region[0][1].stop

        x = np.int(np.mean([x0, x1]))
        y = np.int(np.mean([y0, y1]))

        self.pixel_marker['x'] = x
        self.pixel_marker['y'] = y

        o_display = Display(parent=self)
        o_display.cross_of_pixel_to_fit()
        self.check_status_of_fit_buttons()
        o_display.profile()

    def check_status_of_fit_buttons(self):
        o_event = EventHandler(parent=self)
        o_event.check_status_of_fit_buttons()

    def calculate_profile_of_pixel_selected(self):
        pixel_marker = self.pixel_marker
        x = pixel_marker['x']
        y = pixel_marker['y']
        normalize_projections = self.normalize_projections

        logging.info("calculate profile of pixel selected")
        logging.info(f"-> x:{x}, y:{y}")
        logging.info(f"np.shape(normalize_projections): {np.shape(normalize_projections)}")

        profile = normalize_projections[:, y, x]
        self.profile_of_pixel_selected = profile

    def calculate_profile_of_roi(self):
        # normalize_projections_lambda_x_y = self.normalize_projections.transpose(2, 0, 1)
        normalize_projections_lambda_x_y = self.normalize_projections
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
#                total_counts_of_roi += np.sum(_projection[_y0: _y1 + 1, _x0: _x1 + 1])
                total_counts_of_roi += np.sum(_projection[_x0: _x1 + 1, _y0: _y1 + 1])

            mean_counts_of_roi = total_counts_of_roi / total_number_of_pixels_in_rois
            profile_y_axis.append(mean_counts_of_roi)
        self.profile_y_axis = profile_y_axis

    def number_of_files_to_exclude_slider_changed(self, value):
        left_number = np.int(str(self.ui.left_number_of_files_to_exclude_slider.value()))
        right_number = np.int(str(self.ui.right_number_of_files_to_exclude_slider.value()))
        self.nbr_files_to_exclude_from_plot['left'] = left_number
        self.nbr_files_to_exclude_from_plot['right'] = right_number

        o_display = Display(parent=self)
        o_display.profile()

    def rough_lambda_peak_position_slider_changed(self, value):
        o_display = Display(parent=self)
        o_display.rough_peak_position()

    def lambda_range_changed(self):
        o_get = Get(parent=self)
        rough_peak_position = o_get.rough_peak_position()

        new_lambda_range = self.bragg_peak_range_ui.getRegion()
        left_index = Get.nearest_index(self.o_roi.lambda_array, new_lambda_range[0])
        right_index = Get.nearest_index(self.o_roi.lambda_array, new_lambda_range[1])

        if (rough_peak_position <= new_lambda_range[0]) or (rough_peak_position >= new_lambda_range[1]):
            self.rough_peak_index_position = np.mean([left_index, right_index])

        rough_peak_index_position = Get.nearest_index(self.o_roi.lambda_array, rough_peak_position)

        self.ui.rough_lambda_peak_position_slider.setMinimum(left_index)
        self.ui.rough_lambda_peak_position_slider.setMaximum(right_index)
        self.ui.rough_lambda_peak_position_slider.setValue(rough_peak_index_position)

        self.number_of_files_to_exclude_slider_changed(0)

    def normalization_mode_checkBox_clicked(self):
        is_with_normalization = self.ui.normalization_checkBox.isChecked()
        self.ui.normalization_groupBox.setEnabled(is_with_normalization)

    def select_background_roi_button_clicked(self):
        pass
        # self.o_normalization_roi_gui = NormRoiSelection(parent=self,
        #                                                 main_api=self.o_api)
        # self.o_normalization_roi_gui.show()

    def update_number_of_normalization_roi(self, nbr_of_roi=0):
        p = inflect.engine()
        message = "{} ".format(nbr_of_roi) + p.plural("ROI", nbr_of_roi) + " selected!"
        self.ui.normalization_message_label.setText(message)

    def calculate_mask(self):
        o_event = EventHandler(parent=self)
        o_event.calculate_mask()

    def normalization_radioButton_clicked(self):
        state = self.ui.normal_normalization_radioButton.isChecked()
        self.ui.select_background_roi_pushButton.setEnabled(state)
        self.ui.normalization_message_label.setEnabled(state)

    @wait_cursor
    def fit_pixel_clicked(self):
        o_fit = FitHandler(parent=self)
        o_fit.fit(mode='pixel')

    @wait_cursor
    def fit_full_roi_clicked(self):
        o_fit = FitHandler(parent=self)
        o_fit.fit(mode='full')

    def kernel_dimension_changed(self):
        if self.ui.kernel_dimension_3d_radioButton.isChecked():
            third_dimension_state = True
            self.ui.kernel_size_default_label.setText(self.default_kernel_size_label['3d'])
        else:
            third_dimension_state = False
            self.ui.kernel_size_default_label.setText(self.default_kernel_size_label['2d'])
        self.ui.kernel_size_custom_lambda_label.setVisible(third_dimension_state)
        self.ui.kernel_size_custom_lambda_lineEdit.setVisible(third_dimension_state)

    def step3_settings_button_clicked(self):
        if self.step3_settings_ui:
            self.step3_settings_ui.setFocus()
            self.step3_settings_ui.activateWindow()
        else:
            step3_settings_ui = Step3SettingsHandler(parent=self)
            step3_settings_ui.show()
            self.step3_settings_ui = step3_settings_ui

    def step4_settings_button_clicked(self):
        if self.step4_settings_ui:
            self.step4_settings_ui.setFocus()
            self.step4_settings_ui.activateWindow()
        else:
            step4_settings_ui = Step4SettingsHandler(parent=self)
            step4_settings_ui.show()
            self.step4_settings_ui = step4_settings_ui

    # result tab
    def export_results_pushButton_clicked(self):
        o_export = Export(parent=self)
        o_export.result()

    def apply_clicked(self):
        self.close()

    def closeEvent(self, c):
        logging.info("*** Ending session ***")
