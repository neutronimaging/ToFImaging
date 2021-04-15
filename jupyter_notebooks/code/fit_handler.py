from qtpy import QtGui, QtCore
from qtpy.QtWidgets import QApplication
import numpy as np

from ToFImaging.edgefitting_2D import GaussianBraggEdgeFitting_2D


class FitHandler:

    BRAGG_EDGE_FIT_POSITION_TOLERANCE = 0.1    #  if pos is 4.5, [4.4, 4.6]
    BRAGG_EDGE_FIT_WIDTH_TOLERANCE = 0.01

    def __init__(self, parent=None):
        self.parent = parent

    def fit(self, mode='pixel'):
        self.parent.ui.setEnabled(False)

        # get algorithm selected
        algorithm_selected = self.get_algorithm_selected()

        # get lambda range and full array
        lambda_range = self.parent.bragg_peak_range_ui.getRegion()
        lambda_array = self.parent.o_roi.lambda_array

        # get estimated bragg peak position
        est_position = self.parent.get_rough_peak_position()

        mask = self.parent.mask
        T_mavg = self.parent.T_mavg

        self.parent.ui.statusbar.showMessage("Fitting {} using {} algorithm ... IN PROGRESS".format(
                mode, algorithm_selected))
        QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        QtGui.QGuiApplication.processEvents()

        if mode == 'pixel':
            result = self.fit_pixel(algorithm_selected,
                                    T_mavg,
                                    lambda_array,
                                    lambda_range,
                                    mask,
                                    est_position,
                                    config=self.parent.step3_config)
            if result:
                self.parent.step4_config['estimated_bragg_edge_width_value'] = result['edge_width']
                self.parent.step4_config['estimated_bragg_edge_height_value'] = result['edge_height']
                self.parent.step4_config['estimated_bragg_edge_position_value'] = result['edge_position']
                edge_position = result['edge_position']
                edge_width = result['edge_width']
                pos_bc = [edge_position - self.BRAGG_EDGE_FIT_POSITION_TOLERANCE,
                          edge_position + self.BRAGG_EDGE_FIT_POSITION_TOLERANCE]
                wid_bc = [edge_width - self.BRAGG_EDGE_FIT_WIDTH_TOLERANCE,
                          edge_width + self.BRAGG_EDGE_FIT_WIDTH_TOLERANCE]
                self.parent.step4_config['estimated_bragg_edge_position_range'] = pos_bc
                self.parent.step4_config['estimated_bragg_edge_width_range'] = wid_bc

                self.parent.ui.step4_fit_roi_pushButton.setEnabled(True)
                self.parent.ui.step4_fit_roi_settings_pushButton.setEnabled(True)

        if mode == 'full':
            est_width = self.parent.pixel_fit_result['edge_width']
            est_position = self.parent.pixel_fit_result['edge_position']
            result = self.fit_full_roi(algorithm_selected,
                                       T_mavg,
                                       lambda_array,
                                       lambda_range,
                                       mask,
                                       est_position,
                                       est_width)
            self.parent.full_fit_result = result

        self.parent.ui.setEnabled(True)
        self.parent.ui.statusbar.showMessage("Fitting {} using {} algorithm ... DONE".format(
                mode, algorithm_selected),
                                             15000)
        QApplication.restoreOverrideCursor()

    def fit_full_roi(self, algorithm_selected, T_mavg, lambda_array, lambda_range, mask, config):

        if algorithm_selected == 'gaussian':

            est_position = config['estimated_bragg_edge_position_value']
            est_width = config['estimated_bragg_edge_width_value']
            est_height = config['estimated_bragg_edge_height_value']
            pos_bc = config['estimated_bragg_edge_position_range']
            wid_bc = config['estimated_bragg_edge_width_range']

            fit_result = GaussianBraggEdgeFitting_2D(T_mavg,
                                                     lambda_array,
                                                     lambda_range,
                                                     mask=mask,
                                                     bool_log=False,
                                                     est_position=est_position,
                                                     est_wid=est_width,
                                                     est_h=est_height,
                                                     wid_BC=wid_bc,
                                                     pos_BC=pos_bc,
                                                     )
            return fit_result

        else:
            raise NotImplementedError("algorithm selected has not been implemented yet")

    def fit_pixel(self, algorithm_selected, T_mavg, lambda_array, lambda_range, mask, est_position, config):

        # get pixel coordinates
        pixel_marker = self.parent.pixel_marker
        pixel = [pixel_marker['x'],
                 pixel_marker['y']]

        auto_mask = config['is_automatic_masking']
        mask_thresh = [np.float(config['threshold_low']),
                       np.float(config['threshold_high'])]
        bool_smooth = config['is_perform_savitzky_golay_filtering']
        smooth_w = np.int(config['window_size'])
        smooth_n = np.int(config['order_number'])
        if config['is_interpolation_factor']:
            interp_factor = np.int(config['interpolation_factor_value'])
        else:
            interp_factor = 0
        bool_log = config['is_cross_section_mode']

        if algorithm_selected == 'gaussian':
            fit_result = GaussianBraggEdgeFitting_2D(T_mavg,
                                                     lambda_array,
                                                     lambda_range,
                                                     mask=mask,
                                                     est_pos=est_position,
                                                     pos_BC=lambda_range,
                                                     debug_idx=pixel,
                                                     bool_save=True,
                                                     bool_smooth=bool_smooth,
                                                     smooth_w=smooth_w,
                                                     smooth_n=smooth_n,
                                                     auto_mask=auto_mask,
                                                     mask_thresh=mask_thresh,
                                                     interp_factor=interp_factor,
                                                     bool_log=bool_log)
            return fit_result

        else:
            raise NotImplementedError("algorithm selected has not been implemented yet")

    def get_algorithm_selected(self):
        if self.parent.ui.gaussian_radioButton.isChecked():
            return 'gaussian'
        elif self.parent.ui.advanced_radioButton.isChecked():
            return 'advanced'

        raise NotImplementedError