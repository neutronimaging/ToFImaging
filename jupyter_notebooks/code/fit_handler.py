from qtpy import QtGui, QtCore
from qtpy.QtWidgets import QApplication

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

        mask = self.parent.o_api.mask
        T_mavg = self.parent.o_api.T_mavg

        self.parent.ui.statusbar.showMessage("Fitting {} using {} algorithm ... IN PROGRESS".format(
                mode, algorithm_selected))
        QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        QtGui.QGuiApplication.processEvents()

        if mode == 'pixel':
            result = self.fit_pixel(algorithm_selected, T_mavg, lambda_array, lambda_range, mask,
                                    est_position)
            self.parent.pixel_fit_result = {'edge_width': result['edge_width'],
                                      'edge_height': result['edge_height'],
                                      'edge_position': result['edge_position']}

            if result:
                self.parent.ui.step4_fit_roi_pushButton.setEnabled(True)

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

    def fit_full_roi(self, algorithm_selected, T_mavg, lambda_array, lambda_range, mask, est_position, est_width):

        if algorithm_selected == 'gaussian':

            pos_bc = [est_position - self.BRAGG_EDGE_FIT_POSITION_TOLERANCE,
                      est_position + self.BRAGG_EDGE_FIT_POSITION_TOLERANCE]
            wid_bc = [est_width - self.BRAGG_EDGE_FIT_WIDTH_TOLERANCE,
                      est_width + self.BRAGG_EDGE_FIT_WIDTH_TOLERANCE]

            fit_result = GaussianBraggEdgeFitting_2D(T_mavg,
                                                     lambda_array,
                                                     lambda_range,
                                                     mask=mask,
                                                     bool_log=False,
                                                     )
            return fit_result

        else:
            raise NotImplementedError("algorithm selected has not been implemented yet")

    def fit_pixel(self, algorithm_selected, T_mavg, lambda_array, lambda_range, mask, est_position):

        # get pixel coordinates
        pixel_marker = self.parent.pixel_marker
        pixel = [pixel_marker['x'],
                 pixel_marker['y']]

        if algorithm_selected == 'gaussian':
            fit_result = GaussianBraggEdgeFitting_2D(T_mavg,
                                                     lambda_array,
                                                     lambda_range,
                                                     mask=mask,
                                                     est_pos=est_position,
                                                     pos_BC=lambda_range,
                                                     debug_idx=pixel,
                                                     bool_save=True,
                                                     bool_log=True,
                                                     bool_smooth=True,
                                                     smooth_w=5,
                                                     smooth_n=1)
            return fit_result

        else:
            raise NotImplementedError("algorithm selected has not been implemented yet")

    def get_algorithm_selected(self):
        if self.parent.ui.gaussian_radioButton.isChecked():
            return 'gaussian'
        elif self.parent.ui.advanced_radioButton.isChecked():
            return 'advanced'

        raise NotImplementedError
