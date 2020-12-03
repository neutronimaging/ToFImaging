from qtpy import QtGui, QtCore
from qtpy.QtWidgets import QApplication

from ToFImaging.edgefitting_2D import GaussianBraggEdgeFitting_2D


class FitHandler:

    def __init__(self, parent=None):
        self.parent = parent

    def fit_pixel(self):

        self.parent.ui.setEnabled(False)

        # get pixel coordinates
        pixel_marker = self.parent.pixel_marker
        pixel = [pixel_marker['x'],
                 pixel_marker['y']]

        # get algorithm selected
        algorithm_selected = self.get_algorithm_selected()

        # get lambda range and full array
        lambda_range = self.parent.bragg_peak_range_ui.getRegion()
        lambda_array = self.parent.o_roi.lambda_array

        # get estimated bragg peak position
        est_position = self.parent.get_rough_peak_position()

        mask = self.parent.o_api.mask
        T_mavg = self.parent.o_api.T_mavg

        self.parent.ui.statusbar.showMessage("Fitting pixel using {} algorithm ... IN PROGRESS".format(
                                             algorithm_selected))
        QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        QtGui.QGuiApplication.processEvents()

        if algorithm_selected == 'gaussian':
            GaussianBraggEdgeFitting_2D(T_mavg,
                                        lambda_array,
                                        lambda_range,
                                        mask=mask,
                                        est_pos=est_position,
                                        pos_BC=lambda_range,
                                        debug_idx=pixel)

        self.parent.ui.setEnabled(True)
        self.parent.ui.statusbar.showMessage("Fitting pixel using {} algorithm ... DONE".format(algorithm_selected),
                                             15000)
        QApplication.restoreOverrideCursor()

    def get_algorithm_selected(self):
        if self.parent.ui.gaussian_radioButton.isChecked():
            return 'gaussian'
        elif self.parent.ui.advanced_radioButton.isChecked():
            return 'advanced'

        raise NotImplementedError
