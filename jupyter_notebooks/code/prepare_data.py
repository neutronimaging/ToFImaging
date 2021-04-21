from qtpy import QtGui
import numpy as np

from ToFImaging import reduction_tools

from jupyter_notebooks.code.normalization import Normalization as LocalNormalization
from jupyter_notebooks.code.utilities.get import Get
from jupyter_notebooks.code.parent import Parent
from jupyter_notebooks.code.display import Display
from src.tofimaging import ReductionTools


class PrepareData(Parent):

    def calculate_mask(self):
        self.parent.ui.statusbar.showMessage("Calculate mask ... IN PROGRESS")
        QtGui.QGuiApplication.processEvents()

        list_roi = self.parent.list_roi
        [_, height, width] = np.shape(self.parent.normalize_projections)
        mask = np.zeros((height, width))

        for _roi_key in list_roi.keys():
            _roi = list_roi[_roi_key]
            x0 = _roi['x0']
            y0 = _roi['y0']
            x1 = _roi['x1']
            y1 = _roi['y1']
            mask[y0:y1 + 1, x0:x1 + 1] = 1
        self.parent.mask = mask

        QtGui.QGuiApplication.processEvents()

    def prepare_data(self):

        self.parent.ui.setEnabled(False)

        # collect parameters
        is_with_normalization = self.parent.o_api.normalization_flag_ui.value
        self.parent.is_with_normalization = is_with_normalization

        o_get = Get(parent=self.parent)
        normalization_mode = o_get.normalization_mode()
        kernel_dimension = o_get.kernel_dimensions()
        kernel_size = o_get.kernel_size()
        kernel_type = o_get.kernel_type()

        if self.parent.ui.activate_moving_average_checkBox.isChecked():

            if not self.parent.can_we_use_buffered_data(kernel_dimension=kernel_dimension,
                                                        kernel_size=kernel_size,
                                                        kernel_type=kernel_type):
                self.calculate_moving_average(kernel_dimension=kernel_dimension,
                                              kernel_size=kernel_size,
                                              kernel_type=kernel_type)

        o_norm = LocalNormalization(parent=self.parent)
        o_norm.run(flag=is_with_normalization,
                   mode=normalization_mode)

        self.calculate_mask()

        o_display = Display(parent=self.parent)
        o_display.prepare_data_preview_image()
        o_display.fit_data_tab()

        self.parent.ui.statusbar.showMessage("Prepare data ... Done!", 5000)
        QtGui.QGuiApplication.processEvents()
        self.parent.ui.toolBox.setItemEnabled(1, True)
        # self.parent.ui.toolBox.setCurrentIndex(1)  # switch to next tab
        self.parent.ui.export_prepared_data_pushButton.setEnabled(True)
        self.parent.ui.setEnabled(True)

    def calculate_moving_average(self, kernel_dimension=None, kernel_size=None, kernel_type=None):

        self.parent.ui.statusbar.showMessage("Moving Average of Sample ... IN PROGRESS")
        QtGui.QGuiApplication.processEvents()

        x = kernel_size['x']
        y = kernel_size['y']

        if kernel_dimension == '2d':
            kernel = np.ones((y, x))

        elif kernel_dimension == '3d':
            l = kernel_size['lambda']
            kernel = np.ones((y, x, l))

        else:
            raise NotImplementedError("kernel dimension does not exist!")

        # self.parent.sample_projections = reduction_tools.moving_average_2D(self.parent.sample_projections,
        #                                                                    custom_kernel=kernel)
        if kernel_type == 'box':
            ReductionTools.DataFiltering(self.parent.sample_projections,
                                         BoxKernel=kernel)
        elif kernel_type == 'gaussian':
            ReductionTools.DataFiltering(self.parent.sample_projections,
                                         GaussianKernel=kernel)
        else:
            raise ValueError(f"kernel type {kernel_type} not supported!")

        if self.parent.is_with_normalization:
            self.parent.ui.statusbar.showMessage("Moving Average of OB ... IN PROGRESS")
            QtGui.QGuiApplication.processEvents()

            # self.parent.ob_projections = reduction_tools.moving_average_2D(self.parent.ob_projections,
            #                                                                custom_kernel=kernel)
            if kernel_type == 'box':
                ReductionTools.DataFiltering(self.parent.ob_projections,
                                             BoxKernel=kernel)
            elif kernel_type == 'gaussian':
                ReductionTools.DataFiltering(self.parent.ob_projections,
                                             GaussianKernel=kernel)
            else:
                raise ValueError(f"kernel type {kernel_type} not supported!")

        self.parent.ui.statusbar.showMessage("Moving Average ... DONE!")
        QtGui.QGuiApplication.processEvents()

        # T_mavg = reduction_tools.moving_average_2D(self.normalize_projections,
        #                                            custom_kernel=kernel)
        # self.T_mavg = T_mavg

        moving_average_config = {'kernel_type'     : kernel_type,
                                 'kernel_size'     : kernel_size,
                                 'kernel_dimension': kernel_dimension}
        self.parent.moving_average_config = moving_average_config

        QtGui.QGuiApplication.processEvents()
