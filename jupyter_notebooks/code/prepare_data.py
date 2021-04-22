from qtpy import QtGui
import numpy as np
import logging

from src.tofimaging.ReductionTools import KernelType, data_filering

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
        logging.info("Preparing data ...")
        self.parent.ui.setEnabled(False)

        # collect parameters
        is_with_normalization = self.parent.o_api.normalization_flag_ui.value
        logging.info(f"> with normalization:  {is_with_normalization}")
        self.parent.is_with_normalization = is_with_normalization

        o_get = Get(parent=self.parent)
        normalization_mode = o_get.normalization_mode()
        kernel_dimension = o_get.kernel_dimensions()
        kernel_size = o_get.kernel_size()
        kernel_type = o_get.kernel_type()
        logging.info(f"> normalization mode: {normalization_mode}")
        logging.info(f"> kernel size: {kernel_size}")
        logging.info(f"> kernel type: {kernel_type}")

        is_with_moving_average = self.parent.ui.activate_moving_average_checkBox.isChecked()
        logging.info(f"> moving average: {is_with_moving_average}")
        if is_with_moving_average:

            if not self.can_we_use_buffered_data(kernel_dimension=kernel_dimension,
                                                 kernel_size=kernel_size,
                                                 kernel_type=kernel_type):
                self.calculate_moving_average(kernel_dimension=kernel_dimension,
                                              kernel_size=kernel_size,
                                              kernel_type=kernel_type)
            else:
                logging.info(f"-> using buffered data")

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

    def can_we_use_buffered_data(self, kernel_dimension=None, kernel_size=None, kernel_type=None):
        if self.parent.moving_average_config is None:
            return False

        buffered_kernel_dimension = self.parent.moving_average_config.get('kernel_dimension', None)
        if buffered_kernel_dimension is None:
            return False
        if not (kernel_dimension == buffered_kernel_dimension):
            return False

        buffered_kernel_size = self.parent.moving_average_config.get('kernel_size', None)
        if buffered_kernel_size is None:
            return False
        if not (buffered_kernel_size['x'] == kernel_size['x']):
            return False
        if not (buffered_kernel_size['y'] == kernel_size['y']):
            return False
        if kernel_dimension == '3d':
            if not (buffered_kernel_size['lambda'] == kernel_size['lambda']):
                return False

        buffered_kernel_type = self.parent.moving_average_config.get('kernel_type', None)
        if buffered_kernel_type is None:
            return False
        if not (buffered_kernel_type == kernel_type):
            return False

        return True

    def calculate_moving_average(self, kernel_dimension=None, kernel_size=None, kernel_type=None):

        logging.info("-> calculate moving average ... in progress")
        self.parent.ui.statusbar.showMessage("Moving Average of Sample ... IN PROGRESS")
        QtGui.QGuiApplication.processEvents()

        x = kernel_size['x']
        y = kernel_size['y']
        kernel = [y, x]
        if kernel_dimension == '3d':
            l = kernel_size['lambda']
            kernel.append(l)

        logging.debug(f"--> kernel dimension: {kernel_dimension}")
        logging.debug(f"--> kernel shape: {np.shape(kernel)}")
        logging.debug(f"--> len(sample_projections): {len(self.parent.sample_projections)}")
        logging.debug(f"--> kernel: {kernel}")

        # self.parent.sample_projections = reduction_tools.moving_average_2D(self.parent.sample_projections,
        #                                                                    custom_kernel=kernel)
        self.parent.sample_projections = ReductionTools.data_filering(self.parent.sample_projections,
                                                                      kernel=kernel,
                                                                      kernel_type=kernel_type)

        if self.parent.is_with_normalization:
            self.parent.ui.statusbar.showMessage("Moving Average of OB ... IN PROGRESS")
            QtGui.QGuiApplication.processEvents()

            # self.parent.ob_projections = reduction_tools.moving_average_2D(self.parent.ob_projections,
            #                                                                custom_kernel=kernel)
            self.parent.ob_projections = ReductionTools.data_filering(self.parent.ob_projections,
                                                                      kernel=kernel,
                                                                      kernel_type=kernel_type)

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
        logging.info("-> calculate moving average ... done!")

