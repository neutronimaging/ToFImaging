import numpy as np
import glob
from collections import Counter
import os
from os.path import expanduser


class Get:

    def __init__(self, parent=None):
        self.parent = parent

    def kernel_type(self):
        if self.parent.ui.kernel_type_box_radioButton.isChecked():
            return 'box'
        elif self.parent.ui.kernel_type_gaussian_radioButton.isChecked():
            return 'gaussian'
        else:
            raise NotImplementedError("kernel type not implemented!")

    def kernel_size(self):
        if self.parent.ui.kernel_size_default_radioButton.isChecked():
            return self.parent.default_kernel_size
        else:
            y = np.int(str(self.parent.ui.kernel_size_custom_y_lineEdit.text()))
            x = np.int(str(self.parent.ui.kernel_size_custom_x_lineEdit.text()))
            l = np.int(str(self.parent.ui.kernel_size_custom_lambda_lineEdit.text()))
            return {'y': y, 'x': x, 'lambda': l}

    def kernel_dimensions(self):
        if self.parent.ui.kernel_dimension_3d_radioButton.isChecked():
            return '3d'
        else:
            return '2d'

    def normalization_mode(self):
        if self.parent.ui.normalization_by_roi_radioButton.isChecked():
            return 'by ROI'
        elif self.parent.ui.normalization_pixel_by_pixel_radioButton.isChecked():
            return 'pixel by pixel'
        else:
            raise NotImplementedError("normalization mode not implemented!")

    def log_file_name(self):
        log_file_name = self.parent.config['log_file_name']
        full_log_file_name = Get.full_home_file_name(log_file_name)
        return full_log_file_name

    def rough_peak_position(self):
        return self.parent.rough_peak_ui.pos()[0]

    @staticmethod
    def nearest_index(array, value):
        idx = (np.abs(np.array(array) - value)).argmin()
        return idx

    @staticmethod
    def list_of_most_dominant_extension_from_folder(folder='', files=None):
        """
        This will return the list of files from the most dominant file extension found in the folder
        as well as the most dominant extension used
        """

        if folder:
            list_of_input_files = glob.glob(os.path.join(folder, '*'))
        else:
            list_of_input_files = files

        list_of_input_files.sort()
        list_of_base_name = [os.path.basename(_file) for _file in list_of_input_files]

        # work with the largest common file extension from the folder selected

        counter_extension = Counter()
        for _file in list_of_base_name:
            [_base, _ext] = os.path.splitext(_file)
            counter_extension[_ext] += 1

        dominant_extension = ''
        dominant_number = 0
        for _key in counter_extension.keys():
            if counter_extension[_key] > dominant_number:
                dominant_extension = _key
                dominant_number = counter_extension[_key]

        list_of_input_files = glob.glob(os.path.join(folder, '*' + dominant_extension))
        list_of_input_files.sort()

        list_of_input_files = [os.path.abspath(_file) for _file in list_of_input_files]

        return [list_of_input_files, dominant_extension]

    @staticmethod
    def full_home_file_name(base_file_name):
        home_folder = expanduser("~")
        full_log_file_name = os.path.join(home_folder, base_file_name)
        return full_log_file_name
