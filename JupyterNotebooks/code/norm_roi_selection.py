import os
from JupyterNotebooks.code import load_ui
import numpy as np

from JupyterNotebooks.code.roi_selection import Interface as RoiInterface


class NormRoiSelection(RoiInterface):

    def __init__(self, parent=None, main_api=None):
        super(RoiInterface, self).__init__(parent)
        self.parent = parent

        self.sample_projections_lambda_x_y = main_api.sample_projections  # lambda, y, x
        # self.sample_projections_lambda_x_y = self.sample_projections.transpose(2, 0, 1)  # lambda, x, y
        self.tof_array = main_api.tof_array
        self.lambda_array = main_api.lambda_array

        ui_full_path = os.path.join(os.path.dirname(__file__), os.path.join('ui', 'ui_roi_selection.ui'))
        self.ui = load_ui(ui_full_path, baseinstance=self)
        self.setWindowTitle("Select region to fit!")

        self.init_widgets()
        self.display_image()

    def apply_clicked(self):
        mask = np.zeros(np.shape(self.live_image))
        for _index_roi in self.list_roi.keys():
            _roi = self.list_roi[_index_roi]
            _x0 = _roi['x0']
            _y0 = _roi['y0']
            _x1 = _roi['x1']
            _y1 = _roi['y1']
            mask[_y0: _y1 + 1, _x0: _x1 + 1] = 1
        self.mask = mask

        self.parent.update_number_of_normalization_roi(nbr_of_roi=len(self.list_roi))
        self.parent.ob_list_roi = self.list_roi

        self.close()
