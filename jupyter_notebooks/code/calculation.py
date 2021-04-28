import numpy as np

from jupyter_notebooks.code.parent import Parent


class Calculation(Parent):

    def profile_of_roi(self):
        # normalize_projections_lambda_x_y = self.normalize_projections.transpose(2, 0, 1)
        normalize_projections_lambda_x_y = self.parent.normalize_projections
        list_roi = self.parent.o_roi.list_roi
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
        self.parent.profile_y_axis = profile_y_axis

  def profile_of_pixel_selected(self):
        pixel_marker = self.parent.pixel_marker
        x = pixel_marker['x']
        y = pixel_marker['y']
        normalize_projections = self.parent.normalize_projections
        profile = normalize_projections[:, y, x]
        self.parent.profile_of_pixel_selected = profile
