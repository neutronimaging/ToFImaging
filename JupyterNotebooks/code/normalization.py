from qtpy import QtGui
import copy
import numpy as np

from NeuNorm.normalization import Normalization as NormalizationLibrary


class Normalization:

    def __init__(self, parent=None):
        self.parent = parent

    def run(self, flag=None, mode=None):

        self.parent.ui.statusbar.showMessage("Normalization ... IN PROGRESS")
        QtGui.QGuiApplication.processEvents()

        if flag:
            sample_projections = self.parent.sample_projections
            ob_projections = self.parent.ob_projections

            working_sample_projections = sample_projections.transpose(2, 1, 0)  # lambda, y, x
            working_ob_projections = ob_projections.transpose(2, 1, 0)  # lambda, y, x

            # working_sample_projections = sample_projections
            # working_ob_projections = ob_projections

            # normalize_projections = list()

            if self.parent.ui.normalization_pixel_by_pixel_radioButton.isChecked():

                normalize_projections = \
                    Normalization.normalization_pixel_by_pixel(working_ob_projections,
                                                               working_sample_projections)

            elif self.parent.ui.normalization_by_roi_radioButton.isChecked():

                list_roi = self.parent.o_roi.list_roi
                normalize_projections = Normalization.normalization_by_roi(list_roi,
                                                                           working_ob_projections,
                                                                           working_sample_projections)

            elif self.parent.ui.normal_normalization_radioButton.isChecked():

                list_roi = self.parent.ob_list_roi
                normalize_projections = Normalization.normalization_by_roi_of_ob(list_roi,
                                                                                 working_ob_projections,
                                                                                 working_sample_projections)

            self.parent.normalize_projections = normalize_projections.transpose(2, 1, 0)  # x, y, lambda

        else:  # no normalization
            self.parent.normalize_projections = copy.deepcopy(self.parent.sample_projections)

        QtGui.QGuiApplication.processEvents()

    @staticmethod
    def normalization_pixel_by_pixel(working_ob_projections,
                                     working_sample_projections):

        _sample = working_sample_projections
        _ob = working_ob_projections

        o_norm = NormalizationLibrary()
        o_norm.data['sample']['data'] = _sample
        o_norm.data['ob']['data'] = _ob

        # create fake list of sample and ob
        list_filename = ['N/A' for _ in np.arange(len(_sample))]
        o_norm.data['sample']['file_name'] = list_filename
        o_norm.data['ob']['file_name'] = list_filename

        o_norm.normalization()

        return np.array(o_norm.data['normalized'])

    @staticmethod
    def normalization_by_roi(list_roi,
                             working_ob_projections,
                             working_sample_projections):

        normalize_projections = list()

        total_number_of_pixels_in_rois = 0
        for _index_roi in list_roi.keys():
            _roi = list_roi[_index_roi]
            _x0 = _roi['x0']
            _y0 = _roi['y0']
            _x1 = _roi['x1']
            _y1 = _roi['y1']
            total_number_of_pixels_in_rois += (_y1 - _y0 + 1) * (_x1 - _x0 + 1)
        for _sample, _ob in zip(working_sample_projections, working_ob_projections):

            mean_ob_value = 0
            for _index_roi in list_roi.keys():
                _roi = list_roi[_index_roi]
                _x0 = _roi['x0']
                _y0 = _roi['y0']
                _x1 = _roi['x1']
                _y1 = _roi['y1']

                mean_ob_value += np.sum(_ob[_y0: _y1 + 1, _x0: _x1 + 1])

            mean_ob = mean_ob_value / total_number_of_pixels_in_rois
            _normalize = _sample / mean_ob

            normalize_projections.append(_normalize)
        normalize_projections = np.array(normalize_projections)
        return normalize_projections

    @staticmethod
    def normalization_by_roi_of_ob(ob_list_roi,
                                   working_ob_projections,
                                   working_sample_projections):

        _sample = working_sample_projections
        _ob = working_ob_projections

        o_norm = NormalizationLibrary()
        o_norm.data['sample']['data'] = _sample
        o_norm.data['ob']['data'] = _ob

        # create fake list of sample and ob
        list_filename = ['N/A' for _ in np.arange(len(_sample))]
        o_norm.data['sample']['file_name'] = list_filename
        o_norm.data['ob']['file_name'] = list_filename

        o_norm.normalization(roi=ob_list_roi)
        return np.array(o_norm.data['normalized'])
