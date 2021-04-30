from qtpy import QtGui, QtCore
from qtpy.QtWidgets import QApplication
import numpy as np
import logging
import copy

from src.tofimaging.EdgeFitting import (GaussianBraggEdgeFitting2D, AdvancedBraggEdgeFitting2D,
                                        AdvancedDirectBraggEdgeFitting2D)
from jupyter_notebooks.code.utilities.get import Get
from jupyter_notebooks.code.display import Display


class FitHandler:

    BRAGG_EDGE_FIT_POSITION_TOLERANCE = 0.1    #  if pos is 4.5, [4.4, 4.6]
    BRAGG_EDGE_FIT_WIDTH_TOLERANCE = 0.01

    def __init__(self, parent=None):
        self.parent = parent

    def fit(self, mode='pixel'):
        logging.info("Fitting started ...")
        logging.info(f"-> mode: {mode}")
        self.parent.ui.setEnabled(False)

        # get algorithm selected
        o_get = Get(parent=self.parent)
        algorithm_selected = o_get.algorithm_selected()
        logging.info(f"-> algorithm selected: {algorithm_selected}")

        # get lambda range and full array
        lambda_range = self.parent.bragg_peak_range_ui.getRegion()
        lambda_array = self.parent.o_roi.lambda_array
        logging.info(f"-> len(lambda_array): {len(lambda_array)}")
        logging.info(f"-> lambda_range: {lambda_range}")

        # get estimated bragg peak position
        o_get = Get(parent=self.parent)
        est_position = o_get.rough_peak_position()
        logging.info(f"-> estimated bragg edge position: {est_position}")

        mask = self.parent.mask
        mask = np.transpose(mask)
        logging.info(f"-> np.shape(mask): {np.shape(mask)}")
        normalize_projections = self.parent.normalize_projections
        logging.info(f"-> np.shape(normalize_projections): {np.shape(normalize_projections)}")

        self.parent.ui.statusbar.setStyleSheet("color: blue")
        self.parent.ui.statusbar.showMessage("Fitting {} using {} algorithm ... IN PROGRESS".format(
                mode, algorithm_selected))
        QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        QtGui.QGuiApplication.processEvents()

        if mode == 'pixel':
            result = self.fit_pixel(algorithm_selected,
                                    normalize_projections,
                                    lambda_array,
                                    lambda_range,
                                    mask,
                                    est_position,
                                    config=self.parent.step3_config)
            logging.debug(f"just after self.fit_pixel, result is {result}")

            if result:

                logging.info("-> Result of fitting (pixel_mode):")
                if algorithm_selected == 'gaussian':
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
                    logging.info(f"--> edge width: {edge_width}")
                    logging.info(f"--> edge height: {result['edge_height']}")
                    logging.info(f"--> edge position: {edge_position}")
                    logging.info(f"--> post_bc: {pos_bc}")
                    logging.info(f"--> wid_bc: {wid_bc}")

                elif algorithm_selected == 'advanced':
                    logging.info(f"-> result: {result}")
                    self.parent.step4_config['alpha'] = result['alpha']
                    self.parent.step4_config['sigma'] = result['sigma']
                    self.parent.step4_config['a1'] = result['a1']
                    self.parent.step4_config['a2'] = result['a2']
                    self.parent.step4_config['a5'] = result['a5']
                    self.parent.step4_config['a6'] = result['a6']
                    self.parent.ui.step4_fit_roi_pushButton.setEnabled(True)
                    self.parent.ui.step4_fit_roi_settings_pushButton.setEnabled(True)

        elif mode == 'full':

            if algorithm_selected == 'gaussian':

                result = self.fit_full_roi(algorithm_selected,
                                           normalize_projections,
                                           lambda_array,
                                           lambda_range,
                                           mask,
                                           config=self.parent.step4_config)

                self.parent.full_fit_result = result
                self.parent.ui.toolBox.setItemEnabled(2, True)

                logging.info(f"result of full mode:")
                logging.info(f"-> shape(edge_position): {np.shape(result['edge_position'])}")
                logging.info(f"-> shape(edge_height): {np.shape(result['edge_height'])}")
                logging.info(f"-> shape(edge_width): {np.shape(result['edge_width'])}")
                logging.info(f"-> shape(edge_slope): {np.shape(result['edge_slope'])}")
                logging.info(f"-> shape(median_image): {np.shape(result['median_image'])}")

                o_display = Display(parent=self.parent)
                o_display.result_full_mode(algorithm_selected=algorithm_selected,
                                           input_image=self.parent.live_process_data,
                                           edge_position=np.transpose(result['edge_position']),
                                           edge_height=np.transpose(result['edge_height']),
                                           edge_width=np.transpose(result['edge_width']),
                                           edge_slope=np.transpose(result['edge_slope']),
                                           image_median=np.transpose(result['median_image']))

            elif algorithm_selected == 'advanced':
                alpha = self.parent.step4_config['alpha']
                sigma = self.parent.step4_config['sigma']

                result = self.fit_full_roi(algorithm_selected,
                                           normalize_projections,
                                           lambda_array,
                                           lambda_range,
                                           mask,
                                           alpha=alpha,
                                           sigma=sigma,
                                           config=self.parent.step4_config)

                logging.info(f"result of full mode")
                logging.info(f"-> shape(edge_position): {np.shape(result['edge_position'])}")
                logging.info(f"-> shape(edge_height): {np.shape(result['edge_height'])}")
                logging.info(f"-> shape(edge_width): {np.shape(result['edge_width'])}")
                logging.info(f"-> shape(median_image): {np.shape(result['median_image'])}")

                self.parent.advanced_full_fit_result = result
                self.parent.ui.toolBox.setItemEnabled(2, True)

                o_display = Display(parent=self.parent)
                o_display.result_full_mode(algorithm_selected=algorithm_selected,
                                           input_image=self.parent.live_process_data,
                                           edge_position=np.transpose(result['edge_position']),
                                           edge_height=np.transpose(result['edge_height']),
                                           edge_width=np.transpose(result['edge_width']))

            else:
                raise NotImplementedError("algorithm selected has not been implemented yet")

        self.parent.ui.setEnabled(True)
        self.parent.ui.statusbar.showMessage("Fitting {} using {} algorithm ... DONE".format(
                mode, algorithm_selected),
                                             15000)
        self.parent.ui.statusbar.setStyleSheet("color: green")
        QApplication.restoreOverrideCursor()

    def fit_full_roi(self, algorithm_selected, T_mavg, lambda_array, lambda_range, mask, config,
                     alpha=None,
                     sigma=None):

        logging.info(f"fitting full roi")
        logging.debug("entering fit_full_roi")

        est_position = config['estimated_bragg_edge_position_value']
        est_width = config['estimated_bragg_edge_width_value']
        est_height = config['estimated_bragg_edge_height_value']
        pos_bc = config['estimated_bragg_edge_position_range']
        wid_bc = config['estimated_bragg_edge_width_range']

        logging.info(f"-> estimated bragg edge position value: {est_position}")
        logging.info(f"-> estimated bragg edge width value: {est_width}")
        logging.info(f"-> estimated bragg edge height value: {est_height}")
        logging.info(f"-> estimated bragg edge position range: {pos_bc}")
        logging.info(f"-> estimated bragg edge width range: {wid_bc}")

        T_mavg = copy.deepcopy(T_mavg.transpose(2, 1, 0))  # lambda, y, x -> y, x, lambda
        logging.info(f"-> np.shape(T_mavg): {np.shape(T_mavg)}")
        mask = np.transpose(mask)

        logging.info(f"-> algorithm selected {algorithm_selected}")
        if algorithm_selected == 'gaussian':

            logging.info(f"--> about to run GaussianBraggEdgeFitting2D")
            fit_result = GaussianBraggEdgeFitting2D(T_mavg,
                                                    lambda_array,
                                                    lambda_range,
                                                    mask=mask,
                                                    bool_log=False,
                                                    est_pos=est_position,
                                                    est_wid=est_width,
                                                    est_h=est_height,
                                                    wid_BC=wid_bc,
                                                    pos_BC=pos_bc,
                                                    )
            logging.info(f"--> done with GaussianBraggEdgeFitting2D")
            return fit_result

        elif algorithm_selected == 'advanced':

            logging.info(f"--> about to run AdvancedBraggEdgeFitting2D")
            fit_result = AdvancedBraggEdgeFitting2D(T_mavg,
                                                    lambda_array,
                                                    lambda_range,
                                                    mask=mask,
                                                    est_pos=est_position,
                                                    est_alpha=alpha,
                                                    est_sigma=sigma,
                                                    bool_smooth=True,
                                                    smooth_w=5,
                                                    smooth_n=1,
                                                    )
            logging.info(f"--> done with AdvancedBraggEdgeFitting2D")
            return fit_result

        else:
            raise NotImplementedError("algorithm selected has not been implemented yet")

    def fit_pixel(self, algorithm_selected, T_mavg, lambda_array, lambda_range, mask, est_position, config):

        logging.debug("entering fix_pixel method")

        # get pixel coordinates
        pixel_marker = self.parent.pixel_marker
        pixel = [pixel_marker['x'],
                 pixel_marker['y']]
        logging.info(f"-> pixel: {pixel}")

        auto_mask = config['is_automatic_masking']
        logging.info(f"-> auto_mask: {auto_mask}")
        mask_thresh = [np.float(config['threshold_low']),
                       np.float(config['threshold_high'])]
        logging.info(f"-> mask_thresh: {mask_thresh}")
        bool_smooth = config['is_perform_savitzky_golay_filtering']
        logging.info(f"-> bool_smooth: {bool_smooth}")
        smooth_w = np.int(config['window_size'])
        logging.info(f"-> smooth_w: {smooth_w}")
        smooth_n = np.int(config['order_number'])
        logging.info(f"-> smooth_n: {smooth_n}")
        if config['is_interpolation_factor']:
            interp_factor = np.int(config['interpolation_factor_value'])
        else:
            interp_factor = 0
        logging.info(f"-> interp_factor: {interp_factor}")
        bool_log = config['is_cross_section_mode']
        logging.info(f"-> bool_log: {bool_log}")

        T_mavg = copy.deepcopy(T_mavg.transpose(2, 1, 0))  # lambda, y, x -> y, x, lambda
        logging.info(f"-> np.shape(T_mavg): {np.shape(T_mavg)}")
        mask = np.transpose(mask)

        logging.info(f"-> algorithm selected: {algorithm_selected}")
        if algorithm_selected == 'gaussian':
            logging.info(f"--> about to run GaussianBraggEdgeFitting2D")
            fit_result = GaussianBraggEdgeFitting2D(T_mavg,
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
            logging.info(f"--> done with GaussianBraggEdgeFitting2D")
            return fit_result

        elif algorithm_selected == 'advanced':
            logging.info(f"--> about to run AdvancedBraggEdgeFitting2D")
            fit_result = AdvancedBraggEdgeFitting2D(T_mavg,
                                                    lambda_array,
                                                    lambda_range,
                                                    mask=mask,
                                                    est_pos=est_position,
                                                    smooth_w=smooth_w,
                                                    smooth_n=smooth_n,
                                                    debug_idx=pixel,
                                                    bool_print=True)
            logging.info(f"--> done with AdvancedBraggEdgeFitting2D")
            return fit_result

        elif algorithm_selected == 'advanced direct':
            logging.info(f"--> about to run AdvancedDirectBraggEdgeFitting2D")
            fit_result = AdvancedDirectBraggEdgeFitting2D(T_mavg,
                                                          lambda_array,
                                                          lambda_range,
                                                          mask=mask,
                                                          est_pos=est_position,
                                                          bool_smooth=True,
                                                          smooth_w=smooth_w,
                                                          smooth_n=smooth_n,
                                                          debug_idx=pixel,
                                                          bool_print=True)
            logging.info(f"--> done with AdvancedDirectBraggEdgeFitting2D")
            return fit_result

        else:
            raise NotImplementedError("algorithm selected has not been implemented yet")
