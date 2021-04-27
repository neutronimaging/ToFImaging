import numpy as np
import pyqtgraph as pg
import logging

from jupyter_notebooks.code.parent import Parent
from jupyter_notebooks.code.utilities.get import Get
from jupyter_notebooks.code.event_handler import EventHandler

MARKER_HEIGHT, MARKER_WIDTH = 20, 20
COLOR_LAMBDA_RANGE = [250, 128, 247]
COLOR_ROUGH_LAMBDA = [50, 255, 50]


class Display(Parent):

    def prepare_data_raw_image(self):
        # live_image = np.transpose(self.parent.live_image)
        live_image = self.parent.live_image
        self.parent.ui.raw_image_view.setImage(live_image)

    def fit_data_tab(self):
        self.clear_previous_items()
        self.initialize_pixel_marker()
        self.image()
        self.profile()
        self.roi()
        self.cross_of_pixel_to_fit()
        self.box_around_pixel_to_fit()
        self.lambda_range_to_fit()
        self.init_rough_peak_slider()
        self.rough_peak_position()
        o_event = EventHandler(parent=self.parent)
        o_event.check_status_of_fit_buttons()

    def clear_previous_items(self):
        if self.parent.bragg_peak_range_ui:
            self.parent.ui.plot_view.removeItem(self.parent.bragg_peak_range_ui)
            self.parent.bragg_peak_range_ui = None
        if self.parent.roi_id:
            self.parent.ui.image_view.removeItem(self.parent.roi_id)
            self.parent.roi_id = None
        if self.parent.cross_of_pixel_to_fit:
            self.parent.ui.image_view.removeItem(self.parent.cross_of_pixel_to_fit)
            self.parent.cross_of_pixel_to_fit = None
        if self.parent.pixel_marker_item:
            self.parent.ui.image_view.removeItem(self.parent.pixel_marker_item)
            self.parent.pixel_marker_item = None
        if self.parent.rough_peak_ui:
            self.parent.ui.plot_view.removeItem(self.parent.rough_peak_ui)
            self.parent.rough_peak_ui = None

    def init_rough_peak_slider(self):
        lambda_range = self.parent.bragg_peak_range_ui.getRegion()
        minimum_value = Get.nearest_index(self.parent.o_roi.lambda_array, lambda_range[0])
        maximum_value = Get.nearest_index(self.parent.o_roi.lambda_array, lambda_range[1])

        current_rough_peak_value = self.parent.ui.rough_lambda_peak_position_slider.value()
        if (current_rough_peak_value <= minimum_value) or (current_rough_peak_value >= maximum_value):
            current_rough_peak_value = np.int(np.mean([minimum_value, maximum_value]))
            self.parent.rough_peak_index_position = current_rough_peak_value

        self.parent.ui.rough_lambda_peak_position_slider.setMinimum(minimum_value)
        self.parent.ui.rough_lambda_peak_position_slider.setMaximum(maximum_value)
        self.parent.ui.rough_lambda_peak_position_slider.setValue(current_rough_peak_value)

    def initialize_pixel_marker(self):
        logging.debug("initialize_pixel_marker")
        for _index_roi in self.parent.o_roi.list_roi:
            logging.debug(f"-> _index_roi: {_index_roi}")
            _roi = self.parent.o_roi.list_roi[_index_roi]
            logging.debug(f"-> _roi: {_roi}")
            _x0 = _roi['x0']
            _y0 = _roi['y0']
            _x1 = _roi['x1']
            _y1 = _roi['y1']

            x = np.mean([_x0, _x1])
            y = np.mean([_y0, _y1])

            self.parent.pixel_marker = {'x': np.int(x),
                                        'y': np.int(y)}
            logging.debug(f"-> self.parent.pixel_marker: {self.parent.pixel_marker}")

    def prepare_data_preview_image(self):
        self.parent.ui.process_image_view.clear()
        prepare_data = self.parent.normalize_projections
        self.parent.live_process_data = np.mean(prepare_data, axis=0)
        # live_process_image = np.transpose(self.parent.live_process_data)
        live_process_image = self.parent.live_process_data
        self.parent.ui.process_image_view.setImage(live_process_image)

        self.parent.ui.process_image_view.view.getViewBox().setXLink('raw image')
        self.parent.ui.process_image_view.view.getViewBox().setYLink('raw image')

    def image(self):
        self.parent.live_image = self.parent.o_roi.live_image
        # live_image = np.transpose(self.parent.live_image)
        live_image = self.parent.live_image
        self.parent.ui.image_view.clear()
        self.parent.ui.image_view.setImage(live_image)

    def profile(self):
        self.parent.ui.plot_view.clear()

        self.parent.calculate_profile_of_roi()
        x_axis = self.parent.o_roi.lambda_array

        # remove left and right files
        nbr_left = self.parent.nbr_files_to_exclude_from_plot['left']
        nbr_right = len(x_axis) - self.parent.nbr_files_to_exclude_from_plot['right']

        x_axis = x_axis[nbr_left: nbr_right]
        profile_y_axis = self.parent.profile_y_axis[nbr_left: nbr_right]
        self.parent.ui.plot_view.plot(x_axis, profile_y_axis)

        self.parent.calculate_profile_of_pixel_selected()

        profile_of_pixel_selected = self.parent.profile_of_pixel_selected[nbr_left: nbr_right]
        self.parent.ui.plot_view.plot(x_axis, profile_of_pixel_selected, pen=COLOR_LAMBDA_RANGE)

        if self.parent.bragg_peak_range_ui:
            self.parent.ui.plot_view.addItem(self.parent.bragg_peak_range_ui)

        self.rough_peak_position()

    def roi(self):
        for _index_roi in self.parent.o_roi.list_roi:
            _roi = self.parent.o_roi.list_roi[_index_roi]
            _id = _roi['id']

            pos = _id.pos()
            size = _id.size()

            new_id = pg.ROI(pos, size, movable=False)
            self.parent.roi_id = new_id
            self.parent.ui.image_view.addItem(new_id)

    def cross_of_pixel_to_fit(self):

        if self.parent.cross_of_pixel_to_fit:
            self.parent.ui.image_view.removeItem(self.parent.cross_of_pixel_to_fit)

        # x, y = self.parent.pixel_marker['x'], self.parent.pixel_marker['y']
        y, x = self.parent.pixel_marker['x'], self.parent.pixel_marker['y']

        pos = []
        adj = []

        # vertical guide
        pos.append([x, y - MARKER_HEIGHT / 2])
        pos.append([x, y + MARKER_HEIGHT / 2])
        adj.append([0, 1])

        # horizontal guide
        pos.append([x - MARKER_WIDTH / 2, y])
        pos.append([x + MARKER_WIDTH / 2, y])
        adj.append([2, 3])

        pos = np.array(pos)
        adj = np.array(adj)

        line_color = (255, 0, 0, 255, 1)
        lines = np.array([line_color for _ in np.arange(len(pos))],
                         dtype=[('red', np.ubyte), ('green', np.ubyte),
                                ('blue', np.ubyte), ('alpha', np.ubyte),
                                ('width', float)])
        self.parent.cross_of_pixel_to_fit = pg.GraphItem()
        self.parent.ui.image_view.addItem(self.parent.cross_of_pixel_to_fit)
        self.parent.cross_of_pixel_to_fit.setData(pos=pos,
                                                  adj=adj,
                                                  pen=lines,
                                                  symbol=None,
                                                  pxMode=False)

    def box_around_pixel_to_fit(self):
        # x, y = self.parent.pixel_marker['x'] - MARKER_WIDTH/2, self.parent.pixel_marker['y'] - MARKER_HEIGHT/2
        y, x = self.parent.pixel_marker['x'] - MARKER_WIDTH/2, self.parent.pixel_marker['y'] - MARKER_HEIGHT/2

        self.parent.pixel_marker_item = pg.ROI([x, y],
                                               [MARKER_WIDTH, MARKER_HEIGHT],
                                               scaleSnap=True)
        self.parent.ui.image_view.addItem(self.parent.pixel_marker_item)
        self.parent.pixel_marker_item.sigRegionChanged.connect(self.parent.pixel_marker_changed)

    def lambda_range_to_fit(self):
        if self.parent.default_bragg_peak_range is None:
            lambda_array = self.parent.o_roi.lambda_array
            nbr_lambda = len(lambda_array)
            self.parent.default_bragg_peak_range = [lambda_array[np.int(nbr_lambda / 2)],
                                                    lambda_array[np.int(nbr_lambda / 2) + 10]]

        self.parent.bragg_peak_range_ui = pg.LinearRegionItem(values=self.parent.default_bragg_peak_range,
                                                              orientation='vertical',
                                                              movable=True)
        self.parent.bragg_peak_range_ui.sigRegionChanged.connect(self.parent.lambda_range_changed)

        self.parent.bragg_peak_range_ui.setZValue(-10)
        self.parent.ui.plot_view.addItem(self.parent.bragg_peak_range_ui)
        
    def rough_peak_position(self):
        if self.parent.rough_peak_ui:
            self.parent.ui.plot_view.removeItem(self.parent.rough_peak_ui)

        rough_peak = self.parent.o_roi.lambda_array[self.parent.ui.rough_lambda_peak_position_slider.value()]
        self.parent.rough_peak_ui = pg.InfiniteLine(rough_peak,
                                                    pen=COLOR_ROUGH_LAMBDA,
                                                    label='Estimated Bragg Peak')
        self.parent.ui.plot_view.addItem(self.parent.rough_peak_ui)
