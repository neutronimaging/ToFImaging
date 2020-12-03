import os
from qtpy.QtWidgets import QMainWindow, QVBoxLayout
from jupyter_notebooks.code import load_ui
import pyqtgraph as pg
import numpy as np

from jupyter_notebooks.code.fit_handler import FitHandler
from jupyter_notebooks.code.utilities import find_nearest_index

MARKER_HEIGHT, MARKER_WIDTH = 20, 20
COLOR_LAMBDA_RANGE = [250, 128, 247]
COLOR_ROUGH_LAMBDA = [50, 255, 50]


class Interface(QMainWindow):

    o_roi = None
    o_api = None
    live_image = None

    pixel_marker = {'x': 0,
                    'y': 0}
    cross_of_pixel_to_fit = None

    default_bragg_peak_range = None
    bragg_peak_range_ui = None

    nbr_files_to_exclude_from_plot = {'left': 40,
                                      'right': 40}
    rough_peak_index_position = 0
    rough_peak_ui = None

    pixel_fit_result = None
    full_fit_result = None

    # list_roi = OrderedDict()
    # default_roi = {'x0': 0, 'y0': 0, 'x1': 50, 'y1': 50, 'id': None}
    # live_image = None
    # roi_width = 0.01
    # integrated_image_size = {'width': -1, 'height': -1}
    # splitter_2_has_been_resized = False

    def __init__(self, parent=None, o_roi=None, o_api=None):
        super(Interface, self).__init__(parent)

        self.o_roi = o_roi
        self.o_api = o_api

        # self.sample_projections = main_api.sample_projections  # x, y, lambda
        # self.sample_projections_lambda_x_y = self.sample_projections.transpose(2, 0, 1)  # lambda, x, y
        # self.tof_array = main_api.tof_array
        # self.lambda_array = main_api.lambda_array

        ui_full_path = os.path.join(os.path.dirname(__file__), os.path.join('ui', 'ui_fit.ui'))
        self.ui = load_ui(ui_full_path, baseinstance=self)
        self.setWindowTitle("Fit Interface")

        self.init_widgets()
        self.initialize_pixel_marker()
        self.display_image()
        self.display_profile()
        self.display_roi()
        self.display_cross_of_pixel_to_fit()
        self.display_box_around_pixel_to_fit()
        self.display_lambda_range_to_fit()
        self.init_rough_peak_slider()
        self.display_rough_peak_position()

    def init_rough_peak_slider(self):
        lambda_range = self.bragg_peak_range_ui.getRegion()
        minimum_value = find_nearest_index(self.o_roi.lambda_array, lambda_range[0])
        maximum_value = find_nearest_index(self.o_roi.lambda_array, lambda_range[1])

        current_rough_peak_value = self.ui.rough_lambda_peak_position_slider.value()
        if (current_rough_peak_value <= minimum_value) or (current_rough_peak_value >= maximum_value):
            current_rough_peak_value = np.int(np.mean([minimum_value, maximum_value]))
            self.rough_peak_index_position = current_rough_peak_value

        self.ui.rough_lambda_peak_position_slider.setMinimum(minimum_value)
        self.ui.rough_lambda_peak_position_slider.setMaximum(maximum_value)
        self.ui.rough_lambda_peak_position_slider.setValue(current_rough_peak_value)

    def initialize_pixel_marker(self):
        for _index_roi in self.o_roi.list_roi:
            _roi = self.o_roi.list_roi[_index_roi]
            _x0 = _roi['x0']
            _y0 = _roi['y0']
            _x1 = _roi['x1']
            _y1 = _roi['y1']

            x = np.mean([_x0, _x1])
            y = np.mean([_y0, _y1])

            self.pixel_marker = {'x': np.int(x),
                                 'y': np.int(y)}

    def init_widgets(self):
        # pyqtgraphs
        self.ui.image_view = pg.ImageView(view=pg.PlotItem())
        self.ui.image_view.ui.roiBtn.hide()
        self.ui.image_view.ui.menuBtn.hide()
        verti_layout1 = QVBoxLayout()
        verti_layout1.addWidget(self.ui.image_view)
        self.ui.widget_image.setLayout(verti_layout1)

        self.ui.plot_view = pg.PlotWidget()
        verti_layout2 = QVBoxLayout()
        verti_layout2.addWidget(self.ui.plot_view)
        self.ui.widget_plot.setLayout(verti_layout2)

        # splitters
        self.ui.splitter.setSizes([550, 50])
        self.ui.splitter_2.setSizes([200, 2])

        # sliders
        half_number_of_files = np.int(len(self.o_roi.lambda_array)/2)
        self.ui.right_number_of_files_to_exclude_slider.setMaximum(half_number_of_files)
        self.ui.left_number_of_files_to_exclude_slider.setMaximum(half_number_of_files)
        self.ui.right_number_of_files_to_exclude_slider.setValue(self.nbr_files_to_exclude_from_plot['left'])
        self.ui.left_number_of_files_to_exclude_slider.setValue(self.nbr_files_to_exclude_from_plot['right'])

    def pixel_marker_changed(self):
        region = self.pixel_marker_item.getArraySlice(self.live_image,
                                                      self.ui.image_view.imageItem)

        x0 = region[0][0].start
        y0 = region[0][1].start
        self.pixel_marker['x'] = x0
        self.pixel_marker['y'] = y0

        self.display_cross_of_pixel_to_fit()
        self.check_status_of_fit_buttons()
        self.display_profile()

    def check_status_of_fit_buttons(self):

        # we need to make sure the pixel selected is inside one of the ROI
        x_pixel, y_pixel = self.pixel_marker['x'], self.pixel_marker['y']

        list_roi = self.o_roi.list_roi
        for _index_roi in list_roi.keys():
            _roi = list_roi[_index_roi]
            _x0 = _roi['x0']
            _y0 = _roi['y0']
            _x1 = _roi['x1']
            _y1 = _roi['y1']

            if (x_pixel >= _x0) and (x_pixel <= _x1) and \
                (y_pixel >= _y0) and (y_pixel <= _y1):
                self.ui.step3_fit_pixel_pushButton.setEnabled(True)
                self.ui.statusbar.showMessage("")
                return
        self.ui.step3_fit_pixel_pushButton.setEnabled(False)
        self.ui.statusbar.showMessage("Pixel must be inside one of the ROI selected!")
        self.ui.statusbar.setStyleSheet("color: red")

    def display_lambda_range_to_fit(self):
        if self.default_bragg_peak_range is None:
            lambda_array = self.o_roi.lambda_array
            nbr_lambda = len(lambda_array)
            self.default_bragg_peak_range = [lambda_array[np.int(nbr_lambda/2)],
                                             lambda_array[np.int(nbr_lambda/2) + 10]]

        self.bragg_peak_range_ui = pg.LinearRegionItem(values=self.default_bragg_peak_range,
                                                       orientation=None,
                                                       movable=True)
        self.bragg_peak_range_ui.sigRegionChanged.connect(self.lambda_range_changed)

        self.bragg_peak_range_ui.setZValue(-10)
        self.ui.plot_view.addItem(self.bragg_peak_range_ui)

    def display_box_around_pixel_to_fit(self):
        x, y = self.pixel_marker['x'], self.pixel_marker['y']

        self.pixel_marker_item = pg.ROI([x, y],
                                   [MARKER_WIDTH, MARKER_HEIGHT],
                                   scaleSnap=True)
        self.ui.image_view.addItem(self.pixel_marker_item)
        self.pixel_marker_item.sigRegionChanged.connect(self.pixel_marker_changed)

    def display_cross_of_pixel_to_fit(self):

        if self.cross_of_pixel_to_fit:
            self.ui.image_view.removeItem(self.cross_of_pixel_to_fit)

        x, y = self.pixel_marker['x'], self.pixel_marker['y']

        pos = []
        adj = []

        # vertical guide
        pos.append([x + MARKER_WIDTH / 2, y - MARKER_HEIGHT / 2])
        pos.append([x + MARKER_WIDTH / 2, y + MARKER_HEIGHT + MARKER_HEIGHT / 2])
        adj.append([0, 1])

        # horizontal guide
        pos.append([x - MARKER_WIDTH / 2, y + MARKER_HEIGHT / 2])
        pos.append([x + MARKER_WIDTH + MARKER_WIDTH / 2, y + MARKER_HEIGHT / 2])
        adj.append([2, 3])

        pos = np.array(pos)
        adj = np.array(adj)

        line_color = (255, 0, 0, 255, 1)
        lines = np.array([line_color for _ in np.arange(len(pos))],
                         dtype=[('red', np.ubyte), ('green', np.ubyte),
                                ('blue', np.ubyte), ('alpha', np.ubyte),
                                ('width', float)])
        self.cross_of_pixel_to_fit = pg.GraphItem()
        self.ui.image_view.addItem(self.cross_of_pixel_to_fit)
        self.cross_of_pixel_to_fit.setData(pos=pos,
                                  adj=adj,
                                  pen=lines,
                                  symbol=None,
                                  pxMode=False)

    def display_roi(self):
        for _index_roi in self.o_roi.list_roi:
            _roi = self.o_roi.list_roi[_index_roi]
            _id = _roi['id']

            pos = _id.pos()
            size = _id.size()

            new_id = pg.ROI(pos, size, movable=False)
            self.ui.image_view.addItem(new_id)

    def display_image(self):
        self.live_image = self.o_roi.live_image
        live_image = np.transpose(self.live_image)
        self.ui.image_view.setImage(live_image)

    def display_profile(self):
        self.ui.plot_view.clear()

        self.calculate_profile_of_roi()
        x_axis = self.o_roi.lambda_array

        # remove left and right files
        nbr_left = self.nbr_files_to_exclude_from_plot['left']
        nbr_right = len(x_axis) - self.nbr_files_to_exclude_from_plot['right']

        x_axis = x_axis[nbr_left: nbr_right]
        profile_y_axis = self.profile_y_axis[nbr_left: nbr_right]
        self.ui.plot_view.plot(x_axis, profile_y_axis)

        self.calculate_profile_of_pixel_selected()

        profile_of_pixel_selected = self.profile_of_pixel_selected[nbr_left: nbr_right]
        self.ui.plot_view.plot(x_axis, profile_of_pixel_selected, pen=COLOR_LAMBDA_RANGE)

        if self.bragg_peak_range_ui:
            self.ui.plot_view.addItem(self.bragg_peak_range_ui)

        self.display_rough_peak_position()

    def display_rough_peak_position(self):
        if self.rough_peak_ui:
            self.ui.plot_view.removeItem(self.rough_peak_ui)

        rough_peak = self.o_roi.lambda_array[self.ui.rough_lambda_peak_position_slider.value()]
        self.rough_peak_ui = pg.InfiniteLine(rough_peak,
                                             pen=COLOR_ROUGH_LAMBDA,
                                             label='Estimated Bragg Peak')
        self.ui.plot_view.addItem(self.rough_peak_ui)

    def calculate_profile_of_pixel_selected(self):
        pixel_marker = self.pixel_marker
        x = pixel_marker['x']
        y = pixel_marker['y']
        normalize_projections = self.o_api.T_mavg

        profile = normalize_projections[y, x, :]
        self.profile_of_pixel_selected = profile

    def calculate_profile_of_roi(self):
        normalize_projections_lambda_x_y = self.o_api.normalize_projections.transpose(2, 0, 1)
        list_roi = self.o_roi.list_roi
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
                total_counts_of_roi += np.sum(_projection[_y0: _y1 + 1, _x0: _x1 + 1])

            mean_counts_of_roi = total_counts_of_roi / total_number_of_pixels_in_rois
            profile_y_axis.append(mean_counts_of_roi)
        self.profile_y_axis = profile_y_axis

    def number_of_files_to_exclude_slider_changed(self, value):
        left_number = np.int(str(self.ui.left_number_of_files_to_exclude_slider.value()))
        right_number = np.int(str(self.ui.right_number_of_files_to_exclude_slider.value()))
        self.nbr_files_to_exclude_from_plot['left'] = left_number
        self.nbr_files_to_exclude_from_plot['right'] = right_number
        self.display_profile()

    def rough_lambda_peak_position_slider_changed(self, value):
        self.display_rough_peak_position()

    def get_rough_peak_position(self):
        return self.rough_peak_ui.pos()[0]

    def lambda_range_changed(self):
        rough_peak_position = self.get_rough_peak_position()

        new_lambda_range = self.bragg_peak_range_ui.getRegion()
        left_index = find_nearest_index(self.o_roi.lambda_array, new_lambda_range[0])
        right_index = find_nearest_index(self.o_roi.lambda_array, new_lambda_range[1])

        if (rough_peak_position <= new_lambda_range[0]) or (rough_peak_position >= new_lambda_range[1]):
            self.rough_peak_index_position = np.mean([left_index, right_index])

        rough_peak_index_position = find_nearest_index(self.o_roi.lambda_array, rough_peak_position)

        self.ui.rough_lambda_peak_position_slider.setMinimum(left_index)
        self.ui.rough_lambda_peak_position_slider.setMaximum(right_index)
        self.ui.rough_lambda_peak_position_slider.setValue(rough_peak_index_position)

        self.number_of_files_to_exclude_slider_changed(0)

    def fit_pixel_clicked(self):
        o_fit = FitHandler(parent=self)
        o_fit.fit(mode='pixel')

    def fit_full_roi_clicked(self):
        o_fit = FitHandler(parent=self)
        o_fit.fit(mode='full')

    def apply_clicked(self):
        self.close()
