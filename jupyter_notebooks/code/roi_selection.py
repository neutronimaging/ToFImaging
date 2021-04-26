import os
from qtpy.QtWidgets import QMainWindow, QVBoxLayout, QTableWidgetItem, QTableWidgetSelectionRange
from qtpy.QtGui import QColor, QPen
from IPython.core.display import HTML
from IPython.display import display
from jupyter_notebooks.code import load_ui
import pyqtgraph as pg
import numpy as np
from collections import OrderedDict


class Interface(QMainWindow):

    list_roi = OrderedDict()
    default_roi = {'x0': 0, 'y0': 0, 'x1': 50, 'y1': 50, 'id': None}
    live_image = None
    roi_width = 0.01
    integrated_image_size = {'width': -1, 'height': -1}
    splitter_2_has_been_resized = False

    def __init__(self, parent=None, main_api=None):
        super(Interface, self).__init__(parent)

        self.sample_projections = main_api.sample_projections  # lambda, y, x
        # self.sample_projections_lambda_x_y = self.sample_projections.transpose(2, 0, 1)  # lambda, x, y
        self.sample_projections_lambda_x_y = self.sample_projections

        self.o_api = main_api
        self.tof_array = main_api.tof_array
        self.lambda_array = main_api.lambda_array

        ui_full_path = os.path.join(os.path.dirname(__file__), os.path.join('ui', 'ui_roi_selection.ui'))
        self.ui = load_ui(ui_full_path, baseinstance=self)
        self.setWindowTitle("Select region to fit!")

        self.init_widgets()
        self.display_image()

    def init_widgets(self):
        # pyqtgraphs
        self.ui.image_view = pg.ImageView()
        self.ui.image_view.ui.roiBtn.hide()
        self.ui.image_view.ui.menuBtn.hide()
        verti_layout1 = QVBoxLayout()
        verti_layout1.addWidget(self.ui.image_view)
        self.ui.widget_image.setLayout(verti_layout1)

        self.ui.plot_view = pg.PlotWidget()
        verti_layout2 = QVBoxLayout()
        verti_layout2.addWidget(self.ui.plot_view)
        self.ui.widget_plot.setLayout(verti_layout2)

        # table
        nbr_columns = self.ui.table_roi.columnCount()
        for _col in range(nbr_columns):
            self.ui.table_roi.setColumnWidth(_col, 70)

        # splitters
        self.ui.splitter_2.setSizes([400, 0])

        # x_axis buttons
        self.ui.tof_radioButton.setText(u"TOF (\u03BCs)")
        self.ui.lambda_radioButton.setText(u"\u03BB (\u212B)")

    def display_image(self):
        sample_projections = self.sample_projections
        self.live_image = np.mean(sample_projections, axis=0)
        # live_image = np.transpose(self.live_image)
        live_image = self.live_image
        self.ui.image_view.setImage(live_image)
        self.integrated_image_size['height'], self.integrated_image_size['width'] = np.shape(live_image)

    def check_roi_validity(self, value, x_axis=True):
        """Make sure the ROI selected or defined stays within the image size"""
        min_value = 0

        value = np.int(value)

        if x_axis:
            max_value = self.integrated_image_size['width']
        else:
            max_value = self.integrated_image_size['height']

        if value < 0:
            return min_value

        if value > max_value:
            return max_value

        return value

    def update_table_roi(self, item):
        """Using the table_roi_ui as reference, will update the list_roi dictionary"""
        self.ui.table_roi.blockSignals(True)

        nbr_row = self.ui.table_roi.rowCount()
        new_list_roi = OrderedDict()
        old_list_roi = self.list_roi
        for _row in np.arange(nbr_row):
            _roi = {}

            # checking that x0, y0, x1 and y1 stay within the range of the image
            _x0 = self.check_roi_validity(self._get_item_value(_row, 0))
            _y0 = self.check_roi_validity(self._get_item_value(_row, 1), x_axis=False)

            _x1 = self.check_roi_validity(self._get_item_value(_row, 2))
            _y1 = self.check_roi_validity(self._get_item_value(_row, 3), x_axis=False)

            # updating table content (in case some of the roi were out of scope
            self._set_item_value(_row, 0, _x0)
            self._set_item_value(_row, 1, _y0)
            self._set_item_value(_row, 2, _x1)
            self._set_item_value(_row, 3, _y1)

            _roi['x0'] = _x0
            _roi['y0'] = _y0
            _roi['x1'] = _x1
            _roi['y1'] = _y1
            _roi['id'] = old_list_roi[_row]['id']

            new_list_roi[_row] = _roi

        self.list_roi = new_list_roi
        self.update_image_view_item()
        self.update_plot_view()
        self.ui.table_roi.blockSignals(False)

    def clear_roi_on_image_view(self):
        """remove all ROI from image_view"""
        list_roi = self.list_roi
        for _row in list_roi.keys():
            _roi = list_roi[_row]
            roi_id = _roi['id']
            self.ui.image_view.removeItem(roi_id)

    def add_roi_button_clicked(self):
        self.clear_roi_on_image_view()

        self.ui.table_roi.blockSignals(True)
        _selection = self.ui.table_roi.selectedRanges()
        if _selection:
            row = _selection[0].topRow()
        else:
            row = 0

        # init new row with default value
        self.ui.table_roi.insertRow(row)
        _default_roi = self.default_roi

        _item = QTableWidgetItem(str(_default_roi['x0']))
        self.ui.table_roi.setItem(row, 0, _item)

        _item = QTableWidgetItem(str(_default_roi['y0']))
        self.ui.table_roi.setItem(row, 1, _item)

        _item = QTableWidgetItem(str(_default_roi['x1']))
        self.ui.table_roi.setItem(row, 2, _item)

        _item = QTableWidgetItem(str(_default_roi['y1']))
        self.ui.table_roi.setItem(row, 3, _item)

        # save new list_roi dictionary
        nbr_row = self.ui.table_roi.rowCount()
        list_roi = OrderedDict()
        for _row in np.arange(nbr_row):
            _roi = {}

            _x0 = self._get_item_value(_row, 0)
            _roi['x0'] = np.int(_x0)

            _y0 = self._get_item_value(_row, 1)
            _roi['y0'] = np.int(_y0)

            _x1 = self._get_item_value(_row, 2)
            _roi['x1'] = np.int(_x1)

            _y1 = self._get_item_value(_row, 3)
            _roi['y1'] = np.int(_y1)

            x0_int = int(_x0)
            y0_int = int(_y0)
            width_int = np.abs(x0_int - int(_x1))
            height_int = np.abs(y0_int - int(_y1))

            _roi_id = self.init_roi(x0=x0_int,
                                    y0=y0_int,
                                    width=width_int,
                                    height=height_int)
            _roi['id'] = _roi_id
            list_roi[_row] = _roi

        self.list_roi = list_roi
        self.ui.table_roi.blockSignals(False)
        self.check_add_remove_button_widgets_status()
        self.update_plot_view()

        if not _selection:
            _new_selection = QTableWidgetSelectionRange(0, 0, 0, 3)
            self.ui.table_roi.setRangeSelected(_new_selection, True)

    def check_add_remove_button_widgets_status(self):
        nbr_row = self.ui.table_roi.rowCount()
        if nbr_row > 0:
            self.ui.remove_roi_button.setEnabled(True)
        else:
            self.ui.remove_roi_button.setEnabled(False)

    def remove_row_entry(self, row):
        _roi_id = self.list_roi[row]['id']
        self.ui.image_view.removeItem(_roi_id)
        del self.list_roi[row]

        # rename row
        new_list_roi = {}
        new_row_index = 0
        for _previous_row_index in self.list_roi.keys():
            new_list_roi[new_row_index] = self.list_roi[_previous_row_index]
            new_row_index += 1
        self.list_roi = new_list_roi

    def remove_roi_button_clicked(self):
        self.ui.table_roi.blockSignals(True)

        _selection = self.ui.table_roi.selectedRanges()
        if not _selection:
            return

        row = _selection[0].topRow()
        old_nbr_row = self.ui.table_roi.rowCount()

        # remove entry from list of roi
        self.remove_row_entry(row)

        # update table of rois
        self.update_table_roi_ui()
        self.ui.table_roi.blockSignals(False)
        self.check_add_remove_button_widgets_status()

        # update selection
        new_nbr_row = self.ui.table_roi.rowCount()
        if new_nbr_row == 0:
            self.update_plot_view()
            return

        if row == (old_nbr_row - 1):
            row = new_nbr_row - 1

        _new_selection = QTableWidgetSelectionRange(row, 0, row, 3)
        self.ui.table_roi.setRangeSelected(_new_selection, True)
        self.update_plot_view()

    def clear_table(self):
        nbr_row = self.ui.table_roi.rowCount()
        for _row in np.arange(nbr_row):
            self.ui.table_roi.removeRow(0)

    def update_table_roi_ui(self):
        """Using list_roi as reference, repopulate the table_roi_ui"""

        self.ui.table_roi.blockSignals(True)
        list_roi = self.list_roi

        self.clear_table()

        _index_row = 0
        for _roi_key in list_roi.keys():
            _roi = list_roi[_roi_key]

            self.ui.table_roi.insertRow(_index_row)
            self._set_item_value(_index_row, 0, _roi['x0'])
            self._set_item_value(_index_row, 1, _roi['y0'])
            self._set_item_value(_index_row, 2, _roi['x1'])
            self._set_item_value(_index_row, 3, _roi['y1'])
            _index_row += 1

        self.ui.table_roi.blockSignals(False)

    def _set_item_value(self, row=0, column=0, value=-1):
        _item = QTableWidgetItem(str(value))
        self.ui.table_roi.setItem(row, column, _item)

    def _get_item_value(self, row, column):
        _item = self.ui.table_roi.item(row, column)
        if _item:
            return str(_item.text())
        else:
            return ''

    def init_roi(self, x0=0, y0=0, width=0, height=0):
        _color = QColor(62, 13, 244)
        _pen = QPen()
        _pen.setColor(_color)
        _pen.setWidthF(self.roi_width)
        _roi_id = pg.ROI([y0, x0], [height, width], pen=_pen, scaleSnap=True)
        _roi_id.addScaleHandle([1, 1], [0, 0])
        _roi_id.addScaleHandle([0, 0], [1, 1])
        self.ui.image_view.addItem(_roi_id)
        # add connection to roi
        _roi_id.sigRegionChanged.connect(self.roi_manually_moved)
        return _roi_id

    def roi_manually_moved(self):
        list_roi = self.list_roi

        for _row in list_roi.keys():

            _roi = list_roi[_row]

            roi_id = _roi['id']
            # region = roi_id.getArraySlice(np.transpose(self.live_image),
            #                               self.ui.image_view.imageItem)

            region = roi_id.getArraySlice(self.live_image,
                                          self.ui.image_view.imageItem)
            x0 = region[0][0].start
            x1 = region[0][0].stop
            y0 = region[0][1].start
            y1 = region[0][1].stop

            # y0 = region[0][0].start
            # y1 = region[0][0].stop
            # x0 = region[0][1].start
            # x1 = region[0][1].stop

            _roi['x0'] = x0
            _roi['x1'] = x1
            _roi['y0'] = y0
            _roi['y1'] = y1

            list_roi[_row] = _roi

        self.list_roi = list_roi
        self.update_table_roi_ui()
        self.update_plot_view()

    def update_image_view_item(self):
        self.clear_roi_on_image_view()

        list_roi = self.list_roi
        for _row in list_roi.keys():
            _roi = list_roi[_row]

            _x0 = np.int(_roi['x0'])
            _y0 = np.int(_roi['y0'])
            _x1 = np.int(_roi['x1'])
            _y1 = np.int(_roi['y1'])

            _width = np.abs(_x1 - _x0)
            _height = np.abs(_y1 - _y0)

            _roi_id = self.init_roi(x0=_x0, y0=_y0,
                                    width=_width, height=_height)
            _roi['id'] = _roi_id
            list_roi[_row] = _roi

        self.list_roi = list_roi

    def x_axis_units_changed(self):
        self.update_plot_view()

    def update_plot_view(self):

        self.ui.plot_view.clear()

        # only run if splitter never been resized
        if (not self.splitter_2_has_been_resized) and \
                (not self.ui.splitter_2.sizes()[1] == 0):
            self.splitter_2_has_been_resized = True

        if not self.splitter_2_has_been_resized:
            self.ui.splitter_2.setSizes([200, 200])

        list_roi = self.list_roi
        if not list_roi:
            self.x_axis_units_widgets_enabled(state=False)
        else:
            self.x_axis_units_widgets_enabled(state=True)

        y_axis = []
        total_number_of_pixels_in_roi = 0
        for _projection_index, _projection in enumerate(self.sample_projections_lambda_x_y):

            total_counts_of_roi = 0
            for _index_roi in self.list_roi.keys():
                _roi = self.list_roi[_index_roi]
                _x0 = _roi['x0']
                _y0 = _roi['y0']
                _x1 = _roi['x1']
                _y1 = _roi['y1']

                if _projection_index == 0:
                    total_number_of_pixels_in_roi += (_y1 - _y0 + 1) * (_x1 - _x0 + 1)
                # total_counts_of_roi += np.sum(_projection[_y0: _y1+1, _x0: _x1+1])
                total_counts_of_roi += np.sum(_projection[_x0: _x1 + 1, _y0: _y1 + 1])

            mean_counts_of_roi = total_counts_of_roi / total_number_of_pixels_in_roi
            y_axis.append(mean_counts_of_roi)

        x_axis = self.calculate_x_axis()
        self.ui.plot_view.plot(x_axis, y_axis)

    def calculate_x_axis(self):
        if self.ui.file_index_radioButton.isChecked():
            return np.arange(len(self.sample_projections_lambda_x_y))
        elif self.ui.tof_radioButton.isChecked():
            return self.tof_array * 1e6
        elif self.ui.lambda_radioButton.isChecked():
            return self.lambda_array
        else:
            return []

    def x_axis_units_widgets_enabled(self, state=True):
        list_ui = [self.ui.file_index_radioButton,
                   self.ui.tof_radioButton,
                   self.ui.lambda_radioButton]
        for _ui in list_ui:
            _ui.setEnabled(state)

    def cancel_clicked(self):
        if not self.list_roi:
            display(HTML('<span style="font-size: 15px; color:Red">You need to select at least 1 ROI to run the rest '
                         'of the notebook!</span>'))
        self.close()

    def apply_clicked(self):
        mask = np.zeros(np.shape(self.live_image))
        for _index_roi in self.list_roi.keys():
            _roi = self.list_roi[_index_roi]
            _x0 = _roi['x0']
            _y0 = _roi['y0']
            _x1 = _roi['x1']
            _y1 = _roi['y1']
            mask[_y0: _y1+1, _x0: _x1+1] = 1
        self.mask = mask

        self.close()
        display(HTML('<span style="font-size: 15px; color:green">You selected ' + str(len(self.list_roi)) +
                     ' ROI(s)</span>'))
