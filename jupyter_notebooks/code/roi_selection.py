import os
from qtpy.QtWidgets import QMainWindow, QVBoxLayout, QTableWidgetItem, QTableWidgetSelectionRange
from qtpy.QtGui import QColor, QPen
from jupyter_notebooks.code import load_ui
import pyqtgraph as pg
import numpy as np
from collections import OrderedDict


class Interface(QMainWindow):

    list_roi = OrderedDict()
    default_roi = {'x0': 0, 'y0': 0, 'x1': 50, 'y1': 50, 'id': None}
    live_image = None
    roi_width = 0.01

    def __init__(self, parent=None, working_dir="", sample_projections=None, spectra_file=None):
        super(Interface, self).__init__(parent)

        self.sample_projections = sample_projections

        ui_full_path = os.path.join(os.path.dirname(__file__), os.path.join('ui', 'ui_roi_selection.ui'))
        self.ui = load_ui(ui_full_path, baseinstance=self)
        self.setWindowTitle("Region of Interest Tool")

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

    def display_image(self):
        sample_projections = self.sample_projections
        self.live_image = np.transpose(np.mean(sample_projections, axis=2))
        self.ui.image_view.setImage(self.live_image)

    def update_table_roi(self):
        pass

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
            return

        if row == (old_nbr_row - 1):
            row = new_nbr_row - 1

        _new_selection = QTableWidgetSelectionRange(row, 0, row, 3)
        self.ui.table_roi.setRangeSelected(_new_selection, True)

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
        _roi_id = pg.ROI([x0, y0], [width, height], pen=_pen, scaleSnap=True)
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
            region = roi_id.getArraySlice(self.live_image,
                                          self.ui.image_view.imageItem)

            x0 = region[0][0].start
            x1 = region[0][0].stop
            y0 = region[0][1].start
            y1 = region[0][1].stop

            _roi['x0'] = x0
            _roi['x1'] = x1
            _roi['y0'] = y0
            _roi['y1'] = y1

            list_roi[_row] = _roi

        self.list_roi = list_roi
        self.update_table_roi_ui()


    def cancel_clicked(self):
        self.close()

    def apply_clicked(self):
        print('create mask!')
