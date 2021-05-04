from JupyterNotebooks.code.parent import Parent


class EventHandler(Parent):

    def check_status_of_fit_buttons(self):

        # we need to make sure the pixel selected is inside one of the ROI
        # x_pixel, y_pixel = self.parent.pixel_marker['x'], self.parent.pixel_marker['y']
        y_pixel, x_pixel = self.parent.pixel_marker['x'], self.parent.pixel_marker['y']

        list_roi = self.parent.o_roi.list_roi
        for _index_roi in list_roi.keys():
            _roi = list_roi[_index_roi]
            _x0 = _roi['x0']
            _y0 = _roi['y0']
            _x1 = _roi['x1']
            _y1 = _roi['y1']

            if (x_pixel >= _x0) and (x_pixel <= _x1) and \
                    (y_pixel >= _y0) and (y_pixel <= _y1):
                self.parent.ui.step3_fit_pixel_pushButton.setEnabled(True)
                self.parent.ui.statusbar.showMessage("")
                return

        self.parent.ui.step3_fit_pixel_pushButton.setEnabled(False)
        self.parent.ui.statusbar.showMessage("Pixel must be inside one of the ROI selected!")
        self.parent.ui.statusbar.setStyleSheet("color: red")
