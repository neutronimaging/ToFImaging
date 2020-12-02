class FitHandler:

    def __init__(self, parent=None):
        self.parent = parent

    def fit_pixel(self):
        # get pixel coordinates
        pixel_marker = self.parent.pixel_maker
        pixel = [pixel_marker['x'],
                 pixel_marker['y']]

        # get algorithm selected
        algorithm_selected = self.get_algorithm_selected()

        # get lambda range
        lambda_range = self.parent.bragg_peak_range_ui.getRegion()



    def get_algorithm_selected(self):
        if self.ui.gaussian_radioButton.isChecked():
            return 'gaussian'
        elif self.ui.advanced_radioButton.isChecked():
            return 'advanced'

        raise NotImplementedError