import os
from qtpy.QtWidgets import QMainWindow
from JupyterNotebooks.code import load_ui


class Step3SettingsHandler(QMainWindow):

    def __init__(self, parent=None):
        super(QMainWindow, self).__init__()
        self.parent = parent
        ui_full_path = os.path.join(os.path.dirname(__file__), os.path.join('ui', 'ui_step3_settings.ui'))
        self.ui = load_ui(ui_full_path, baseinstance=self)
        self.setWindowTitle("Step3 - Settings")
        self.init_widgets()

    def init_widgets(self):
        if self.parent.debugging_mode:
            bragg_peak_estimated = "3.5"
        else:
            bragg_peak_estimated = self.parent.o_roi.lambda_array[
                self.parent.ui.rough_lambda_peak_position_slider.value()]
        self.ui.estimated_position_bragg_edge_label.setText("{:.3f}".format(bragg_peak_estimated))
        self.ui.estimated_units_label.setText(u"\u212B")

        config = self.parent.step3_config
        self.ui.cross_section_mode_radioButton.setChecked(config['is_cross_section_mode'])
        self.ui.automatic_masking_checkBox.setChecked(config['is_automatic_masking'])
        self.ui.low_threshold_lineEdit.setText(str(config['threshold_low']))
        self.ui.high_threshold_lineEdit.setText(str(config['threshold_high']))
        self.ui.perform_filtering_algorithm_checkBox.setChecked(config['is_perform_savitzky_golay_filtering'])
        self.ui.window_size_lineEdit.setText(str(config['window_size']))
        self.ui.order_number_lineEdit.setText(str(config['order_number']))
        self.ui.interpolation_factor_checkBox.setChecked(config['is_interpolation_factor'])
        self.ui.interpolation_factor_lineEdit.setText(str(config['interpolation_factor_value']))

        self.automatic_masking_clicked(config['is_automatic_masking'])
        self.perform_filtering_algorithm_clicked(config['is_perform_savitzky_golay_filtering'])
        self.interpolation_factor_clicked(config['is_interpolation_factor'])

    def automatic_masking_clicked(self, status):
        self.ui.threshold_groupBox.setEnabled(status)

    def perform_filtering_algorithm_clicked(self, status):
        self.ui.perform_filtering_algorithm_groupBox.setEnabled(status)

    def interpolation_factor_clicked(self, status):
        self.ui.interpolation_factor_lineEdit.setEnabled(status)

    def save_settings(self):
        step3_config = {'is_automatic_masking': self.ui.automatic_masking_checkBox.isChecked(),
                        'threshold_low': str(self.ui.low_threshold_lineEdit.text()),
                        'threshold_high': str(self.ui.high_threshold_lineEdit.text()),
                        'is_perform_savitzky_golay_filtering':
                            self.ui.perform_filtering_algorithm_checkBox.isChecked(),
                        'window_size': str(self.ui.window_size_lineEdit.text()),
                        'order_number': str(self.ui.order_number_lineEdit.text()),
                        'is_interpolation_factor': self.ui.interpolation_factor_checkBox.isChecked(),
                        'interpolation_factor_value': str(self.ui.interpolation_factor_lineEdit.text()),
                        'is_cross_section_mode': self.ui.cross_section_mode_radioButton.isChecked()}
        self.parent.step3_config = step3_config

    def validate_changes_clicked(self):
        self.save_settings()
        self.close()

    def closeEvent(self, event=None):
        self.parent.step3_settings_ui = None
        self.close()