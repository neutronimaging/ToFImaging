import os
from qtpy.QtWidgets import QMainWindow
from jupyter_notebooks.code import load_ui

ANGSTROMS = u"\u212B"


class Step4SettingsHandler(QMainWindow):

    def __init__(self, parent=None):
        super(QMainWindow, self).__init__()
        self.parent = parent
        ui_full_path = os.path.join(os.path.dirname(__file__), os.path.join('ui', 'ui_step4_settings.ui'))
        self.ui = load_ui(ui_full_path, baseinstance=self)
        self.setWindowTitle("Step4 - Settings")
        self.init_widgets()

    def init_widgets(self):
        config = self.parent.step4_config

        self.ui.position_label.setText(ANGSTROMS)
        self.ui.width_label.setText(ANGSTROMS)
        self.ui.position_from_label.setText(ANGSTROMS)
        self.ui.position_to_label.setText(ANGSTROMS)
        self.ui.width_from_label.setText(ANGSTROMS)
        self.ui.width_to_label.setText(ANGSTROMS)

        self.ui.estimated_bragg_edge_position_lineEdit.setText("{:.2f}".format(config[
                                                                                  'estimated_bragg_edge_position_value']))
        self.ui.estimated_bragg_edge_width_lineEdit.setText("{:.2f}".format(config['estimated_bragg_edge_width_value']))
        self.ui.estimated_bragg_edge_height_lineEdit.setText("{:.2f}".format(config[
                                                                                'estimated_bragg_edge_height_value']))

        position_from = config['estimated_bragg_edge_position_range'][0]
        position_to = config['estimated_bragg_edge_position_range'][1]
        self.ui.boundary_position_from_lineEdit.setText("{:.2f}".format(position_from))
        self.ui.boundary_position_to_lineEdit.setText("{:.2f}".format(position_to))

        width_from = config['estimated_bragg_edge_width_range'][0]
        width_to = config['estimated_bragg_edge_width_range'][1]
        self.ui.boundary_width_from_lineEdit.setText("{:.2f}".format(width_from))
        self.ui.boundary_width_to_lineEdit.setText("{:.2f}".format(width_to))

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

    def save_settings(self):
        pass

    def validate_changes_clicked(self):
        self.save_settings()
        self.close()

    def automatic_masking_clicked(self, status):
        self.ui.threshold_groupBox.setEnabled(status)

    def perform_filtering_algorithm_clicked(self, status):
        self.ui.perform_filtering_algorithm_groupBox.setEnabled(status)

    def interpolation_factor_clicked(self, status):
        self.ui.interpolation_factor_lineEdit.setEnabled(status)

    def closeEvent(self, event=None):
        self.parent.step4_settings_ui = None
        self.close()
