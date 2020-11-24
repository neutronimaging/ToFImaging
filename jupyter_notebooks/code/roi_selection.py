import os
from qtpy.QtWidgets import QMainWindow
from jupyter_notebooks.code import load_ui


class Interface(QMainWindow):

    def __init__(self, parent=None, working_dir="", o_bragg=None, spectra_file=None):
        super(Interface, self).__init__(parent)

        ui_full_path = os.path.join(os.path.dirname(__file__), os.path.join('ui', 'ui_roi_selection.ui'))
        self.ui = load_ui(ui_full_path, baseinstance=self)
        self.setWindowTitle("Region of Interest Tool")

    def update_table_roi(self):
        pass

    def add_roi_button_clicked(self):
        pass

    def remove_roi_button_clicked(self):
        pass

    def cancel_clicked(self):
        self.close()

    def apply_clicked(self):
        print('create mask!')
