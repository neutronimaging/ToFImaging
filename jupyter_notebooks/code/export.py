from qtpy.QtWidgets import QFileDialog

from NeuNorm.normalization import Normalization


class Export:

    def __init__(self, parent=None):
        self.parent = parent

    def run(self):
        working_dir = self.parent.o_api.working_dir
        export_folder = QFileDialog.getExistingDirectory(self.parent,
                                                         caption="Select output folder",
                                                         directory=working_dir)

        if export_folder:
            print('export data')
