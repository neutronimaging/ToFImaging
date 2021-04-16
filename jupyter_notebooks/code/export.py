from qtpy.QtWidgets import QFileDialog
from qtpy import QtGui
from pathlib import Path

from NeuNorm.normalization import Normalization

from jupyter_notebooks.code.utilities.file import make_or_reset_folder


class Export:

    def __init__(self, parent=None):
        self.parent = parent

    def run(self):
        working_dir = str(Path(self.parent.o_api.working_dir).parent)
        export_folder = QFileDialog.getExistingDirectory(self.parent,
                                                         caption="Select output folder",
                                                         directory=working_dir)

        QtGui.QGuiApplication.processEvents()

        if export_folder:
            normalize_projections = self.parent.normalize_projections
            list_sample_filename = self.parent.o_api.list_sample_projections_filename
            sample_folder_name = str(Path(list_sample_filename[0]).parent.name)
            output_file_name = str(Path(str(export_folder), sample_folder_name + "_prepared_data"))
            make_or_reset_folder(output_file_name)

            o_norm = Normalization()
            o_norm.load(data=normalize_projections)
            o_norm.data['sample']['file_name'] = list_sample_filename
            o_norm.export(output_file_name, data_type='sample')

            self.parent.ui.statusbar.showMessage("Export to {} ... Done!".format(output_file_name), 15000)
