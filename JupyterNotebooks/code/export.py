from qtpy.QtWidgets import QFileDialog
from qtpy import QtGui
from pathlib import Path
import logging

from NeuNorm.normalization import Normalization

from JupyterNotebooks.code.utilities import make_or_reset_folder


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
            logging.info(f"Exporting prepared data to {export_folder}")
            self.parent.ui.setEnabled(False)
            self.parent.ui.statusbar.showMessage("Export of images ... IN PROGRESS")
            self.parent.ui.statusbar.setStyleSheet("color: blue")
            QtGui.QGuiApplication.processEvents()

            normalize_projections = self.parent.normalize_projections
            list_sample_filename = self.parent.o_api.list_sample_projections_filename
            sample_folder_name = str(Path(list_sample_filename[0]).parent.name)
            output_file_name = str(Path(str(export_folder), sample_folder_name + "_prepared_data"))
            make_or_reset_folder(output_file_name)

            o_norm = Normalization()
            o_norm.load(data=normalize_projections)
            o_norm.data['sample']['file_name'] = list_sample_filename
            o_norm.export(output_file_name, data_type='sample')

            self.parent.ui.statusbar.showMessage("Export to {} - Done!".format(output_file_name), 15000)
            self.parent.ui.statusbar.setStyleSheet("color: green")
            QtGui.QGuiApplication.processEvents()
            self.parent.ui.setEnabled(True)

    def result(self):
        working_dir = str(Path(self.parent.o_api.working_dir).parent)
        export_folder = QFileDialog.getExistingDirectory(self.parent,
                                                         caption="Select output folder",
                                                         directory=working_dir)

        QtGui.QGuiApplication.processEvents()

        if export_folder:

            logging.info("Exporting result data:")
            list_sample_filename = self.parent.o_api.list_sample_projections_filename
            sample_folder_name = str(Path(list_sample_filename[0]).parent.name)
            output_file_name = str(Path(str(export_folder), sample_folder_name + "_results"))
            make_or_reset_folder(output_file_name)
            logging.info(f"-> export folder: {output_file_name}")

            result_dict = self.parent.full_fit_result
            list_key_with_data = ["edge_position",
                                  "edge_height",
                                  "edge_width",
                                  "edge_slope",
                                  "median_image"]
            for _key in list_key_with_data:
                _data = result_dict[_key]
                o_data = Normalization()
                o_data.load(data=_data)
                file_name = _key + ".tiff"
                o_data.data['sample']['file_name'][0] = file_name
                logging.info(f"-> exporting {_key} with file name {file_name}")
                o_data.export(output_file_name, data_type='sample')

            logging.info("Exporting result data ... DONE!")
            self.parent.ui.statusbar.showMessage("Exported results to {} ... Done!".format(output_file_name), 15000)
            self.parent.ui.statusbar.setStyleSheet("color: green")
