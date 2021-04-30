from qtpy.QtWidgets import QFileDialog
from qtpy import QtGui
from pathlib import Path
import logging
import numpy as np

from NeuNorm.normalization import Normalization

from jupyter_notebooks.code.display import Display
from jupyter_notebooks.code.prepare_data import PrepareData

from jupyter_notebooks.code.utilities.file import make_or_reset_folder


class Import:

    def __init__(self, parent=None):
        self.parent = parent

    def run(self):
        working_dir = str(Path(self.parent.o_api.working_dir).parent)
        import_folder = QFileDialog.getExistingDirectory(self.parent,
                                                         caption="Select input folder",
                                                         directory=working_dir)

        QtGui.QGuiApplication.processEvents()

        if import_folder:
            logging.info(f"Importing prepared data from {import_folder}")
            self.parent.ui.setEnabled(False)
            self.parent.ui.statusbar.showMessage("Import of images ... IN PROGRESS")
            self.parent.ui.statusbar.setStyleSheet("color: blue")
            QtGui.QGuiApplication.processEvents()

            logging.info(f"np.shape(untouched_sample_projections): "
                         f"{np.shape(self.parent.untouched_sample_projections)}")
            logging.info(f"type(untouched_sample_projections): {type(self.parent.untouched_sample_projections)}")

            o_norm = Normalization()
            o_norm.load(folder=import_folder)
            self.parent.normalize_projections = np.array(o_norm.data['sample']['data'])

            logging.info(f"np.shape(self.parent.normalize_projections): "
                         f"{np.shape(self.parent.normalize_projections)}")
            logging.info(f"type(normalize_projections): {type(self.parent.normalize_projections)}")

            o_display = Display(parent=self.parent)
            o_display.prepare_data_preview_image()
            o_display.fit_data_tab()

            o_prepare = PrepareData(parent=self.parent)
            o_prepare.calculate_mask()

            QtGui.QGuiApplication.processEvents()
            self.parent.ui.toolBox.setItemEnabled(1, True)
            self.parent.ui.export_prepared_data_pushButton.setEnabled(True)
            self.parent.ui.setEnabled(True)

            self.parent.ui.statusbar.showMessage("Import of images from {} ... Done!".format(import_folder),
                                                 15000)
            self.parent.ui.statusbar.setStyleSheet("color: green")
            logging.info(f"Imported prepared data - DONE!")
