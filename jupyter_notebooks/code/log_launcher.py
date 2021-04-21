from qtpy.QtWidgets import QMainWindow
import os
from qtpy.QtGui import QIcon
from qtpy import QtGui
import logging

from . import load_ui
from .utilities.get import Get
from .utilities.file_utilities import read_ascii, write_ascii
from . import refresh_image


class LogLauncher:

    def __init__(self, parent=None):
        self.parent = parent

        if self.parent.log_id is None:
            log_id = Log(parent=self.parent)
            log_id.show()
            self.parent.log_id = log_id
        else:
            self.parent.log_id.activateWindow()
            self.parent.log_id.setFocus()


class Log(QMainWindow):

    def __init__(self, parent=None):
        self.parent = parent
        QMainWindow.__init__(self, parent=parent)
        ui_full_path = os.path.join(os.path.dirname(os.path.dirname(__file__)),
                                    os.path.join('ui',
                                                 'log.ui'))
        self.ui = load_ui(ui_full_path, baseinstance=self)
        self.setWindowTitle("Log")
        self.ui.log_text.setReadOnly(True)

        refresh_icon = QIcon(refresh_image)
        self.ui.refresh_pushButton.setIcon(refresh_icon)

        o_get = Get(parent=self.parent)
        self.log_file_name = o_get.get_log_file_name()
        self.loading_logging_file()

        # jump to end of file
        self.ui.log_text.moveCursor(QtGui.QTextCursor.End)

    def closeEvent(self, c):
        self.parent.log_id = None

    def loading_logging_file(self):
        try:
            log_text = read_ascii(self.log_file_name)
            self.ui.log_text.setPlainText(log_text)
            self.ui.log_text.moveCursor(QtGui.QTextCursor.End)
        except FileNotFoundError:
            self.ui.log_text.setPlainText("")

    def clear_clicked(self):
        if os.path.exists(self.log_file_name):
            write_ascii(text="", filename=self.log_file_name)
            logging.info("log file has been cleared by user")
        self.loading_logging_file()
