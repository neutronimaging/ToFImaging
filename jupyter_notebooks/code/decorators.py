from qtpy.QtWidgets import QApplication
from qtpy import QtCore, QtGui


def wait_cursor(function):
    """
    Add a wait cursor during the running of the function
    """

    def wrapper(self, *args, **kwargs):
        QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        QtGui.QGuiApplication.processEvents()
        # function(*args, **kwargs)
        function(self)
        QApplication.restoreOverrideCursor()
        QtGui.QGuiApplication.processEvents()

    return wrapper
