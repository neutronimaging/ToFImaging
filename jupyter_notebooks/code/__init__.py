from qtpy.uic import loadUi

__all__ = ['load_ui']

def load_ui(ui_filename, baseinstance):
    return loadUi(ui_filename, baseinstance=baseinstance)
