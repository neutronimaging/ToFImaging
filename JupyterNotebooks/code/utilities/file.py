import os
import shutil
import ntpath
from pathlib import Path

def make_or_reset_folder(folder_name):
    if os.path.exists(folder_name):
        shutil.rmtree(folder_name)
    os.makedirs(folder_name)


def remove_file_ending_by(list=None, ending=""):
    clean_list = [_file for _file in list if not _file.endswith(ending)]
    return clean_list


def read_ascii(filename=''):
    '''return contain of an ascii file'''
    with open(filename, 'r') as f:
        text = f.read()
    return text


def write_ascii(text="", filename=''):
    with open(filename, 'w') as f:
        f.write(text)


def path_leaf(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)


def validate_path(path):
    """this method will check to see if the path exists, if it doesn't it will recursively check the parent"""
    path = Path(path)
    if os.path.exists(path):
        return path
    else:
        return validate_path(path.parent)
