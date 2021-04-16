import os
import shutil


def make_or_reset_folder(folder_name):
    if os.path.exists(folder_name):
        shutil.rmtree(folder_name)
    os.makedirs(folder_name)


def remove_file_ending_by(list=None, ending=""):
    clean_list = [_file for _file in list if not _file.endswith(ending)]
    return clean_list
