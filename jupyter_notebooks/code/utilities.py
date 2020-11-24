from collections import Counter
import os
import glob
import shutil


def retrieve_list_of_most_dominant_extension_from_folder(folder='', files=None):
    """
    This will return the list of files from the most dominant file extension found in the folder
    as well as the most dominant extension used
    """

    if folder:
        list_of_input_files = glob.glob(os.path.join(folder, '*'))
    else:
        list_of_input_files = files

    list_of_input_files.sort()
    list_of_base_name = [os.path.basename(_file) for _file in list_of_input_files]

    # work with the largest common file extension from the folder selected

    counter_extension = Counter()
    for _file in list_of_base_name:
        [_base, _ext] = os.path.splitext(_file)
        counter_extension[_ext] += 1

    dominant_extension = ''
    dominant_number = 0
    for _key in counter_extension.keys():
        if counter_extension[_key] > dominant_number:
            dominant_extension = _key
            dominant_number = counter_extension[_key]

    list_of_input_files = glob.glob(os.path.join(folder, '*' + dominant_extension))
    list_of_input_files.sort()

    list_of_input_files = [os.path.abspath(_file) for _file in list_of_input_files]

    return [list_of_input_files, dominant_extension]


def make_or_reset_folder(folder_name):
    if os.path.exists(folder_name):
        shutil.rmtree(folder_name)
    os.makedirs(folder_name)
