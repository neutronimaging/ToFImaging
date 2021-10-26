from ipywidgets import widgets
from IPython.core.display import HTML
from IPython.display import display, clear_output
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path, PurePath
import shutil
import logging

import ipywe.fileselector
import ipywe.fileselector
from NeuNorm.normalization import Normalization
from neutronbraggedge.experiment_handler import TOF, Experiment
from JupyterNotebooks.code.neutronimaging.detector_correction import read_shutter_count
from JupyterNotebooks.code.neutronimaging.detector_correction import read_shutter_time
from JupyterNotebooks.code.neutronimaging.detector_correction import read_spectra
from JupyterNotebooks.code.neutronimaging.detector_correction import merge_meta_data
from JupyterNotebooks.code.neutronimaging.detector_correction import load_images
from JupyterNotebooks.code.neutronimaging.detector_correction import correct_images
from JupyterNotebooks.code.neutronimaging.detector_correction import skipping_meta_data
from JupyterNotebooks.code.utilities.get import Get
from JupyterNotebooks.code.utilities import file as file_utilities

from tofimaging import ReductionTools as reduction_tools

DEBUG = True
DEBUG_PATH = "/Volumes/G-DRIVE/IPTS/IPTS-26171-testing_ibeatles/"


class StrainMappingAPIForNotebook:

    normalization_flag_ui = None

    working_dir = DEBUG_PATH if DEBUG else "./"
    is_working_with_raw_data_default = False
    do_you_want_to_normalize_the_data = False

    sample_projections = None     # x, y, lambda
    ob_projections = None         # x, y, lambda
    normalize_projections = None  # x, y, lambda

    list_sample_projections_filename = None
    list_ob_projections_filename = None

    tof_array = None
    lambda_array = None
    list_roi = None

    message = []

    def __init__(self):
        pass

    def general_settings(self):
        box1 = widgets.HBox([widgets.Checkbox(value=self.is_working_with_raw_data_default,
                                              layout=widgets.Layout(width="200px")),
                             widgets.HTML("Are you working with <b>raw data</b> (i.e. not MCP corrected)?",
                                          layout=widgets.Layout(width="400px")),
                             ])
        self.is_mcp_corrected_ui = box1.children[0]

        box2 = widgets.HBox([widgets.Checkbox(value=self.do_you_want_to_normalize_the_data,
                                              layout=widgets.Layout(width="200px")),
                             widgets.Label("Do you want to normalize your data?",
                                           layout=widgets.Layout(width="400px")),
                                         ])
        vertical_box_1 = widgets.VBox([box1, box2])
        self.normalization_flag_ui = box2.children[0]
        display(vertical_box_1)

        box3 = widgets.HBox([widgets.Label("distance source-detector",
                                           layout=widgets.Layout(width="200px")),
                             widgets.Text(str(16.08),
                                          layout=widgets.Layout(width='10%')),
                             widgets.Label("m")])

        box4 = widgets.HBox([widgets.Label("detector offset",
                                           layout=widgets.Layout(width="200px")),
                             widgets.Text(str(3700),
                                          layout=widgets.Layout(width='10%')),
                             widgets.Label(u"\u00B5s")])
        vertical_box = widgets.VBox([box3, box4])
        display(vertical_box)

        self.dsd = box3.children[1]
        self.doff = box4.children[1]

    def load_files(self, input_folders=None, data_type='sample'):
        if self.is_mcp_corrected_ui.value:
            input_folders = StrainMappingAPIForNotebook.mcp_correct_input_folders(input_folders=input_folders)
            clear_output(wait=False)
            display(HTML('<span style="font-size: 15px; color:blue">MCP detector correction: </span>'
                         '<span style="font-size: 15px; color:green">DONE</span>'))
            display(HTML('<span style="font-size: 15px; color:blue">List of corrected ' + data_type + ' folders:</span>'))
            for _folder in input_folders:
                display(HTML('<span style="font-size: 15px; color:blue">- ' + _folder + '</span>'))

        if len(input_folders) == 1:
            # only one folder
            self.working_dir = input_folders[0]
            list_files, ext = Get.list_of_most_dominant_extension_from_folder(folder=input_folders[0])
            list_files = file_utilities.remove_file_ending_by(list=list_files, ending="_SummedImg.fits")
            o_norm = Normalization()
            if data_type == 'sample':
                o_norm.load(file=list_files, notebook=True)
                projections = np.array(o_norm.data['sample']['data'])
                self.list_sample_projections_filename = o_norm.data['sample']['file_name']
            elif data_type == 'ob':
                o_norm.load(file=list_files, notebook=True, data_type='ob')
                projections = np.array(o_norm.data['ob']['data'])
                self.list_ob_projections_filename = o_norm.data['ob']['file_name']
        else:
            # calculate mean of projections of the input folders
            list_projections = []
            for _folder in input_folders:
                list_files, ext = Get.list_of_most_dominant_extension_from_folder(folder=_folder)
                o_norm = Normalization()
                o_norm.load(file=list_files, notebook=True)
                projection = np.array(o_norm.data['sample']['data'])
                list_projections.append(projection)
            projections = reduction_tools.mean_of_tof_arrays(list_array_3d=list_projections)

        if data_type == 'sample':
            # self.sample_projections = projections.transpose(1, 2, 0)  # x, y, lambda
            self.sample_projections = projections
            self.locate_and_load_spectra_file(input_folder=input_folders[0])
        elif data_type == 'ob':
            self.ob_projections = projections
            # self.ob_projections = projections.transpose(1, 2, 0)  # x, y, lambda
        else:
            raise NotImplementedError("Data type not implemented!")

    def load_data(self, input_folders):
        self.load_sample(input_folders=input_folders)
        self.load_ob(input_folders=input_folders)

    def load_sample(self, input_folders=None):
        self.load_files(input_folders=input_folders, data_type='sample')
        display(HTML('<span style="font-size: 15px; color:blue">Sample loaded successfully!</span>'))
        self.select_ob()

    def load_ob(self, input_folders=None):
        if self.normalization_flag_ui.value:
            self.load_files(input_folders=input_folders, data_type='ob')
            display(HTML('<span style="font-size: 15px; color:blue">OB loaded successfully!</span>'))
        else:
            display(HTML('<span style="font-size: 15px; color:blue">Normalizaton: OFF!</span>'))

    def locate_and_load_spectra_file(self, input_folder=None):
        spectra_file_name = glob.glob(input_folder + "/*_Spectra.txt")
        if not spectra_file_name:
            self.locate_spectra_file(input_folder=input_folder)
        else:
            self.load_spectra_file(spectra_file_name[0])

    def locate_spectra_file(self, input_folder=None):
        sfil_ui = ipywe.fileselector.FileSelectorPanel(instruction="Select Time Spectra File ...",
                                                       start_dir=input_folder,
                                                       multiple=False,
                                                       next=self.load_spectra_file)
        sfil_ui.show()

    def load_spectra_file(self, spectra_file):
        # make sure the name is right
        base_file_name = str(PurePath(spectra_file).name)
        if not ("_Spectra.txt" in base_file_name):
            self.locate_spectra_file(input_folder=str(Path(spectra_file).parent))
        else:
            tof_handler = TOF(filename=spectra_file)
            exp_handler = Experiment(tof=tof_handler.tof_array,
                                     distance_source_detector_m=np.float(self.dsd.value),
                                     detector_offset_micros=np.float(self.doff.value))
            self.lambda_array = exp_handler.lambda_array * 1e10  # to be in Angstroms
            self.tof_array = tof_handler.tof_array
            display(HTML('<span style="font-size: 15px; color:blue">Spectra file has been loaded successfully!</span>'))

    def select_projections(self, next_method=None, instruction='Select data folder ...'):
        fsel_ui = ipywe.fileselector.FileSelectorPanel(instruction=instruction,
                                                       start_dir=self.working_dir,
                                                       type='directory',
                                                       next=next_method,
                                                       multiple=True)
        fsel_ui.show()

    def select_data(self):
        self.select_sample()

    def select_sample(self):
        next_method = self.load_sample
        self.select_projections(next_method=next_method,
                                instruction="Select sample data folder ...")

    def select_ob(self):
        if self.normalization_flag_ui.value:
            next_method = self.load_ob
            self.working_dir = os.path.dirname(self.working_dir)
            self.select_projections(next_method=next_method,
                                    instruction="Select open beam folder ...")
        else:
            display(HTML('<span style="font-size: 15px; color:blue">No OB needed. Normalization process will be '
                         'skipped!</span>'))

    def preview_sample(self):
        self.display_integrated_signal(self.sample_projections)

    def display_integrated_signal(self, projections):
        plt.figure(figsize=[15, 15])
        plt.imshow(np.nanmean(projections, axis=2))
        plt.title('Integrated sample image'), plt.colorbar()
        plt.show()
        plt.close()

    def select_spectra_file(self):
        ts_ui = ipywe.fileselector.FileSelectorPanel(instruction='Select time spectra file ...',
                                                     default_filter='ASCII',
                                                     filters={'ASCII': ['*.txt']},
                                                     start_dir=self.working_dir)
        ts_ui.show()

    def display_message(self):
        clear_output(wait=False)
        message = self.message
        message = "\n".join(message)
        print(message)

    def calculate_moving_average(self, plot=False):
        # moving average with custom kernel to increase neutron statistics
        # custom_kernel = np.zeros((10,10))
        # custom_kernel[:,3:7] = 1
        self.message.append("Calculate moving average ... IN PROGRESS")
        self.display_message()
        custom_kernel = np.ones((5, 5))

        if plot:
            plt.imshow(custom_kernel)
            plt.title('Moving average kernel (10x4)'),
            plt.show()
            plt.close()

        T_mavg = reduction_tools.moving_average_2D(self.normalize_projections,
                                                   custom_kernel=custom_kernel)
        self.T_mavg = T_mavg
        self.message[-1] = "Calculate moving average ... Done"
        self.display_message()

    def calculate_mask(self):

        self.message.append("Calculate mask ... IN PROGRESS")
        self.display_message()

        list_roi = self.list_roi
        [height, width, _] = np.shape(self.T_mavg)
        mask = np.zeros((height, width))

        for _roi_key in list_roi.keys():
            _roi = list_roi[_roi_key]
            x0 = _roi['x0']
            y0 = _roi['y0']
            x1 = _roi['x1']
            y1 = _roi['y1']
            mask[y0:y1+1, x0:x1+1] = 1

        self.mask = mask

        self.message[-1] = "Calculate mask ... Done"
        self.display_message()

    @staticmethod
    def mcp_correct_input_folders(input_folders=None):
        new_input_folders = []
        for _index_folder, _input_folder in enumerate(input_folders):
            print(f"MCP detector correction with folder {_index_folder+1}/{len(input_folders)} ... IN PROGRESS",
                  end="")
            os.system('mcp_detector_correction')

            _short_input_folder = str(PurePath(_input_folder).name) + "_corrected"
            _output_folder = str(Path(_input_folder).parent / _short_input_folder)
            file_utilities.make_or_reset_folder(_output_folder)
            try:
                _shutter_count_file = glob.glob(_input_folder + "/*_ShutterCount.txt")[0]
                _shutter_time_file = glob.glob(_input_folder + "/*_ShutterTimes.txt")[0]
                _spectra_file = glob.glob(_input_folder + '/*_Spectra.txt')[0]
            except IndexError:
                raise IndexError(f"Unable to find _ShutterCount and/or _ShutterTimes files in {_input_folder}")

            df_meta = merge_meta_data(read_shutter_count(_shutter_count_file),
                                      read_shutter_time(_shutter_time_file),
                                      read_spectra(_spectra_file))
            o_norm = load_images(_input_folder)

            # running correction
            img_corrected = correct_images(o_norm,
                                           df_meta,
                                           skip_first_and_last=True)

            # exporting data
            o_norm.data['sample']['data'] = img_corrected
            o_norm.export(folder=_output_folder, data_type='sample')

            out_shutter_count = os.path.join(_output_folder,
                                             os.path.basename(_shutter_count_file))
            out_shutter_time = os.path.join(_output_folder,
                                            os.path.basename(_shutter_time_file))
            out_spectra_file = os.path.join(_output_folder,
                                            os.path.basename(_spectra_file))
            shutil.copyfile(_shutter_count_file, out_shutter_count)
            shutil.copyfile(_shutter_time_file, out_shutter_time)
            # handle proper spectra parsing
            skipping_meta_data(df_meta).to_csv(
                    out_spectra_file,
                    columns=['shutter_time', 'counts'],
                    index=False,
                    index_label=False,
                    header=False,
                    sep='\t',
            )

            new_input_folders.append(_output_folder)
            clear_output(wait=False)
            print(f"MCP detector correction with folder {_index_folder + 1}/{len(input_folders)} ... DONE")
        return new_input_folders
