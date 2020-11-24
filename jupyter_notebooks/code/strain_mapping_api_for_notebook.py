from ipywidgets import widgets
from IPython.core.display import HTML
from IPython.display import display, clear_output
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path, PurePath
import shutil

import ipywe.fileselector
import ipywe.fileselector
from NeuNorm.normalization import Normalization
from neutronbraggedge.experiment_handler import TOF, Experiment

from ToFImaging import reduction_tools
from jupyter_notebooks.code import detector_correction

DEBUG = True
DEBUG_PATH = "/Users/j35/IPTS/IPTS-strain-mapping/raw"


class StrainMappingAPIForNotebook:

    working_dir = DEBUG_PATH if DEBUG else "./"
    is_working_with_raw_data_default = False

    def __init__(self):
        pass

    def general_settings(self):
        box1 = widgets.HBox([widgets.Checkbox(value=self.is_working_with_raw_data_default,
                                              layout=widgets.Layout(width="120px")),
                             widgets.HTML("Are you working with <b>raw data</b> (not MCP corrected)?",
                                          layout=widgets.Layout(width="400px")),
                             ])
        self.is_mcp_corrected_ui = box1.children[0]
        display(box1)

    def mcp_correct_input_folders(self, input_folders=None):
        new_input_folders = []
        for _index_folder, _input_folder in enumerate(input_folders):
            print(f"MCP detector correction with folder {_index_folder+1}/{len(input_folders)} ... IN PROGRESS",
                  end="")
            _short_input_folder = str(PurePath(_input_folder).name) + "_corrected"
            _output_folder = str(Path(_input_folder).parent / _short_input_folder)
            StrainMappingAPIForNotebook.make_or_reset_folder(_output_folder)
            _shutter_count_file = glob.glob(_input_folder + "/*_ShutterCount.txt")[0]
            _shutter_time_file = glob.glob(_input_folder + "/*_ShutterTimes.txt")[0]
            _spectra_file = glob.glob(_input_folder + '/*_Spectra.txt')[0]
            df_meta = detector_correction.merge_meta_data(detector_correction.read_shutter_count(_shutter_count_file),
                                                          detector_correction.read_shutter_time(_shutter_time_file),
                                                          detector_correction.read_spectra(_spectra_file))
            o_norm = detector_correction.load_images(_input_folder)
            img_corrected = detector_correction.correct_images(o_norm,
                                                               df_meta,
                                                               skip_first_and_last=True)
            o_norm.data['sample']['data'] = img_corrected
            o_norm.export(folder=_output_folder, data_type='sample')
            new_input_folders.append(_output_folder)
            clear_output(wait=False)
            print(f"MCP detector correction with folder {_index_folder + 1}/{len(input_folders)} ... DONE")
        return new_input_folders

    def load_sample(self, input_folders):

        if self.is_mcp_corrected_ui.value:
            input_folders = self.mcp_correct_input_folders(input_folders=input_folders)
            clear_output(wait=False)
            display(HTML('<span style="font-size: 15px; color:blue">MCP detector correction: </span>'
                         '<span style="font-size: 15px; color:green">DONE</span>'))
            display(HTML('<span style="font-size: 15px; color:blue">List of corrected sample folders:</span>'))
            for _folder in input_folders:
                display(HTML('<span style="font-size: 15px; color:blue">- ' + _folder + '</span>'))

        if len(input_folders) == 1:
            # only one folder

            self.working_dir = input_folders[0]
            list_tif = glob.glob(input_folders[0] + "/*.tif")
            list_tif.sort()
            o_norm = Normalization()
            o_norm.load(file=list_tif, notebook=True)
            self.projections = np.array(o_norm.data['sample']['data'])
        else:
            # calculate mean of projections of the input folders

            list_projections = []
            for _folder in input_folders:
                list_tif = glob.glob(_folder[0] + "/*.tif")
                list_tif.sort()
                o_norm = Normalization()
                o_norm.load(file=list_tif, notebook=True)
                projection = np.array(o_norm.data['sample']['data'])
                list_projections.append(projection)
            self.projections = reduction_tools.mean_of_tof_arrays(list_array_3d=list_projections)

        self.projections = self.projections.transpose(1, 2, 0)   # x, y, lambda

    def select_projections(self, next_method=None, instruction='Select data folder ...'):
        fsel_ui = ipywe.fileselector.FileSelectorPanel(instruction=instruction,
                                                       start_dir=self.working_dir,
                                                       type='directory',
                                                       next=next_method,
                                                       multiple=True)
        fsel_ui.show()

    def select_sample(self):
        next_method = self.load_sample
        self.select_projections(next_method=next_method,
                                instruction="Select sample data folder ...")


    def display_integrated_signal(self):
        plt.imshow(np.nanmean(self.projections, axis=2))
        plt.title('White beam transmission image'), plt.colorbar()
        plt.show()
        plt.close()

    def define_instrument(self):
        box2 = widgets.HBox([widgets.Label("dSD (m)",
                                           layout=widgets.Layout(width="150px")),
                             widgets.Text(str(16.08),
                                          layout=widgets.Layout(width='30%'))])

        box3 = widgets.HBox([widgets.Label("detector offset (microns)",
                                           layout=widgets.Layout(width="150px")),
                             widgets.Text(str(3700),
                                          layout=widgets.Layout(width='30%'))])

        vertical_box = widgets.VBox([box2, box3])
        display(vertical_box)

        self.dsd = box2.children[1]
        self.doff = box3.children[1]

    def select_spectra_file(self):
        ts_ui = ipywe.fileselector.FileSelectorPanel(instruction='Select time spectra file ...',
                                                     default_filter='ASCII',
                                                     filters={'ASCII': ['*.txt']},
                                                     start_dir=self.working_dir)
        ts_ui.show()

    def load_spectra_file(self, spectra_file):
        tof_handler = TOF(filename=spectra_file)
        exp_handler = Experiment(tof=tof_handler.tof_array,
                                 distance_source_detector_m=np.float(self.dsd.value),
                                 detector_offset_micros=np.float(self.doff.value))
        self.lambda_array = exp_handler.lambda_array * 1e10  # to be in Angstroms
        display(HTML('<span style="font-size: 20px; color:blue">Spectra file laoded!</span>'))

    def prepare_fitting_mask(self, mask=None, plot=False):
        # mask[mask<0.7]=0
        # mask[mask>0]=1

        if plot:
            plt.imshow(mask)
            plt.title('Sample mask')
            plt.colorbar()
            plt.show()
            plt.close()
        self.mask = mask

    def calculate_moving_average(self, plot=False):
        # moving average with custom kernel to increase neutron statistics
        # custom_kernel = np.zeros((10,10))
        # custom_kernel[:,3:7] = 1
        custom_kernel = np.ones((5, 5))

        if plot:
            plt.imshow(custom_kernel)
            plt.title('Moving average kernel (10x4)'),
            plt.show()
            plt.close()

        T_mavg = reduction_tools.moving_average_2D(self.projections,
                                                   custom_kernel=custom_kernel)
        self.T_mavg = T_mavg

    @staticmethod
    def make_or_reset_folder(folder_name):
        if os.path.exists(folder_name):
            shutil.rmtree(folder_name)
        os.makedirs(folder_name)
