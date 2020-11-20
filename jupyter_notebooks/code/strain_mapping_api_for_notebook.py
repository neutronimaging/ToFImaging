from ipywidgets import widgets
from IPython.core.display import HTML
from IPython.display import display
import os
import glob
import numpy as np
import matplotlib.pyplot as plt

import ipywe.fileselector
import ipywe.fileselector
from NeuNorm.normalization import Normalization
from neutronbraggedge.experiment_handler import TOF, Experiment

from ToFImaging import reduction_tools

DEBUG = True
DEBUG_PATH = "/Users/j35/IPTS/SNAP/Si_normalized/"


class StrainMappingAPIForNotebook:

    working_dir = DEBUG_PATH if DEBUG else "./"

    def __init__(self):
        pass

    def load_sample(self, input_folders):
        if len(input_folders) == 1:
            self.working_dir = input_folders[0]
            list_tif = glob.glob(input_folders[0] + "/*.tif")
            list_tif.sort()
            o_norm = Normalization()
            o_norm.load(file=list_tif, notebook=True)
            self.projections = np.array(o_norm.data['sample']['data'])
            self.projections = self.projections.transpose(1, 2, 0)   # x, y, lambda
        else:
            print(f"load all folders and combine them")

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
