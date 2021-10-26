from pathlib import Path
from astropy.io import fits
import numpy as np
# import glob
import pytest

from NeuNorm.normalization import Normalization

from tofimaging import ReductionTools


class TestReductionTools:

    def setup_method(self):
        data_path = Path(__file__).parent
        self.images_path = Path(data_path).parent / 'test_data' / 'images'
        fits_file1 = str(self.images_path / 'image1.fits')
        fits_file2 = str(self.images_path / 'image2.fits')
        fits_file3 = str(self.images_path / 'image3.fits')

        self.image1 = fits.getdata(fits_file1)
        self.image2 = fits.getdata(fits_file2)
        self.image3 = fits.getdata(fits_file3)

        image12_3d = np.array([self.image1, self.image2])
        self.image12_3d = image12_3d.transpose(1, 2, 0)

        image123_3d = np.array([self.image1, self.image2, self.image3])
        self.image123_3d = image123_3d.transpose(1, 2, 0)

    def test_data_filtering_parameters_error(self):
        with pytest.raises(ValueError):
            ReductionTools.data_filtering()

        with pytest.raises(ValueError):
            my_signal = np.zeros((10))
            ReductionTools.data_filtering(mysignal=my_signal)

        with pytest.raises(ValueError):
            box_kernel = np.ones((10, 10))
            my_signal = np.zeros((10, 10, 10, 10))
            ReductionTools.data_filtering(mysignal=my_signal, kernel=box_kernel)

        with pytest.raises(ValueError):
            my_signal = np.zeros((10, 10))
            ReductionTools.data_filtering(mysignal=my_signal)

        with pytest.raises(ValueError):
            my_signal = np.zeros((10, 10))
            kernel = [3]
            ReductionTools.data_filtering(mysignal=my_signal, kernel=kernel)

        with pytest.raises(ValueError):
            my_signal = np.zeros((10, 10))
            kernel = [2,3,3]
            ReductionTools.data_filtering(mysignal=my_signal, kernel=kernel)
