from pathlib import Path
from astropy.io import fits
import numpy as np

from ToFImaging.reduction_tools import *


class TestReductionTools:

    def setup_method(self):
        data_path = Path(__file__).parent
        fits_path = Path(data_path).parent / 'test_data' / 'images'
        fits_file1 = str(fits_path / 'image1.fits')
        fits_file2 = str(fits_path / 'image2.fits')

        self.image1 = fits.getdata(fits_file1)
        self.image2 = fits.getdata(fits_file2)

    def test_average_image(self):
        mean_image = average_image(self.image1)
        assert np.shape(mean_image) == np.shape(self.image1)
        assert (mean_image == self.image1).all

        image_3d = np.array([self.image1, self.image2])
        image_3d = image_3d.transpose(1, 2, 0)

        mean_image = average_image(image_3d)
        assert np.shape(mean_image) == np.shape(self.image1)

    def test_combine_images(self):
        """assert that an array of images are correctly combined"""
        pass
        # fits_file1 = self.fits_path + 'image1.fits'
        # fits_file2 = self.fits_path + 'image2.fits'

