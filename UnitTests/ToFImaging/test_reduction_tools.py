from pathlib import Path
from astropy.io import fits

from ToFImaging.reduction_tools import *


class TestReductionTools:

    def setup_method(self):
        data_path = Path(__file__).parent
        fits_path = str(Path(data_path) / 'test_data' / 'images')
        fits_file1 = fits_path + 'image1.fits'
        fits_file2 = fits_path + 'image2.fits'

        self.image1 = fits.getdata(fits_file1)
        self.image2 = fits.getdata(fits_file2)

    def test_averageimage(self):
        mean_image = averageimage([self.image1])

        print(mean_image)
        assert False

    def test_combine_images(self):
        """assert that an array of images are correctly combined"""
        fits_file1 = self.fits_path + 'image1.fits'
        fits_file2 = self.fits_path + 'image2.fits'

