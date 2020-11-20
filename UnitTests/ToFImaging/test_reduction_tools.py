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
        fits_file3 = str(fits_path / 'image3.fits')

        self.image1 = fits.getdata(fits_file1)
        self.image2 = fits.getdata(fits_file2)
        self.image3 = fits.getdata(fits_file3)

        image12_3d = np.array([self.image1, self.image2])
        self.image12_3d = image12_3d.transpose(1, 2, 0)

        image123_3d = np.array([self.image1, self.image2, self.image3])
        self.image123_3d = image123_3d.transpose(1, 2, 0)

    def test_mean_image(self):
        mean_image = average_image(self.image1)
        assert np.shape(mean_image) == np.shape(self.image1)
        assert (mean_image == self.image1).all

        mean_image = average_image(self.image12_3d)
        assert np.shape(mean_image) == np.shape(self.image1)
        expected_mean_image = np.zeros((4, 4))
        expected_mean_image[0:2, 0:2] = 5
        expected_mean_image[0:2, 2:] = 6
        expected_mean_image[2:, 0:2] = 7
        expected_mean_image[2:, 2:] = 8
        assert (expected_mean_image == mean_image).all()

    def test_median_image(self):
        mean_image = median_image(self.image1)
        assert np.shape(mean_image) == np.shape(self.image1)
        assert (mean_image == self.image1).all

        md_image = median_image(self.image123_3d)
        assert np.shape(md_image) == np.shape(self.image1)
        expected_md_image = np.zeros((4, 4))
        expected_md_image[0:2, 0:2] = 10
        expected_md_image[0:2, 2:] = 11
        expected_md_image[2:, 0:2] = 12
        expected_md_image[2:, 2:] = 13
        assert (expected_md_image == md_image).all()

    def test_combine_images(self):
        """assert that an array of images are correctly combined"""
        mean_image = combine_images(images=self.image12_3d, algorithm='mean')
        assert np.shape(mean_image) == np.shape(self.image1)
        expected_mean_image = np.zeros((4, 4))
        expected_mean_image[0:2, 0:2] = 5
        expected_mean_image[0:2, 2:] = 6
        expected_mean_image[2:, 0:2] = 7
        expected_mean_image[2:, 2:] = 8
        assert (expected_mean_image == mean_image).all()

