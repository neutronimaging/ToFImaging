from pathlib import Path
from astropy.io import fits
import numpy as np
import glob

from NeuNorm.normalization import Normalization

from ToFImaging import reduction_tools


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

    def test_mean_image(self):
        mean_image = reduction_tools.average_image(self.image1)
        assert np.shape(mean_image) == np.shape(self.image1)
        assert (mean_image == self.image1).all

        mean_image = reduction_tools.average_image(self.image12_3d)
        assert np.shape(mean_image) == np.shape(self.image1)
        expected_mean_image = np.zeros((4, 4))
        expected_mean_image[0:2, 0:2] = 5
        expected_mean_image[0:2, 2:] = 6
        expected_mean_image[2:, 0:2] = 7
        expected_mean_image[2:, 2:] = 8
        assert (expected_mean_image == mean_image).all()

    def test_median_image(self):
        mean_image = reduction_tools.median_image(self.image1)
        assert np.shape(mean_image) == np.shape(self.image1)
        assert (mean_image == self.image1).all

        md_image = reduction_tools.median_image(self.image123_3d)
        assert np.shape(md_image) == np.shape(self.image1)
        expected_md_image = np.zeros((4, 4))
        expected_md_image[0:2, 0:2] = 10
        expected_md_image[0:2, 2:] = 11
        expected_md_image[2:, 0:2] = 12
        expected_md_image[2:, 2:] = 13
        assert (expected_md_image == md_image).all()

    # def test_weighted_average_images(self):
    #     weighted_image = reduction_tools.weighted_average_image(self.image1)
    #     assert np.shape(weighted_image) == np.shape(self.image1)
    #     assert (weighted_image == self.image1).all()
    #
    #     weighted_image = weighted_image(self.image123_3d)
    #     assert np.shape(weighted_image) == np.shape(self.image123_3d)

    def test_combine_images(self):
        """assert that an array of images are correctly combined"""
        mean_image = reduction_tools.combine_images(images=self.image12_3d, algorithm='mean')
        assert np.shape(mean_image) == np.shape(self.image1)
        expected_mean_image = np.zeros((4, 4))
        expected_mean_image[0:2, 0:2] = 5
        expected_mean_image[0:2, 2:] = 6
        expected_mean_image[2:, 0:2] = 7
        expected_mean_image[2:, 2:] = 8
        assert (expected_mean_image == mean_image).all()

        md_image = reduction_tools.combine_images(self.image123_3d, algorithm='median')
        assert np.shape(md_image) == np.shape(self.image1)
        expected_md_image = np.zeros((4, 4))
        expected_md_image[0:2, 0:2] = 10
        expected_md_image[0:2, 2:] = 11
        expected_md_image[2:, 0:2] = 12
        expected_md_image[2:, 2:] = 13
        assert (expected_md_image == md_image).all()

    def test_mean_of_tof_array(self):
        images_path = str(self.images_path)
        full_list_tif = np.array(glob.glob(images_path + '/*.tif'))
        print(full_list_tif)
        list_tif_1_index = [0, 1, 2, 4]
        list_tif_2_index = [1, 2, 3, 0]
        list_tif_3_index = [2, 3, 4, 1]

        list_tif_1 = full_list_tif[list_tif_1_index]
        list_tif_2 = full_list_tif[list_tif_2_index]
        list_tif_3 = full_list_tif[list_tif_3_index]

        o_norm = Normalization()
        data = []
        for _file in list_tif_1:
            o_norm.load(file=_file, notebook=True)
        data.append(o_norm.data['sample']['data'])
        del o_norm

        o_norm1 = Normalization()
        for _file in list_tif_2:
            o_norm1.load(file=_file, notebook=True)
        data.append(o_norm1.data['sample']['data'])
        del o_norm1

        o_norm2 = Normalization()
        for _file in list_tif_3:
            o_norm2.load(file=_file, notebook=True)
        data.append(o_norm2.data['sample']['data'])
        del o_norm2

        mean_of_tof_arrays_returned = reduction_tools.mean_of_tof_arrays(data)
        print()
        assert False


