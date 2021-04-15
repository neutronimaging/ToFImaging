from pathlib import Path
from astropy.io import fits
import numpy as np
import glob
import pytest

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
        assert pytest.approx(mean_of_tof_arrays_returned[0][0][0], 1e-5) == 1.1293992


class TestMovingAverage1D:

    def test_wrong_input_raises_errors(self):
        input_array = np.ones((3, 3))
        with pytest.raises(ValueError):
            reduction_tools.moving_average_1D(input_array=input_array)

        with pytest.raises(ValueError):
            reduction_tools.moving_average_1D()

    def test_using_kernel_size(self):
        input_array = np.array([1, 10, 11, 9, 2, 1, 1, 30, 15, 0])
        output_array = reduction_tools.moving_average_1D(input_array=input_array, kernel_size=1)
        assert len(input_array) == len(output_array)
        for _input, _output in zip(input_array, output_array):
            assert _input == _output

        input_array = np.array([1, 10, 11, 9, 2, 1, 1, 30, 15, 0])
        output_array = reduction_tools.moving_average_1D(input_array=input_array, kernel_size=2)
        expected_array = np.array([1 / 2, 11 / 2, 21 / 2, 20 / 2, 11 / 2, 3 / 2, 2 / 2, 31 / 2, 45 / 2, 15 / 2])
        assert len(output_array) == len(expected_array)
        for _input, _output in zip(expected_array, output_array):
            assert _input == _output

        input_array = np.array([1, 10, 11, 9, 2, 1, 1, 30, 15, 0])
        output_array = reduction_tools.moving_average_1D(input_array=input_array, kernel_size=3)
        expected_array = np.array([(1 + 10) / 3, (1 + 10 + 11) / 3, (10 + 11 + 9) / 3, (11 + 9 + 2) / 3,
                                   (9 + 2 + 1) / 3, (2 + 1 + 1) / 3, (1 + 1 + 30) / 3,
                                   (1 + 30 + 15) / 3, (30 + 15 + 0) / 3, (15 + 0) / 3])
        assert len(output_array) == len(expected_array)
        for _input, _output in zip(expected_array, output_array):
            assert _input == _output

    def test_using_custom_kernel(self):
        input_array = np.array([1, 10, 11, 9, 2, 1, 1, 30, 15, 0])
        output_array = reduction_tools.moving_average_1D(input_array=input_array, custom_kernel=np.array([1, 1]))
        expected_array = np.array([1 / 2, 11 / 2, 21 / 2, 20 / 2, 11 / 2, 3 / 2, 2 / 2, 31 / 2, 45 / 2, 15 / 2])
        assert len(output_array) == len(expected_array)
        print(expected_array)
        for _input, _output in zip(expected_array, output_array):
            assert _input == _output

        input_array = np.array([1, 10, 11, 9, 2, 1, 1, 30, 15, 0])
        output_array = reduction_tools.moving_average_1D(input_array=input_array, custom_kernel=np.array([1, 1, 1]))
        expected_array = np.array([(1 + 10) / 3, (1 + 10 + 11) / 3, (10 + 11 + 9) / 3, (11 + 9 + 2) / 3,
                                   (9 + 2 + 1) / 3, (2 + 1 + 1) / 3, (1 + 1 + 30) / 3,
                                   (1 + 30 + 15) / 3, (30 + 15 + 0) / 3, (15 + 0) / 3])
        assert len(output_array) == len(expected_array)
        for _input, _output in zip(expected_array, output_array):
            assert _input == _output

        input_array = np.array([1, 10, 11, 9, 2, 1, 1, 30, 15, 0])
        output_array = reduction_tools.moving_average_1D(input_array=input_array, custom_kernel=np.array([1, 2, 1]))
        expected_array = np.array([(2 * 1 + 10) / 4, (1 + 2 * 10 + 11) / 4, (10 + 2 * 11 + 9) / 4, (11 + 2 * 9 + 2) / 4,
                                   (9 + 2 * 2 + 1) / 4, (2 + 2 * 1 + 1) / 4, (1 + 2 * 1 + 30) / 4,
                                   (1 + 2 * 30 + 15) / 4, (30 + 2 * 15 + 0) / 4, (15 + 2 * 0) / 4])
        assert len(output_array) == len(expected_array)
        for _input, _output in zip(expected_array, output_array):
            assert _input == _output


class TestMovingAverage2D:

    def test_wrong_input_raises_errors(self):
        input_array = np.array([1, 10, 11, 9, 2, 1, 1, 30, 15, 0])
        with pytest.raises(ValueError):
            reduction_tools.moving_average_2D(input_array=input_array)

        input_array = np.ones((3, 4, 5, 6))
        with pytest.raises(ValueError):
            reduction_tools.moving_average_2D(input_array=input_array)

    def test_using_custom_kernel(self):
        input_array = np.array([[[1, 10, 11], [9, 2, 1], [1, 30, 15]],
                               [[1, 10, 11], [9, 2, 1], [1, 30, 15]],
                               [[1, 10, 11], [9, 2, 1], [1, 30, 15]],
                               [[1, 10, 11], [9, 2, 1], [1, 30, 15]]])
        output_array = reduction_tools.moving_average_2D(input_array=input_array,
                                                         custom_kernel=np.ones((2, 2)))
        expected_array = [[[0.25, 2.75, 5.25],
                           [2.5, 5.5, 6],
                           [2.5, 10.5, 12]],
                          [[0.25, 2.75, 5.25],
                           [2.5, 5.5, 6.],
                           [2.5, 10.5, 12.]],
                          [[0.25, 2.75, 5.25],
                           [2.5, 5.5, 6],
                           [2.5, 10.5, 12]],
                          [[0.25, 2.75, 5.25],
                           [2.5, 5.5, 6],
                           [2.5, 10.5, 12.]]]
        assert np.shape(input_array) == np.shape(output_array)
        for _output_row, _expected_row in zip(output_array, expected_array):
            for _output, _expected in zip(_output_row, _expected_row):
                assert _output == pytest.approx(_expected, 1e-5)

    def test_using_custom_kernel_with_tof_data(self):
        # x:3, y:3 and tof:4
        input_array = np.array([[[1, 10, 11], [9, 2, 1], [1, 30, 15]],
                                [[1, 10, 11], [9, 2, 1], [1, 30, 15]],
                                [[1, 10, 11], [9, 2, 1], [1, 30, 15]],
                                [[1, 10, 11], [9, 2, 1], [1, 30, 15]]])
        input_array = input_array.transpose(1, 2, 0)
        output_array = reduction_tools.moving_average_2D(input_array=input_array,
                                                         custom_kernel=np.ones((3, 3)))
        expected_array = np.array([[[2.44, 3.67, 3.67, 2.44],
                                   [4.89, 7.33, 7.33, 4.89],
                                   [4.67, 7., 7., 4.67]],
                                  [[2.44, 3.67, 3.67, 2.44],
                                   [2.67, 4., 4., 2.67],
                                   [0.67, 1, 1, 0.67]],
                                  [[6.89, 10.33, 10.33, 6.89],
                                   [10.22, 15.33, 15.33, 10.22],
                                   [10, 15, 15, 10]]])
        assert np.shape(input_array) == np.shape(output_array)
        for _output_tof, _expected_tof in zip(output_array, expected_array):
            for _output_row, _expected_row in zip(_output_tof, _expected_tof):
                for _output, _expected in zip(_output_row, _expected_row):
                    assert _output == pytest.approx(_expected, 1e-2)

    def test_using_box_kernel(self):
        input_array = np.array([[1, 10, 11], [9, 2, 1], [1, 30, 15]])
        with pytest.raises(ValueError):
            reduction_tools.moving_average_2D(input_array=input_array,
                                              box_kernel=[1, 2, 3])

        input_array = np.array([[[1, 10, 11], [9, 2, 1], [1, 30, 15]],
                                [[1, 10, 11], [9, 2, 1], [1, 30, 15]],
                                [[1, 10, 11], [9, 2, 1], [1, 30, 15]],
                                [[1, 10, 11], [9, 2, 1], [1, 30, 15]]])
        output_array = reduction_tools.moving_average_2D(input_array=input_array,
                                                         box_kernel=[2, 2])
        expected_array = [[[0.25, 2.75, 5.25],
                           [2.5, 5.5, 6],
                           [2.5, 10.5, 12]],
                          [[0.25, 2.75, 5.25],
                           [2.5, 5.5, 6.],
                           [2.5, 10.5, 12.]],
                          [[0.25, 2.75, 5.25],
                           [2.5, 5.5, 6],
                           [2.5, 10.5, 12]],
                          [[0.25, 2.75, 5.25],
                           [2.5, 5.5, 6],
                           [2.5, 10.5, 12.]]]
        for _output_row, _expected_row in zip(output_array, expected_array):
            for _output, _expected in zip(_output_row, _expected_row):
                assert _output == pytest.approx(_expected, 1e-2)

        # input_array = np.array([[1, 10, 11], [9, 2, 1], [1, 30, 15]])
        # output_array = reduction_tools.moving_average_2D(input_array=input_array,
        #                                                  box_kernel=[2, 3])
        # expected_array = np.array([[1/4, 11/4, 21/4],
        #                           [10/4, (1+10+9+2)/4, (10+11+2+1)/4],
        #                           [(9+1)/4, (9+2+1+30)/4, (2+1+30+15)/4]])
        #
        # print(output_array)
        #
        # for _output_row, _expected_row in zip(output_array, expected_array):
        #     for _output, _expected in zip(_output_row, _expected_row):
        #         assert _output == pytest.approx(_expected, 1e-5)


class TestDataFiltering:

    def test_error_raised(self):
        input_array = np.array([1, 2, 3, 4, 5, 6])
        with pytest.raises(ValueError):
            reduction_tools.data_filtering(input_array)

        input_array = np.ones((3, 3, 4, 5))
        with pytest.raises(ValueError):
            reduction_tools.data_filtering(input_array)

        input_array = np.array([[1, 2, 3], [4, 5, 6]])
        with pytest.raises(ValueError):
            reduction_tools.data_filtering(input_array)

        input_array = np.array([[1, 2, 3], [4, 5, 6]])
        box_kernel = np.array([1])
        with pytest.raises(ValueError):
            reduction_tools.data_filtering(input_array,
                                           box_kernel=box_kernel)

        input_array = np.array([[1, 2, 3], [4, 5, 6]])
        gaussian_kernel = np.array([1])
        with pytest.raises(ValueError):
            reduction_tools.data_filtering(input_array,
                                           gaussian_kernel=gaussian_kernel)

        input_array = np.array([[1, 2, 3], [4, 5, 6]])
        box_kernel = np.array([1, 2, 3])
        with pytest.raises(ValueError):
            reduction_tools.data_filtering(input_array,
                                           box_kernel=box_kernel)

        input_array = np.array([[1, 2, 3], [4, 5, 6]])
        gaussian_kernel = np.array([1, 2, 3])
        with pytest.raises(ValueError):
            reduction_tools.data_filtering(input_array,
                                           gaussian_kernel=gaussian_kernel)

    # def test_3d_input_array_with_2d_box_kernel(self):
    #     input_array = np.array([[[1, 10, 11], [9, 2, 1], [1, 30, 15]],
    #                             [[1, 10, 11], [9, 2, 1], [1, 30, 15]],
    #                             [[1, 10, 11], [9, 2, 1], [1, 30, 15]],
    #                             [[1, 10, 11], [9, 2, 1], [1, 30, 15]]])
    #     # input_array = input_array.transpose(1, 2, 0)  # x, y, tof
    #     output_array = reduction_tools.data_filtering(input_array, box_kernel=[2, 2])
    #     print(output_array)
    #     assert False

    def test_2d_input_array_with_2d_box_kernel(self):
        input_array = np.array([[1, 10, 11], [9, 2, 1], [1, 30, 15]])
        output_array = reduction_tools.data_filtering(input_array, box_kernel=[2, 2])
        print(output_array)
        assert False

