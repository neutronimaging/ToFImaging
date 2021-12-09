import pytest
from pathlib import Path
from JupyterNotebooks.code.utilities.file import validate_path


class TestFile:

    def setup_method(self):
        data_path = Path(__file__).parent.parent.parent
        self.existing_path = Path(data_path).parent / 'test_data' / 'images'
        self.partially_not_existing_path = Path(data_path).parent / 'test_data' / 'fake_folder'
        self.full_not_existing_path = Path("/fake_full_path") / 'fake_folder'

    def test_validate_good_path(self):
        path_expected = self.existing_path
        path_returned = validate_path(self.existing_path)
        assert (path_expected == path_returned)

    def test_validate_partially_bad_path(self):
        path_expected = Path(self.partially_not_existing_path).parent
        path_returned = validate_path(self.partially_not_existing_path)
        assert (path_expected == path_returned)

    def test_validate_full_bad_path(self):
        path_expected = Path("/")
        path_returned = validate_path(self.full_not_existing_path)
        assert (path_expected == path_returned)
