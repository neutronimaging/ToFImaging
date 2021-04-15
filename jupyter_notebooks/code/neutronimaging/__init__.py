# type: ignore
__description__ = "Neutron Imaging Reduction"
__url__ = "https://code.ornl.gov/sns-hfir-scse/imaging/neutronimaging"

__author__ = "C.Zhang"
__email__ = "zhangc@ornl.gov"

# __license__ = "GNU GENERAL PUBLIC LICENSE"

from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions
