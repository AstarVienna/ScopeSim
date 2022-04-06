"""The Instrument data simulator for MICADO on the E-ELT"""

################################################################################
#                            TURN OFF WARNINGS                                 #
################################################################################
import sys
import logging
import warnings
import yaml
from astropy.utils.exceptions import AstropyWarning

warnings.simplefilter('ignore', UserWarning)
warnings.simplefilter('ignore', FutureWarning)
warnings.simplefilter('ignore', RuntimeWarning)     # warnings for the developer
warnings.simplefilter('ignore', category=AstropyWarning)
yaml.warnings({'YAMLLoadWarning': False})

################################################################################
#                         PACKAGE GLOBAL VARIABLES                             #
################################################################################

from . import rc



################################################################################
#                         SET BASIC LOGGING LEVEL                              #
################################################################################

root = logging.getLogger()
root.setLevel("DEBUG")            # DEBUG

if rc.__config__["!SIM.logging.log_to_file"] is True:
    file_path = rc.__config__["!SIM.logging.file_path"]
    write_mode = rc.__config__["!SIM.logging.file_open_mode"]
    file_handler = logging.FileHandler(file_path, write_mode)
    file_handler.setLevel(rc.__config__["!SIM.logging.file_level"])  # DEBUG
    formatter = logging.Formatter('%(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    root.addHandler(file_handler)

if rc.__config__["!SIM.logging.log_to_console"] is True:
    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(rc.__config__["!SIM.logging.console_level"])  # WARNING
    formatter = logging.Formatter('%(levelname)s - %(message)s')
    stdout_handler.setFormatter(formatter)
    root.addHandler(stdout_handler)


################################################################################
#                         IMPORT PACKAGE MODULES                               #
################################################################################

# Import all the modules to go under ScopeSim

from . import commands
from . import source
from . import optics
from . import detector
from . import effects
from . import server
from . import utils

# import specific classes from the modules to included in the global namespace

from .utils import bug_report
from .optics.optical_train import OpticalTrain
from .commands.user_commands import UserCommands
from .source.source import Source

from .server.database import (list_packages, download_package,
                              list_example_data,
                              download_example_data)

################################################################################
#                          VERSION INFORMATION                                 #
################################################################################

try:
    from .version import version as __version__
except ImportError:
    __version__ = "Version number is not available"
