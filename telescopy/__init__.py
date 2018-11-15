import warnings
from astropy.utils.exceptions import AstropyWarning

from . import commands
from . import data
from . import detector
from . import optics
from . import source
from . import simulation
from . import utils

from .version import __version__
from .default_keywords import *

# turn off warnings - interesting for development, but not for runtime
#warnings.simplefilter('ignore', UserWarning)   # user should see UserWarnings
# warnings.simplefilter('ignore', FutureWarning)
# warnings.simplefilter('ignore', RuntimeWarning)  # warnings for the developer
# warnings.simplefilter('ignore', category=AstropyWarning)

