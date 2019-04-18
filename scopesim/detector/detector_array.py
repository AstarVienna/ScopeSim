import warnings

from astropy.io import fits

from ..optics.effects.effects_utils import get_all_effects
from ..optics import effects as efs

from . import Detector


class DetectorArray:
    def __init__(self, **kwargs):
        self.meta = {}
        self.meta.update(kwargs)
        self.effects = []
        self.detectors = []
        self.latest_exposure = None

    def readout(self, image_plane, effects=[], **kwargs):

        # .. note:: Detector is what used to be called Chip
        #           DetectorArray is the old Detector

        # 1. make a series of Detectors for each row in a DetectorList object
        # 2. iterate through all Detectors, extract image from image_plane
        # 3. apply all effects (to all Detectors)
        # 4. add necessary header keywords

        # 5. Generate a HDUList with the ImageHDUs and any extras:
        # - add PrimaryHDU with meta data regarding observation in header
        # - add ImageHDUs
        # - add ASCIITableHDU with Effects meta data in final table extension

        self.effects += effects
        self.meta.update(kwargs)

        detector_list = get_detector_list(self.effects)
        self.detectors = [Detector(hdr, **self.meta)
                          for hdr in detector_list.detector_headers()]

        for detector in self.detectors:
            detector.extract(image_plane)

            for effect in self.effects:
                detector = effect.apply_to(detector, **self.meta)

        pri_hdu = make_primary_hdu(self.meta)
        effects_hdu = make_effects_hdu(self.effects)

        hdu_list = fits.HDUList([pri_hdu] +
                                [dtcr.hdu for dtcr in self.detectors] +
                                [effects_hdu])
        self.latest_exposure = hdu_list

        return self.latest_exposure


def make_primary_hdu(meta):
    prihdu = fits.PrimaryHDU()
    prihdu.header.update(meta)

    return prihdu


def make_effects_hdu(effects):
    # .. todo:: decide what goes into the effects table of meta data
    return fits.TableHDU()


def get_detector_list(effects):
    detector_list = get_all_effects(effects, efs.DetectorList)

    if len(detector_list) != 1:
        warnings.warn("None or more than one DetectorList found. Using the"
                      " first instance.{}".format(detector_list))

    return detector_list[0]
