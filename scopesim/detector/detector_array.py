import logging

from astropy.io import fits

from .detector import Detector

from ..effects.effects_utils import get_all_effects
from .. import effects as efs
from .. import utils


class DetectorArray:
    def __init__(self, detector_list=None, **kwargs):
        self.meta = {}
        self.meta.update(kwargs)
        self.detector_list = detector_list
        self.array_effects = []
        self.dtcr_effects = []
        self.detectors = []
        self.latest_exposure = None

    def readout(self, image_planes, array_effects=[], dtcr_effects=[], **kwargs):
        """
        Read out the detector array into a FITS file

        Parameters
        ----------
        image_planes : list of ImagePlane objects
            The correct image plane is automatically chosen from the list

        array_effects : list of Effect objects
            A list of effects related to the detector array

        dtcr_effects : list of Effect objects
            A list of effects related to the detectors

        Returns
        -------
        self.latest_exposure : fits.HDUList

        """

        # .. note:: Detector is what used to be called Chip
        #           DetectorArray is the old Detector

        # 0. Select the relevant image plane to extract images from
        # 1. make a series of Detectors for each row in a DetectorList object
        # 2. iterate through all Detectors, extract image from image_plane
        # 3. apply all effects (to all Detectors)
        # 4. add necessary header keywords

        # 5. Generate a HDUList with the ImageHDUs and any extras:
        # - add PrimaryHDU with meta data regarding observation in header
        # - add ImageHDUs
        # - add ASCIITableHDU with Effects meta data in final table extension

        self.array_effects = array_effects
        self.dtcr_effects = dtcr_effects
        self.meta.update(kwargs)

        # 0. Get the image plane that corresponds to this detector array
        image_plane_id = self.detector_list.meta["image_plane_id"]
        image_plane = [implane for implane in image_planes if
                       implane.id == image_plane_id][0]

        # 0a. Apply detector array effects (apply to the entire image plane)
        for effect in self.array_effects:
            image_plane = effect.apply_to(image_plane, **self.meta)

        # 1. make a series of Detectors for each row in a DetectorList object
        self.detectors = [Detector(hdr, **self.meta)
                          for hdr in self.detector_list.detector_headers()]

        # 2. iterate through all Detectors, extract image from image_plane
        for detector in self.detectors:
            detector.extract_from(image_plane)

            # 3. apply all effects (to all Detectors)
            for effect in self.dtcr_effects:
                detector = effect.apply_to(detector)

            # 4. add necessary header keywords
            # .. todo: add keywords

        # 5. Generate a HDUList with the ImageHDUs and any extras:
        pri_hdu = make_primary_hdu(self.meta)

        # ..todo: effects_hdu unnecessary as long as make_effects_hdu does not do anything
        # effects_hdu = make_effects_hdu(self.array_effects + self.dtcr_effects)

        hdu_list = fits.HDUList([pri_hdu]
                                + [dtcr.hdu for dtcr in self.detectors])
                                # + [effects_hdu])
        self.latest_exposure = hdu_list

        return self.latest_exposure


def make_primary_hdu(meta):
    """Create the primary header from meta data"""
    new_meta = utils.stringify_dict(meta)
    prihdu = fits.PrimaryHDU()
    prihdu.header.update(new_meta)

    return prihdu


def make_effects_hdu(effects):
    # .. todo:: decide what goes into the effects table of meta data
    return fits.TableHDU()


# def get_detector_list(effects):
#     detector_lists = get_all_effects(effects, efs.DetectorList)
#
#     if len(detector_lists) != 1:
#         logging.warning("None or more than one DetectorList found. Using the"
#                       " first instance.{}".format(detector_lists))
#
#     return detector_lists[0]
