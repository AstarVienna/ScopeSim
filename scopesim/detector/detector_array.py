# -*- coding: utf-8 -*-
"""Contains DetectorArray and aux functions."""

from astropy.io import fits

from .detector import Detector
from ..utils import stringify_dict, get_logger


logger = get_logger(__name__)


class DetectorArray:
    """Manages the individual Detectors, mostly used for readout."""

    # TODO: Should this class be called DetectorManager analogous to the
    #       FOVManager and OpticsManager classes?
    # TODO: Either way this could inherit from some collections.abc

    def __init__(self, detector_list=None, cmds=None, **kwargs):
        self.meta = {}
        self.meta.update(kwargs)
        self.cmds = cmds

        # The effect from which the instance is constructed
        self._detector_list = detector_list

        self.array_effects = []
        self.dtcr_effects = []
        self.detectors = []

        # The (initially empty) HDUList produced by the readout
        self.latest_exposure = None

    def readout(self, image_planes, array_effects=None, dtcr_effects=None,
                **kwargs) -> fits.HDUList:
        """
        Read out the detector array into a FITS HDU List.

        1. Select the relevant image plane to extract images from.
        2. Apply detector array effects (apply to the entire image plane)
        3. Make a series of Detectors for each row in a DetectorList object.
        4. Iterate through all Detectors, extract image from image_plane.
        5. Apply all effects (to all Detectors).
        6. Add necessary header keywords (not implemented).
        7. Generate a HDUList with the ImageHDUs and any extras:
          - add ``PrimaryHDU`` with meta data regarding observation in header
          - add ``ImageHDU`` objects
          - add ``ASCIITableHDU`` with Effects meta data in final table
            extension (not implemented)

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
        latest_exposure : fits.HDUList
            Output FITS HDU List.

        """
        # .. note:: Detector is what used to be called Chip
        #           DetectorArray is the old Detector

        self.array_effects = array_effects or []
        self.dtcr_effects = dtcr_effects or []
        self.meta.update(kwargs)

        # 1. Get the image plane that corresponds to this detector array
        # TODO: This silently only returns the first match, is that intended??
        image_plane = next(implane for implane in image_planes if
                           implane.id == self._detector_list.image_plane_id)

        # 2. Apply detector array effects (apply to the entire image plane)
        for effect in self.array_effects:
            image_plane = effect.apply_to(image_plane, **self.meta)

        # 3. make a series of Detectors for each row in a DetectorList object
        self.detectors = [Detector(hdr, cmds=self.cmds, **self.meta)
                          for hdr in self._detector_list.detector_headers()]

        # 4. iterate through all Detectors, extract image from image_plane
        logger.info("Extracting from %d detectors...", len(self.detectors))
        for detector in self.detectors:
            detector.extract_from(image_plane)

            # 5. apply all effects (to all Detectors)
            for effect in self.dtcr_effects:
                detector = effect.apply_to(detector)

            # 6. add necessary header keywords
            # .. todo: add keywords

        # 7. Generate a HDUList with the ImageHDUs and any extras:
        pri_hdu = make_primary_hdu(self.meta)

        # ..todo: effects_hdu unnecessary as long as make_effects_hdu does not do anything
        # effects_hdu = make_effects_hdu(self.array_effects + self.dtcr_effects)

        hdu_list = fits.HDUList([pri_hdu]
                                + [dtcr.hdu for dtcr in self.detectors])
                                # + [effects_hdu])

        for effect in self.array_effects:
            image_plane = effect.apply_to(image_plane, **self.meta)

        self.latest_exposure = hdu_list

        return self.latest_exposure

    def __repr__(self):
        msg = (f"{self.__class__.__name__}"
               f"({self._detector_list!r}, **{self.meta!r})")
        return msg

    def __str__(self):
        return f"{self.__class__.__name__} with {self._detector_list!s}"

    def _repr_pretty_(self, p, cycle):
        """For ipython."""
        if cycle:
            p.text(f"{self.__class__.__name__}(...)")
        else:
            p.text(str(self))


def make_primary_hdu(meta):
    """Create the primary header from meta data."""
    prihdu = fits.PrimaryHDU()
    prihdu.header.update(stringify_dict(meta))
    return prihdu


def make_effects_hdu(effects):
    # .. todo:: decide what goes into the effects table of meta data
    return fits.TableHDU()
