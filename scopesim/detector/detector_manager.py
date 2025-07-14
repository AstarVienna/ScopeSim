# -*- coding: utf-8 -*-
"""Contains DetectorManager and aux functions."""

from collections.abc import Sequence

from astropy.io.fits import HDUList, PrimaryHDU, TableHDU

from .detector import Detector
from ..effects import Effect
from ..utils import stringify_dict, get_logger


logger = get_logger(__name__)


class DetectorManager(Sequence):
    """Manages the individual Detectors, mostly used for readout."""

    def __init__(self, detector_list=None, cmds=None, **kwargs):
        self.meta = {}
        self.meta.update(kwargs)
        self.cmds = cmds

        # The effect from which the instance is constructed
        self._detector_list = detector_list

        # TODO: Typing could get more specific once z_order is replaced with
        #       Effect subclasses, which might also be used to create the two
        #       lists more efficiently...
        self._array_effects: list[Effect] = []
        self._dtcr_effects: list[Effect] = []
        self._detectors: list[Detector] = []
        if self._detector_list is not None:
            self._detectors = [
                Detector(hdr, cmds=self.cmds, **self.meta)
                for hdr in self._detector_list.detector_headers()]
        else:
            logger.warning("No detector effect was passed, cannot fully "
                           "initialize detector manager.")

        self._latest_exposure: HDUList | None = None

    def readout(self, image_planes, array_effects=None, dtcr_effects=None,
                **kwargs) -> HDUList:
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
        #           DetectorManager is the old Detector

        self._array_effects = array_effects or []
        self._dtcr_effects = dtcr_effects or []
        self.meta.update(kwargs)

        # 1. Get the image plane that corresponds to this detector array
        # TODO: This silently only returns the first match, is that intended??
        image_plane = next(implane for implane in image_planes if
                           implane.id == self._detector_list.image_plane_id)

        # 2. Apply detector array effects (apply to the entire image plane)
        for effect in self._array_effects:
            image_plane = effect.apply_to(image_plane, **self.meta)

        # 3. iterate through all Detectors, extract image from image_plane
        logger.info("Extracting from %d detectors...", len(self))
        for detector in self:
            detector.extract_from(image_plane)

            # 5. apply all effects (to all Detectors)
            for effect in self._dtcr_effects:
                detector = effect.apply_to(detector)

            # 6. add necessary header keywords
            # .. todo: add keywords

        # FIXME: Why is this applied twice ???
        for effect in self._array_effects:
            image_plane = effect.apply_to(image_plane, **self.meta)

        self._latest_exposure = self._make_hdulist()

        return self._latest_exposure

    def latest_exposure(self) -> HDUList:
        """Return the (initially empty) HDUList produced by the readout."""
        if self._latest_exposure is None:
            logger.warning("Run readout before accessing .latest_exposure.")
        return self._latest_exposure

    def __getitem__(self, index) -> Detector:
        """x.__getitem__(y) <==> x[y]."""
        return self._detectors[index]

    def __len__(self) -> int:
        """Return len(self)."""
        return len(self._detectors)

    def __repr__(self):
        """Return repr(self)."""
        msg = (f"{self.__class__.__name__}"
               f"({self._detector_list!r}, **{self.meta!r})")
        return msg

    def __str__(self):
        """Return str(self)."""
        return f"{self.__class__.__name__} with {self._detector_list!s}"

    def _repr_pretty_(self, printer, cycle):
        """For ipython."""
        if cycle:
            printer.text(f"{self.__class__.__name__}(...)")
        else:
            printer.text(str(self))

    def _make_primary_hdu(self):
        """Create the primary header from meta data."""
        prihdu = PrimaryHDU()
        prihdu.header.update(stringify_dict(self.meta))
        return prihdu

    def _make_effects_hdu(self):
        # .. todo:: decide what goes into the effects table of meta data
        # effects = self._array_effects + self._dtcr_effects
        return TableHDU()

    def _make_hdulist(self):
        """Generate a HDUList with the ImageHDUs and any extras."""
        # TODO: effects_hdu unnecessary as long as make_effects_hdu does not do anything
        hdu_list = HDUList([
            self._make_primary_hdu(),
            *[dtcr.hdu for dtcr in self],
            # *self._make_effects_hdu(),
        ])
        return hdu_list
