from ..effects.data_container import DataContainer
from ..base_classes import SourceBase, FieldOfViewBase, ImagePlaneBase, \
    DetectorBase
from .. import rc


class Effect(DataContainer):
    """
    The base class for representing the effects (artifacts) in an optical system

    The ``Effect`` class is conceived to independently apply the changes that
    an optical component (or series thereof) has on an incoming 3D description
    of an on-sky object. In other words, **an Effect object should receive a
    derivative of a ``Source`` object, alter it somehow, and return it**.

    The interface for the Effect base-class has been kept very general so that
    it can easily be sub-classed as data for new effects becomes available.
    Essentially, a sub-classed Effects object must only contain the following
    attributes:

    * ``self.meta`` - a dictionary to contain meta data.
    * ``self.apply_to(obj, **kwargs)`` - a method which accepts a
      Source-derivative and returns an instance of the same class as ``obj``
    * ``self.fov_grid(which="", **kwargs)``


    Parameters
    ----------
    See :class:`DataContainer` for input parameters

    Methods


    """

    def __init__(self, **kwargs):
        super(Effect, self).__init__(**kwargs)
        # .. todo:: should this be here, or only in apply_to() and fov_grid()?
        self.update_bang_keywords()
        self.meta["z_order"] = []

    def apply_to(self, obj):
        if not isinstance(obj, (SourceBase, FieldOfViewBase,
                                ImagePlaneBase, DetectorBase)):
            raise ValueError("object must one of the following: "
                             "Source, FieldOfView, ImagePlane, Detector: "
                             "{}".format(type(obj)))

        return obj

    def fov_grid(self, which="", **kwargs):
        """
        Returns the edges needed to generate FieldOfViews for an observation

        Parameters
        ----------
        which : str
            ["waveset", "edges", "shifts"] where:
            * waveset - wavelength bin extremes
            * edges - on sky coordinate edges for each FOV box
            * shifts - wavelength dependent FOV position offsets

        kwargs
        ------
        waverange : list
            [um] list of wavelength

        wave_mid
            [um] wavelength what will be centred on optical axis


        Returns
        -------
        waveset : list
            [um] N+1 wavelengths that set edges of N spectral bins

        edges : list of lists
            [arcsec] Contains a list of footprint lists

        shifts : list of 3 lists
            [wave, dx, dy] Contains lists corresponding to the (dx, dy) offset
            from the optical axis (0, 0) induced for each wavelength in (wave)
            [um, arcsec, arcsec]

        """
        self.update(**kwargs)
        return []

    def update(self, **kwargs):
        self.meta.update(kwargs)
        self.update_bang_keywords()

    def update_bang_keywords(self):
        for key in self.meta:
            if isinstance(self.meta[key], str) and self.meta[key][0] == "!":
                bang_key = self.meta[key]
                self.meta[key] = rc.__currsys__[bang_key]

    def __repr__(self):
        if "name" not in self.meta:
            self.meta["name"] = "<unknown name>"
        return '{}: "{}"'.format(type(self).__name__, self.meta["name"])


