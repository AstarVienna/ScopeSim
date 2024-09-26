# -*- coding: utf-8 -*-
"""
Contains ``SourceField`` and its subclasses.

While the ``Source`` object serves as the high-level interface between target
descriptions and the ScopeSim optical train, the actual information about the
observed objects is stored in the ``SourceField`` classes, which constitute the
members of ``Source.fields`` collection. Any target to be understood by
ScopeSim can be characterized by either a ``Table`` of point sources, a
two-dimensional image (``ImageHDU``) plus a separate (averaged) spectrum, or a
three-dimensional datacube containing spectral information for each spatial
pixel. This threefold abstraction is mirrored by the three final subclasses of
``SourceField``: ``TableSourceField`` for point source tables with a spectrum
reference for each individual point source, ``ImageSourceField`` for a 2D image
with an average spectrum, and finally ``CubeSourceField`` with a full 3D data
cube. The ``ImageSourceField`` and ``CubeSourceField`` also contain a ``WCS``
coordinate information and the wavelength axis of the ``CubeSourceField`` is
available via the ``CubeSourceField.wave`` attribute.

In previous versions of ScopeSim (pre-0.9), the ``Source.fields`` collection
simply held the individual ``Table`` and ``ImageHDU`` (2D or 3D) objects, which
are now stored in the ``.field`` attribute of each source field. This new
distinction of the different cases allows much clearer separation of the logic
required to handle various operations on those objects, such as plotting and
shifting the source, which previously had to incorporate a number of case
differentiations that made the ``Source`` class rather overloaded with logic.
This now also allows for well-structured validation logic of the individual
source field data upon creation of each ``SourceField`` subclass instance.

Creation of the source field classes is usually handled by the ``Source`` class
itself via its various constructions methods, so the user rarely interacts with
these classes directly, except for debugging. They serve more as an internal
abstraction layer to handle the different cases of target object descriptions,
as described above.

The following class diagram illustrates the relationship between the
``SourceField`` subclasses:

```mmd
classDiagram
class SourceField{+field}
class SpectrumSourceField{+spectra}
class HDUSourceField{+wcs}

SourceField <|-- SpectrumSourceField
SourceField <|-- HDUSourceField
SpectrumSourceField <|-- TableSourceField
SpectrumSourceField <|-- ImageSourceField
HDUSourceField <|-- ImageSourceField
HDUSourceField <|-- CubeSourceField
```

"""

from copy import deepcopy
from pathlib import Path
from typing import TextIO, Any
from dataclasses import dataclass, KW_ONLY, field as dataclass_field
# rename it dataclass_field to avoid confusion with source field

import numpy as np

from astropy.table import Table, Column
from astropy.io.registry import IORegistryError
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS, SingularMatrixError, FITSFixedWarning

from synphot import SourceSpectrum


from ..optics import image_plane_utils as imp_utils
from ..utils import (quantify, quantity_from_table, close_loop, get_logger,
                     convert_table_comments_to_dict)


logger = get_logger(__name__)


# TODO: consider making this a metaclass
@dataclass
class SourceField:
    """Base class for source fields, not meant to be instantiated."""

    field: Any
    _: KW_ONLY
    meta: dict = dataclass_field(default_factory=dict, repr=False)

    def _write_stream(self, stream: TextIO) -> None:
        raise NotImplementedError("Subclasses should implement this.")

    def __getitem__(self, key):
        # For backwards-combatibility to allow direct access of
        # Source.fields[x][y] if possible. Maybe in the long run get rid of
        # this and force the use of .field...
        # warn("Direct item access for source fields may become deprecated "
        #      "in the future. Use the .field attribute instead.",
        #      PendingDeprecationWarning, stacklevel=2)
        return self.field.__getitem__(key)

    def __setitem__(self, key, value):
        # For backwards-combatibility to allow direct access of
        # Source.fields[x][y] if possible. Maybe in the long run get rid of
        # this and force the use of .field...
        # warn("Direct item assignment for source fields may become deprecated "
        #      "in the future. Use the .field attribute instead.",
        #      PendingDeprecationWarning, stacklevel=2)
        self.field.__setitem__(key, value)

    @property
    def name(self) -> str:
        """Name of the object (if set)."""
        return self.meta.get("object", "<unknown>")


@dataclass
class SpectrumSourceField(SourceField):
    """Base class for source fields with separate spectra (no cube)."""

    spectra: dict


@dataclass
class TableSourceField(SpectrumSourceField):
    """Source field with table of point source(s)."""

    field: Table

    @classmethod
    def from_file(cls, filename: Path | str,
                  spectra: dict[int, SourceSpectrum],
                  **kwargs):
        """Load source table from file."""
        try:
            tbl = Table.read(filename)
            # There used to be a header combining functionality here...
            # hdr1 = fits.getheader(filename, 1)
            # hdr.update(hdr1)
            # tbl = Table(data, meta=dict(hdr))
            # tbl.meta.update(convert_table_comments_to_dict(tbl))
        except IORegistryError:
            logger.debug("Table format guessing failed, retry with ascii.")
            tbl = Table.read(filename, format="ascii")

        tbl.meta.update(convert_table_comments_to_dict(tbl))
        return cls(tbl, spectra=spectra, meta=kwargs)

    @classmethod
    def from_arrays(cls, x, y, ref, weight,
                    spectra: dict[int, SourceSpectrum],
                    **kwargs):
        """Construct source table from arrays for each column."""
        if weight is None:
            weight = np.ones(len(x))

        x = quantify(x, u.arcsec)
        y = quantify(y, u.arcsec)
        tbl = Table(names=["x", "y", "ref", "weight"],
                    data=[x, y, ref, weight])
        tbl.meta["x_unit"] = "arcsec"
        tbl.meta["y_unit"] = "arcsec"
        return cls(tbl, spectra=spectra, meta=kwargs)

    def __post_init__(self):
        """Validate input."""
        assert self.spectra, "Spectra must be non-empty for table source."
        if "weight" not in self.field.colnames:
            self.field.add_column(
                Column(name="weight", data=np.ones(len(self.field)))
            )
        self.meta.update(self.field.meta)

    def __len__(self) -> int:
        """Return len(self)."""
        return len(self.field)

    def _write_stream(self, stream: TextIO) -> None:
        stream.write(f"Table with {len(self)} rows, referencing "
                     f"spectra {set(self.spectra)}")

    def plot(self, axes, color) -> None:
        """Plot source."""
        axes.plot(self.field["x"], self.field["y"], color+".", label=self.name)

    def shift(self, dx, dy) -> None:
        """Shift source by dx, dy."""
        x = quantity_from_table("x", self.field, u.arcsec)
        x += quantify(dx, u.arcsec)
        self.field["x"] = x

        y = quantity_from_table("y", self.field, u.arcsec)
        y += quantify(dy, u.arcsec)
        self.field["y"] = y


@dataclass
class HDUSourceField(SourceField):
    """Base class for source fields with HDU."""

    field: fits.ImageHDU
    wcs: WCS | None = dataclass_field(default=None, init=False)

    def __new__(cls, *args, **kwargs):
        """Override creation to create subclasses."""
        if issubclass(cls, (CubeSourceField, ImageSourceField)):
            # Allow for direct subclass access
            return super().__new__(cls)

        field = kwargs.get("field", args[0])
        if field.header["NAXIS"] == 3:
            return super().__new__(CubeSourceField)
        if field.header["NAXIS"] == 2:
            return super().__new__(ImageSourceField)

        # If we get here, something went wrong
        raise TypeError(f"{field.header['NAXIS'] = } must be 2 or 3.")

    @property
    def header(self) -> fits.Header:
        """Shortcut for `field.header`."""
        return self.field.header

    @property
    def data(self) -> np.ndarray:
        """Shortcut for `field.data`."""
        return self.field.data

    @data.setter
    def data(self, value):
        self.field.data = value

    @property
    def img_size(self) -> str:
        """Shortcut for `field.data.shape`."""
        if self.data is None:
            return "<empty>"
        return str(self.data.shape)

    def _write_stream(self, stream: TextIO) -> None:
        stream.write(f"ImageHDU with size {self.img_size}, referencing "
                     f"spectrum {self.field.header.get('SPEC_REF', '-')}")

    def plot(self, axes, color) -> None:
        """Plot source."""
        xypts = imp_utils.calc_footprint(self.header)
        convf = u.Unit(self.header["CUNIT1"]).to(u.arcsec)
        outline = np.array(list(close_loop(xypts))) * convf
        axes.plot(*outline.T, color, label=self.name)

    def shift(self, dx, dy) -> None:
        """Shift source by dx, dy."""
        dx = dx << u.arcsec << self.wcs.wcs.cunit[0]
        dy = dy << u.arcsec << self.wcs.wcs.cunit[1]
        self.header["CRVAL1"] += dx.value
        self.header["CRVAL2"] += dy.value
        # TODO: self.wcs should be updated here!


@dataclass
class ImageSourceField(SpectrumSourceField, HDUSourceField):
    """Source field with 2D image and a single (average) spectrum."""

    def __post_init__(self):
        """Validate input."""
        assert self.spectra, "Spectra must be non-empty for 2D image source."
        try:
            self.wcs = WCS(self.field)
        except (SingularMatrixError, FITSFixedWarning):
            # This occurs for BG SRC
            logger.debug("Couldn't create source field WCS.")
            self.wcs = None


@dataclass
class CubeSourceField(HDUSourceField):
    """Source field with 3D data cube."""

    from_hdul: bool = False

    def __post_init__(self):
        """Validate input."""
        if self.wcs is None and not self.from_hdul:
            self.wcs = WCS(self.field)

        try:
            bunit = u.Unit(self.header["BUNIT"])
        except KeyError:
            bunit = u.erg / u.s / u.cm**2 / u.arcsec**2
            logger.warning(
                "Keyword \"BUNIT\" not found, setting to %s by default", bunit)
        except ValueError as error:
            logger.error("\"BUNIT\" keyword is malformed: %s", error)
            raise
        self.field.header["BUNIT"] = str(bunit)

    @classmethod
    def from_hdulist(cls, hdulist: fits.HDUList, ext: int = 0, **kwargs):
        """Load source cube from HDUL."""
        cube = fits.ImageHDU(header=hdulist[ext].header.copy(),
                             data=deepcopy(hdulist[ext].data))
        new_csf = cls(field=cube, meta=kwargs, from_hdul=True)
        new_csf.wcs = WCS(hdulist[ext], fobj=hdulist)
        return new_csf

    def shift(self, dx, dy) -> None:
        """Shift source by dx, dy."""
        logger.warning(
            "Source shift for cubes assumes first two axes are celestial.")
        super().shift(dx, dy)

    @property
    def wave(self) -> u.Quantity:
        """Construct wavelength axis for cube in um."""
        swcs = self.wcs.spectral
        with u.set_enabled_equivalencies(u.spectral()):
            wave = swcs.pixel_to_world(np.arange(swcs.pixel_shape[0])) << u.um
        return wave
