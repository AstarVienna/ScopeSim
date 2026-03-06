# -*- coding: utf-8 -*-
"""Spectral grating efficiencies."""

from typing import ClassVar, Callable

import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from astropy.table import Table

from .effects import Effect
from .ter_curves import TERCurve
from .ter_curves_utils import apply_throughput_to_cube
from ..utils import figure_factory, get_logger, check_keys
from .data_container import DataContainer
from ..optics import echelle


logger = get_logger(__name__)


class SpectralEfficiency(Effect):
    """
    Applies the grating efficiency (blaze function) for a SpectralTraceList.

    Input Data Format
    -----------------
    The efficiency curves are taken from a fits file `filename`with a
    structure similar to the trace definition file (see `SpectralTraceList`).
    The required extensions are:
    - 0 : PrimaryHDU [header]
    - 1 : BinTableHDU or TableHDU[header, data] : Overview table of all traces
    - 2..N : BinTableHDU or TableHDU : Efficiency curves, one per trace. The
             tables must have the two columns `wavelength` and `efficiency`

    Note that there must be one extension for each trace defined in the
    `SpectralTraceList`. Extensions for other traces are ignored.

    EXT 0 : PrimaryHDU
    ++++++++++++++++++
    Required header keywords:

    - ECAT : int : Extension number of overview table, normally 1
    - EDATA : int : Extension number of first Trace table, normally 2

    No data is required in this extension

    EXT 1 : (Bin)TableHDU : Overview of traces
    ++++++++++++++++++++++++++++++++++++++++++
    No special header keywords are required in this extension.

    Required Table columns:
    - description : str : identifier for each trace
    - extension_id : int : which extension is each trace in

    EXT 2 : (Bin)TableHDU : Efficiencies for individual traces
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Required header keywords:
    - EXTNAME : must be identical to the `description` in EXT 1

    Required Table columns:
    - wavelength : float : [um]
    - efficiency : float : number [0..1]

    """

    z_order: ClassVar[tuple[int, ...]] = (630,)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        if "hdulist" in kwargs and isinstance(kwargs["hdulist"], fits.HDUList):
            self._file = kwargs["hdulist"]

        self.efficiencies = self.get_efficiencies()

    def get_efficiencies(self):
        """Read effciencies from file, returns a dictionary."""
        hdul = self._file
        self.ext_data = hdul[0].header["EDATA"]
        self.ext_cat = hdul[0].header["ECAT"]
        self.catalog = Table(hdul[self.ext_cat].data)

        efficiencies = {}
        for row in self.catalog:
            params = {col: row[col] for col in row.colnames}
            params.update(self.meta)
            hdu = self._file[row["extension_id"]]
            name = hdu.header['EXTNAME']

            tbl = Table.read(hdu)
            wavelength = tbl['wavelength'].quantity
            efficiency = tbl['efficiency'].value
            params.pop("filename", None)  # don't pass filename to TERCurve!
            effic_curve = TERCurve(array_dict={"wavelength":wavelength,
                                   "transmission":efficiency},
                                   **params)
            efficiencies[name] = effic_curve

        hdul.close()
        return efficiencies

    def apply_to(self, obj, **kwargs):
        """Interface between FieldOfView and SpectralEfficiency."""
        trace_id = obj.trace_id
        try:
            effic = self.efficiencies[trace_id]
        except KeyError:
            logger.warning("No grating efficiency for trace %s", trace_id)
            return obj

        swcs = WCS(obj.hdu.header).spectral
        with u.set_enabled_equivalencies(u.spectral()):
            wave = swcs.pixel_to_world(np.arange(swcs.pixel_shape[0])) << u.um
        obj.hdu = apply_throughput_to_cube(obj.hdu, effic.throughput, wave)
        return obj

    def plot(self):
        """Plot the grating efficiencies."""
        fig, axes = figure_factory()
        for name, effic in self.efficiencies.items():
            wave = effic.throughput.waveset
            axes.plot(wave.to(u.um), effic.throughput(wave), label=name)

        axes.set_xlabel("Wavelength [um]")
        axes.set_ylabel("Grating efficiency")
        axes.set_title(f"Grating efficiencies {self.display_name}")
        axes.legend()

        return fig


class EchelleSpectralEfficiency(Effect):
    """
    Spectral efficiency list from analytical calculations of the blaze function for ZShooter gratings.
    Requires same input trace parameter table as EchelleSpectralTraceList, supply as kwarg "filename"
    """
    z_order: ClassVar[tuple[int, ...]] = (630,)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.efficiency_generator = self._generate_efficiency_curve_func(DataContainer(filename=kwargs.pop('filename')))
        self.efficiencies = {}

    def _generate_efficiency_curve_func(self, trace_params) -> Callable:
        spectrographs = {}
        for row in trace_params.table:
            prefix = row["prefix"]  # note trance ids are assumed to be prefix_{order}
            min_order = row['m0'] - row['n']
            max_order = row['m0']
            min_wave = row['min_wave'] * u.Unit(trace_params.meta["min_wave_unit"])
            max_wave = row['max_wave'] * u.Unit(trace_params.meta["max_wave_unit"])
            design_res = row['design_res']
            focal_len = row['focal_length'] * u.Unit(trace_params.meta["focal_length_unit"])
            disp_npix = row['n_disp'] - 2 * row['detector_pad']
            xdisp_npix = row['n_xdisp']- 2 * row['detector_pad']
            pix_size = row['pixel_size'] * u.Unit(trace_params.meta["pixel_size_unit"])
            echelle_angle = np.deg2rad(row['echelle_blaze'])*u.rad
            xdisp_beta_center = np.deg2rad(row['xbeta_center'])*u.rad

            xdisp_groove_length = u.Unit(trace_params.meta["xdisp_freq_unit"]) / row['xdisp_freq']
            echelle_groove_length = u.Unit(trace_params.meta["disp_freq_unit"]) / row['disp_freq']
            pix_per_res_elem = row['fwhm']

            spectrograph = echelle.spectrograph_factory(min_wave, max_wave, focal_len,
                                                        design_res, echelle_angle, min_order, max_order,
                                                        echelle_groove_length, pix_per_res_elem, disp_npix, xdisp_npix,
                                                        pix_size, xdisp_groove_length=xdisp_groove_length,
                                                        xdisp_beta_center=xdisp_beta_center)

            spectrographs[prefix] = spectrograph

        def efficiency_curve(trace_id, wavelength):
            """Trace ID MUST be in the form prefix_{order}"""
            prefix, _, order = trace_id.partition('_')
            order = int(order)
            spec = spectrograph[prefix]
            # blaze = spec.blaze(wavelength)[trace_id]  # compute all of them and then subscript
            blaze = spec.grating.blaze(spec.grating.beta(wavelength, order), order)
            xdisp = spec.xdisp_efficiency(wavelength)
            return blaze*xdisp

        return efficiency_curve

    def apply_to(self, obj, **kwargs):
        """Interface between FieldOfView and SpectralEfficiency."""
        trace_id = obj.trace_id

        swcs = WCS(obj.hdu.header).spectral
        with u.set_enabled_equivalencies(u.spectral()):
            wave = swcs.pixel_to_world(np.arange(swcs.pixel_shape[0])) << u.um

        try:
            efficiency = self.efficiency_generator(trace_id, wave)
            params = {"description": trace_id}
            params.update(self.meta)
            effic_curve = TERCurve(array_dict={"wavelength": wave, "transmission": efficiency}, **params)
            self.efficiencies[trace_id] = effic_curve
        except:
            raise ValueError(f"Error generating efficiency curve for trace {trace_id} with wavelength range {wave.min()} - {wave.max()}")

        obj.hdu = apply_throughput_to_cube(obj.hdu, effic_curve.throughput, wave)
        return obj

    def plot(self):
        """Plot the grating efficiencies."""
        fig, axes = figure_factory()
        for name, effic in self.efficiencies.items():
            wave = effic.throughput.waveset
            axes.plot(wave.to(u.um), effic.throughput(wave), label=name)

        axes.set_xlabel("Wavelength [um]")
        axes.set_ylabel("Grating efficiency")
        axes.set_title(f"Grating efficiencies {self.display_name}")
        axes.legend()

        return fig