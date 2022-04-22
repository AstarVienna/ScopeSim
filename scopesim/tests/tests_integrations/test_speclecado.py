import pytest
import os
from astropy import units as u
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from scopesim import UserCommands, OpticalTrain
from scopesim.tests.mocks.py_objects import source_objects as so
from scopesim import effects as efs

from scopesim import rc
MOCK_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                        "../mocks/MICADO_SPEC/"))
rc.__search_path__.insert(0, MOCK_DIR)


PLOTS = False

@pytest.mark.skip(reason="Ignoring old Spectroscopy integration tests")
class TestMicadoSpec:
    def test_full_run_through(self):
        src1 = so._vega_source(mag=20)
        src2 = so._vega_source(mag=22)
        src3 = so._vega_source(mag=25)
        src2.fields[0]["x"][0] += 1
        src3.fields[0]["x"][0] -= 1
        src = src1 + src2 + src3

        cmd = UserCommands(yamls=["mock_MICADO_SPEC.yaml"])
        cmd["!SIM.file.local_packages_path"] = "./"
        opt = OpticalTrain(cmd)
        assert isinstance(opt, OpticalTrain)

        #opt.observe(src)
        # opt.readout(filename="temp_speclecado.fits")
        # opt.image_planes[0].hdu.writeto("temp_implane.fits", overwrite=True)

        if PLOTS:
            plt.imshow(opt.image_planes[0].data, origin="lower", norm=LogNorm())
            plt.show()


# def test_plot_spec_trace_layout():
#
#     spt = efs.SpectralTraceList(filename="TRACE_15arcsec.fits",
#                                 wave_colname="lam", s_colname="xi")
#     det = efs.DetectorList(filename="FPA_array_layout.dat")
#
#     if PLOTS:
#         spt.plot(1.4, 2.5)
#         det.plot()
#         plt.show()
#
#
# def test_reads_psf_scao_file():
#     psf = efs.FieldConstantPSF(filename="PSF_SCAO.fits")
#     assert isinstance(psf, efs.PSF)
#
#     if PLOTS:
#         print(psf._file.info())
#         plt.imshow(psf.get_data(1), norm=LogNorm())
#         plt.show()
#
#
# def test_spanish_vo_filter_timing():
#     filt = efs.SpanishVOFilterCurve(observatory="Paranal", instrument="HAWKI",
#                                     filter_name="H")
#     assert isinstance(filt, efs.SpanishVOFilterCurve)
