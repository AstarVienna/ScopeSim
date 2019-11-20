from astropy.io import fits
import anisocado as aniso


def make_standard_scao_constpsf():
    waves = [1.2, 1.6, 2.15]
    psfs = []
    for wave in waves:
        psf = aniso.AnalyticalScaoPsf(pixelSize=0.004, N=256, wavelength=wave)
        psf.shift_off_axis(0, 5)
        psfs += [psf.hdu]

    keys = {"AUTHOR" : "Kieran Leschinski",
            "DATE_CRE" : "2019-07-30",
            "DATE_MOD" : "2019-07-30",
            "SOURCE" : "AnisoCADO",
            "STATUS" : "Best guess for a standard observations",
            "ETYPE" : "CONSTPSF",
            "ECAT" : (-1, "Field constant. No layer catalogue"),
            "EDATA" : (1, "PSFs begin from EXT 1"),
            "XOFFSET": (0, "[arcsec] offset from NGS"),
            "YOFFSET": (5, "[arcsec] offset from NGS"),
            }

    pri = fits.PrimaryHDU()
    pri.header.update(keys)

    hdus = fits.HDUList([pri] + psfs)
    hdus.writeto("SCAO_ConstPSF_5off.fits", clobber=True)
    print(hdus.info())

