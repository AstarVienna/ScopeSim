from astropy.io import fits
import anisocado as aniso


def make_standard_scao_constpsf():
    waves = [1.2, 1.6, 2.15]
    nmRms = [310, 310, 310]     # achieves the 40-18-6 Strehls for MAORY
    offset = [0, 0]
    psfs = []
    for wave, rms in zip(waves, nmRms):
        psf = aniso.AnalyticalScaoPsf(pixelSize=0.004, N=256, wavelength=wave,
                                      nmRms=rms)
        psf.shift_off_axis(offset[0], offset[1])
        psfs += [psf.hdu]

    keys = {"AUTHOR" : "Kieran Leschinski",
            "DATE_CRE" : "2019-07-30",
            "DATE_MOD" : "2019-07-30",
            "SOURCE" : "AnisoCADO",
            "STATUS" : "Best guess for a MAORY ConstantPSF with AnisoCADO",
            "ETYPE" : "CONSTPSF",
            "ECAT" : (-1, "Field constant. No layer catalogue"),
            "EDATA" : (1, "PSFs begin from EXT 1"),
            "XOFFSET": (offset[0], "[arcsec] offset from NGS"),
            "YOFFSET": (offset[1], "[arcsec] offset from NGS"),
            }

    pri = fits.PrimaryHDU()
    pri.header.update(keys)

    hdus = fits.HDUList([pri] + psfs)
    hdus.writeto("MCAO_ConstPSF_40_18_6.fits", clobber=True)
    for i in range(1, 4):
        print(hdus[i].header["WAVE0"], hdus[i].header["STREHL"])

    print(hdus.info())
