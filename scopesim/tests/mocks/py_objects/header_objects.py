from astropy import wcs


def _basic_fov_header():
    w, h = 150, 150
    skywcs = wcs.WCS(naxis=2)
    skywcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    skywcs.wcs.cdelt = [0.1, 0.1]
    skywcs.wcs.cunit = ["arcsec", "arcsec"]
    skywcs.wcs.crval = [0, 0]
    skywcs.wcs.crpix = [w / 2, h / 2]

    detwcs = wcs.WCS(naxis=2, key="D")
    detwcs.wcs.ctype = ["LINEAR", "LINEAR"]
    detwcs.wcs.cdelt = [1, 1]
    detwcs.wcs.cunit = ["mm", "mm"]
    detwcs.wcs.crval = [0, 0]
    detwcs.wcs.crpix = [w / 2, h / 2]

    skyhdr = skywcs.to_header()
    dethdr = detwcs.to_header()
    skyhdr.update(dethdr)
    skyhdr["NAXIS"] = 2
    skyhdr["NAXIS1"] = w
    skyhdr["NAXIS2"] = h

    return skyhdr


def _implane_header():

    w, h = 150, 150
    detwcs = wcs.WCS(naxis=2, key="D")
    detwcs.wcs.ctype = ["LINEAR", "LINEAR"]
    detwcs.wcs.cdelt = [1, 1]
    detwcs.wcs.cunit = ["mm", "mm"]
    detwcs.wcs.crval = [0, 0]
    detwcs.wcs.crpix = [w / 2, h / 2]

    dethdr = detwcs.to_header()
    dethdr["NAXIS"] = 2
    dethdr["NAXIS1"] = w
    dethdr["NAXIS2"] = h

    return dethdr


def _fov_header():
    w, h = 100, 100
    skywcs = wcs.WCS(naxis=2)
    skywcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    skywcs.wcs.cdelt = [0.2, 0.2]
    skywcs.wcs.cunit = ["arcsec", "arcsec"]
    skywcs.wcs.crval = [0, 0]
    skywcs.wcs.crpix = [w / 2, h / 2]

    detwcs = wcs.WCS(naxis=2, key="D")
    detwcs.wcs.ctype = ["LINEAR", "LINEAR"]
    detwcs.wcs.cdelt = [1, 1]
    detwcs.wcs.cunit = ["mm", "mm"]
    detwcs.wcs.crval = [0, 0]
    detwcs.wcs.crpix = [w / 2, h / 2]

    skyhdr = skywcs.to_header()
    dethdr = detwcs.to_header()
    skyhdr.update(dethdr)
    skyhdr["NAXIS"] = 2
    skyhdr["NAXIS1"] = w
    skyhdr["NAXIS2"] = h

    return skyhdr
