import numpy as np
from matplotlib import pyplot as plt

import scopesim as sim
from scopesim.effects.psfs import GridFieldVaryingPSF
import synphot
from astropy.table import Table
from astropy import units as u

# sim.download_packages(["LFOA"])
# sim.download_packages(["Armazones", "ELT", "MICADO", "MAORY"])

# Create sims and load mds
cmd = sim.UserCommands(use_instrument="MICADO", set_modes=["SCAO", "IMG_1.5mas"])
micado = sim.OpticalTrain(cmd)
micado.cmds["!SIM.sub_pixel.flag"] = True

for effect_name in ["full_detector_array", "micado_adc_3D_shift",
                    "micado_ncpas_psf", "relay_psf"]:
    micado[effect_name].include = False
    print(micado[effect_name])

micado["detector_window"].data["x_cen"] = 0  # [mm] distance from optical axis on the focal plane
micado["detector_window"].data["y_cen"] = 0
micado["detector_window"].data["x_size"] = 256  # [pixel] width of detector
micado["detector_window"].data["y_size"] = 256

params = {
    "filename": "psf_grid.fits",
    "wave_key": "WAVELENG"}

psf = GridFieldVaryingPSF(name="psf_grid", **params)
micado.optics_manager.add_effect(psf)

# Object
pixel_scale = 0.0015 * u.arcsec / u.pixel
vega = synphot.spectrum.SourceSpectrum.from_vega()
sep = 64
i = [*([-sep, 0, sep] * 3)]
j = [*([sep] * 3), *([0] * 3), *([-sep] * 3)]
x = [ii * u.pixel * pixel_scale for ii in i]
y = [jj * u.pixel * pixel_scale for jj in j]
ref = [0] * len(x)
weight = [*([10 ** -1] * 3), *([10 ** -2] * 3), *([10 ** -3] * 3)]
obj = Table(names=["x", "y", "ref", "weight"],
            data=[x, y, ref, weight],
            units=[u.arcsec, u.arcsec, None, None])
src = sim.Source(table=obj, spectra=[vega])

# # Baseline
micado["psf_grid"].include = False
micado.observe(src, update=True)

plt.figure(figsize=(8, 8))
plt.imshow(micado.image_planes[0].data, origin="lower")
plt.title("Baseline")
plt.show()


# GridFieldVaryingPSF
micado["psf_grid"].include = True
print(micado.effects)

micado.observe(src, update=True)
plt.figure(figsize=(8, 8))
img = micado.image_planes[0].data
plt.imshow(img, origin="lower")
plt.title("Custom PSF")
plt.show()


s_img = 256
n = 1000
i = np.random.rand(n) * (s_img-1) - s_img // 2
j = np.random.rand(n) * (s_img-1) - s_img // 2
x = [ii * u.pixel * pixel_scale for ii in i]
y = [jj * u.pixel * pixel_scale for jj in j]
ref = [0] * len(x)
weight = [1e-2]*n
obj_rand = Table(names=["x", "y", "ref", "weight"],
                 data=[x, y, ref, weight],
                 units=[u.arcsec, u.arcsec, None, None])
src_rand = sim.Source(table=obj_rand, spectra=[vega])

# Baseline (random)
micado["psf_grid"].include = False
micado.observe(src_rand, update=True)
micado[effect_name].include = False

plt.figure(figsize=(8, 8))
plt.imshow(micado.image_planes[0].data, origin="lower")
plt.title("Baseline (random)")
plt.show()


# GridFieldVaryingPSF (random)
micado["psf_grid"].include = True
micado.observe(src_rand, update=True)
plt.figure(figsize=(8, 8))
img = micado.image_planes[0].data
plt.imshow(img, origin="lower")
plt.title("Custom PSF (random)")
plt.show()
