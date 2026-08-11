# Sources

## What packages do I need to create sources?

For common astronomical sources use
[ScopeSim Templates](https://scopesim-templates.readthedocs.io/en/latest/):

```bash
pip install scopesim_templates
```

For spectral libraries use
[SpeXtra](https://spextra.readthedocs.io/en/latest/) (many catalogues) or
[Pyckles](https://pyckles.readthedocs.io/en/latest/) (Pickles stellar library).

For hand-crafted sources the `scopesim.Source` class accepts arrays and FITS
images directly.

---

## How do I create a single star?

```python
import scopesim_templates as st

src = st.stellar.star(mag=20, filter_name="Ks", spec_type="A0V")
```

---

## How do I create a star field?

```python
src = st.stellar.star_field(
    n=200,
    mmin=18,
    mmax=24,
    width=60,       # [arcsec] square field side length
    filter_name="Ks",
)
```

---

## How do I create a stellar cluster?

```python
src = st.stellar.cluster(
    mass=1e4,        # [M_sun] total stellar mass
    distance=50000,  # [pc]
    filter_name="V",
    core_radius=0.3, # [pc]
)
```

---

## How do I create a galaxy?

```python
src = st.extragalactic.elliptical(
    total_magnitude=18,
    filter_name="Ks",
    pixel_scale=0.004,  # [arcsec/pix] output pixel scale
    r_eff=0.5,          # [arcsec] effective radius
    n=4,                # Sérsic index
    ellip=0.3,
    theta=45,
)
```

---

## How do I combine multiple sources?

Use the `+` operator:

```python
stars = st.stellar.star_field(100, "V", 18, 24, width=30)
galaxy = st.extragalactic.elliptical(20, "Ks", ...)

combined = stars + galaxy
sim(combined)
```

---

## How do I create a source from a FITS image?

```python
from astropy.io import fits
from scopesim import Source
import spextra as sp

# Load your image
hdu = fits.open("my_image.fits")[0]

# Need a spectrum — use SpeXtra or synphot
spectrum = sp.Spextrum("pickles/a0v")  # example: A0V star spectrum

src = Source(image_hdu=hdu, spectra=[spectrum])
```

The FITS image header must contain WCS keywords (``CDELT``, ``CRPIX``,
``CRVAL``) defining the pixel scale in arcsec/pix. Pixel values are treated
as flux weights multiplied by the given spectrum.

---

## How do I set the position of a source?

Source coordinates are in arcsec relative to the field centre:

```python
src = st.stellar.star(mag=20, filter_name="Ks")
src.shift(dx=2.5, dy=-1.0)   # offset in arcsec
```

For point-source arrays, positions are set via the `x` and `y` arguments of
most template functions.

---

## How do I check what my source looks like?

```python
src.plot()                  # spatial footprint
src.spectra[0].plot()       # first spectrum
```

Or view the source table:

```python
src.source_table            # x, y positions and spectral references
```
