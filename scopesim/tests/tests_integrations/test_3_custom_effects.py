#!/usr/bin/env python
# coding: utf-8

# 3: Writing and including custom Effects
# =======================================
# 
# In this tutorial, we will load the model of MICADO (including Armazones, ELT, MORFEO) and then turn off all effect that modify the spatial extent of the stars. The purpose here is to see in detail what happens to the **distribution of the stars flux on a sub-pixel level** when we add a plug-in astrometric Effect to the optical system.
# 
# For real simulation, we will obviously leave all normal MICADO effects turned on, while still adding the plug-in Effect. Hopefully this tutorial will serve as a refernce for those who want to see **how to create Plug-ins** and how to manipulate the effects in the MICADO optical train model.
# 
# 
# Create and optical model for MICADO and the ELT
# -----------------------------------------------

def run():
    # In[1]:
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    import scopesim as sim
    from scopesim_templates.stellar import stars, star_grid

    # sim.rc.__config__["!SIM.file.local_packages_path"] = "../temp/"
    sim.rc.__config__["!SIM.file.local_packages_path"] = "F:/Work/irdb/"


    # We assume that the MICADO (plus support) packages have been downloaded.
    #
    # We can see which Effects are already included by calling ``micado.effects``:

    # In[2]:


    cmd = sim.UserCommands(use_instrument="MICADO", set_modes=["SCAO", "IMG_1.5mas"])
    micado = sim.OpticalTrain(cmd)

    micado.effects


    # Now we turn off all ``Effects`` that cause spatial aberrations:

    # In[3]:


    for effect_name in ["full_detector_array", "micado_adc_3D_shift",
                        "micado_ncpas_psf", "relay_psf"]:
        micado[effect_name].include = False
        print(micado[effect_name])


    # The normal detector window is set to 1024 pixels square.
    # Let's reduce the size of the detector readout window:

    # In[4]:


    micado["detector_window"].data["x_cen"] = 0.          # [mm] distance from optical axis on the focal plane
    micado["detector_window"].data["y_cen"] = 0.
    micado["detector_window"].data["xhw"] = 0.015 * 21.5  # [mm] x width [pixel]
    micado["detector_window"].data["yhw"] = 0.015 * 21.5


    # By default ScopeSim works on the whole pixel level for saving computation time.
    # However it is capable of integrating sub pixel shift.
    # For this we need to turn on the sub-pixel mode:

    # In[5]:


    micado.cmds["!SIM.sub_pixel.flag"] = True


    # We can test what's happening by making a grid of stars and observing them:

    # In[6]:


    src = star_grid(n=9, mmin=20, mmax=20.0001, separation=0.0015 * 15)
    src.fields[0]["x"] -= 0.00075
    src.fields[0]["y"] -= 0.00075

    micado.observe(src, update=True)

    plt.figure(figsize=(8,8))
    plt.imshow(micado.image_planes[0].data, origin="lower")

    plt.show()

    # In[18]:


    micado["detector_window"].data


    # Writing a custom Effect object
    # ------------------------------
    #
    # The following code snippet creates a new ``Effect`` class.
    #

    # In[ ]:


    import numpy as np
    from astropy.table import Table

    from scopesim.effects import Effect
    from scopesim.base_classes import SourceBase


    class PointSourceJitter(Effect):
        def __init__(self, **kwargs):
            super(PointSourceJitter, self).__init__(**kwargs)   # initialise the underlying Effect class object
            self.z_order = [500]                        # z_order number for knowing when and how to apply the Effect
            self.meta["max_jitter"] = 0.001                     # [arcsec] - a parameter needed by the effect
            self.meta.update(kwargs)                            # add any extra parameters passed when initialising

        def apply_to(self, obj):                                # the function that does the work
            if isinstance(obj, SourceBase):
                for field in obj.fields:
                    if isinstance(field, Table):
                        dx, dy = 2 * (np.random.random(size=(2, len(field))) - 0.5)
                        field["x"] += dx * self.meta["max_jitter"]
                        field["y"] += dy * self.meta["max_jitter"]

            return obj


    # Lets break it down a bit::
    #
    #     class PointSourceJitter(Effect):
    #
    # Here we are subclassing the ``Effect`` object from ScopeSim.
    # This has the basic functionality for reading in ASCII and FITS files, and for communicating with the ``OpticsManager`` class in ScopeSim.
    #
    # The initialisation function looks like this::
    #
    #     def __init__(self, **kwargs):
    #         super(PointSourceJitter, self).__init__(**kwargs)   # initialise the underlying Effect class object
    #         self.z_order = [500]
    #
    # Here we make sure to activate the underlying Effect object.
    # The ``z_order`` keyword in the meta dictionary is used by ScopeSim to determine when and where this Effect should be applied during a simulations run.
    # The exact z-order numbers are described in [insert link here].
    #
    # The main function of any Effect is the ``apply_to`` method::
    #
    #     def apply_to(self, obj):
    #         if isinstance(obj, SourceBase):
    #             ...
    #
    #         return obj
    #
    # It should be noted that what is passed in via (``obj``) must be returned in the same format. The contents of the ``obj`` can change, but the ``obj`` object must be returned.
    #
    # All the code which enacts the results of the physical effect are contained in this method.
    # For example, if we are writing a redshifting Effect, we could write the code to shift the wavelength array of a ``Source`` object by ``z+1`` here.
    #
    # There are 4 main classes that are cycled through during an observation run:
    # * ``SourceBase``: contains the original 2+1D distribtion of light,
    # * ``FieldOfViewBase``: contains a (quasi-)monochromatic cutout from the Source object,
    # * ``ImagePlaneBase``: contains the expectation flux image on the detector plane
    # * ``DetectorBase``: contains the electronic readout image
    #
    # An ``Effect`` object can be applied to any number of objects based on one or more of these base classes.
    # Just remember to segregate the base-class-specific code with ``if`` statements.
    #
    # One further method should be mentioned: ``def fov_grid()``.
    # This method is used by ``FOVManager`` to estimate how many ``FieldOfView`` objects to generate in order to best simulation the observation.
    # If your Effect object might alter this estimate, then you should include this method in your class. See the code base for further details.
    #
    # .. note:: The ``fov_grid`` method will be depreciated in a future release of ScopeSim.
    #     It will most likely be replaced by a ``FOVSetupBase`` class that will be cycled through the ``apply_to`` function.
    #     However this is not yet 100% certain, so please bear with us.
    #

    # Including a custom Effect
    # -------------------------
    #
    # First we need to initialise an instance of the Effect object:

    # In[ ]:


    jitter_effect = PointSourceJitter(max_jitter=0.001, name="random_jitter")


    # Then we can add it to the optical model:

    # In[ ]:


    micado.optics_manager.add_effect(jitter_effect)

    micado.effects


    # When we want to observe, we need to include the ``update=True`` flag so that the optical model is updated to include the instance of our new ``Effect``:
    #

    # In[ ]:


    micado.observe(src, update=True)

    plt.figure(figsize=(8,8))
    plt.imshow(micado.image_planes[0].data, origin="lower")
    plt.show()

    # We can update the parameters of the object on-the-fly by accessing the meta dictionary:

    # In[ ]:


    micado["random_jitter"].meta["max_jitter"] = 0.005

    micado.observe(src, update=True)

    plt.figure(figsize=(8,8))
    plt.imshow(micado.image_planes[0].data, origin="lower")


    # Here we can see that there is a certain amount of sub-pixel jitter being introduced into each observation.
    # However this bare-bones approach is not very realistic.
    # We should therefore turn the PSF back on to get a more realistic observation:

    # In[ ]:


    micado["relay_psf"].include = True

    micado.observe(src, update=True)
    hdus = micado.readout()

    plt.figure(figsize=(8,8))
    plt.imshow(hdus[0][1].data, origin="lower", norm=LogNorm())

    plt.show()


# run()
