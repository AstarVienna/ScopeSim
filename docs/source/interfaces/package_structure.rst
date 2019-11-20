Instrument Package structure
============================

Management of properties associated with optical elements
---------------------------------------------------------
The major idea here is to better manage the optical element specific information
(like ``OBS_DIT``, or ``ATMO_TEMPERATURE``). Instead of having an
``observation_dict`` where everything is in a single KEYWORD : VALUE pair format,
I've introduced the concept of the SystemDict. This is a subclassed ``dict``
which allows nested parameters (as delivered by yaml) to be referenced by a
so-called bang-string. Thus to get the value of
``dict["SIM"]["computing"]["chunk_size]``, it's now possible to just call
``dict["!SIM.computing.chunk_size"]``. A bang-sting must start with a ``!``
though - hence the name "bang"-string (see definition of a shebang #!)

All optical elements are now described by yaml files, and each optical element
**should have an ``alias`` tag**. Essentially the ``alias`` tag could be
anything, but for the sake of short-hand, I've stuck with the following for
MICADO: ``ATMO``, ``TEL``, ``RO`` (relay optics), ``INST``, ``DET``, ``OBS``,
and ``SIM``.

The ``UserCommands object`` now contains a ``SystemDict`` (in ``.cmds``) and all
yaml files are passed through it. Only values that are stored in the yaml
``properties`` dictionary are kept in the ``SystemDict`` and are referencable
by bang-string, i.e ``INST.filter_name``

Furthermore, whenever a ``UserCommand`` object is passed to an ``OpticalTrain``,
the ``OpticalTrain`` updates the global variable ``scopesim.rc.__currsys__`` to
reference the ``SystemDict`` inside the ``UserCommands``. This means that the
``properties`` values of any optical element can be accessed by any other
element (e.g. ``Effect``-objects) by simply calling
``rc.__currsys__["!<ALIAS>.<property_name>]``. This makes the management of
these properties, which are often needed by objects which aren't in the same
optical element (e.g. ADC needs atmospheric properties) much simpler.


Package format
--------------
Not directly dealt with here explicitly, but worth noting as the new
``UserCommands`` object has been designed for this. An instrument package should
look like this:

**Necessary files**::

    Top-level folder <package-name>
    |
    | - default.yaml
    | - <package-name>.yaml
    | - <package-name>_<mode-name>.yaml
    | - .... misc_support_data.files


**Optional files**::

    | - <package-name>_<extra-optics>.yaml
    | - <package-name>_<detector-name>.yaml


``default.yaml``
++++++++++++++++

The default yaml file contains a list of everything needed for an
``OpticalTrain`` to be successfully built. This means it needs to reference any
other packages which are needed for the full optical system. In the case of
MICADO this means including the the location (Armazones), the telescope (ELT),
and possibly an external relay optics (MAORY). These packages willl need to be
downloaded separately but should still be references in the ``default.yaml``
file with the ``packages`` list.

The default file should also contain a list of all the yaml files needed for the
default configuration of the instrument, including those yaml files which are
contained in other packages. This is only necessary if the package
is a "primary package" - i.e. one that will be used to create detector readout
images, e.g. MICADO. If the package is a "support package" (e.g. ELT, MAORY),
then this list isn't needed. **This list is referenced using the ``yamls``
keyword**

The last keyword in the default yaml file should be the ``properties`` keyword.
Here we put any global values relating to the observation such as
``filter_name``, or ``dit`` and ``ndit``, or ``airmass``.

Here is an example ot the ``OBS`` optical element yaml dictionary in a
``default.yaml`` file::

    ### default observation parameters for a MICADO simulation
    name : micado_default_params
    alias : OBS

    packages:
    - Armazones
    - ELT
    - MAORY
    - MICADO

    yamls:
    - Armazones.yaml
    - ELT.yaml
    - MAORY.yaml
    - MICADO.yaml
    - MICADO_IMG_wide.yaml
    - MICADO_H4RG.yaml

    properties:
        filter_name : "Ks"
        airmass : 1.2
        dit : 20
        ndit : 100


Generally it shouldn't be needed, but if there are any ``scopesim`` simulation
parameters which need to be specially configured for the package, these can
be added in a second yaml dictionary below the ``OBS`` dictionary.

.. Note::
    Yaml dictionaries must be separated by three dashes ``---``. This signals to
    ``pyyaml`` to treat each as an individual nested dictionary object.

An example of an additional ``SIM`` dictionary is::

    ---

    ### MICADO-specific simulation parameter
    alias: SIM

    properties :
        computing :
            chunk_size : 512

        spectral :
            wave_min : 0.7
            wave_mid : 1.2
            wave_max : 2.5

All default parameters can be found in the ``simcado.rc.__config__`` dictionary.

<package-name>.yaml vs <package-name>_<mode-name>.yaml
++++++++++++++++++++++++++++++++++++++++++++++++++++++

The file ``<package_name>.yaml`` contains the list of ``Effect`` objects and
their default parameters for the optics which are always static. E.g the
entrance window transmission curve, the number of static mirrors, etc

Any optics which can be moved into or out of the optical path and belong to a
specific mode configuration should be described in a separate yaml file.
For MICADO these include the removable optics mirrors or grating, and the
spectral order trace files. Here is where properties like ``pixel_scale`` should
be kept, as this is a property of a specific mode configuration.

Settings like the filter or slit choice should be kept in the main
``<package-name>.yaml`` in the ``properties`` section with a bang-string
referencing a dynamic value in the ``OBS`` dictionary. This way the value of the
variable can be changed without having to dig deeply into the description of the
instrument.

An example of a ``<package_name>.yaml`` file, note the ``filter_name`` property::

    ### MICADO INSTRUMENT WIDE FIELD MODE
    object : instrument
    alias : INST
    name : MICADO
    description : base configuration for MICADO

    properties :
        temperature : -190

    effects :
    -   name: micado_static_surfaces
        description : surfaces list for wide field optics
        class: SurfaceList
        kwargs:
            filename: LIST_MICADO_mirrors_static.dat

    -   name: micado_filter
        description : transmission curce for filter
        class: TERCurve
        kwargs:
            filename: "!OBS.filter_name"


Quick note on exponent notation floats in yaml files
----------------------------------------------------
.. note:: Preface floats like `1e6` in yaml files with the ``!!float`` keyword

    ``pyyaml`` is generally pretty good at recognising variable types.
    However, if you want to specify a value in exponent notation, i.e. ``4.2e1``
    instead of ``42``, ``pyyaml`` will assume that you have written a string,
    not a number. To override this, make sure to preface the number with the
    ``pyyaml`` keyword ``!!float``. E.g.

    ``power_in_watts: !!float 1.21e9``
