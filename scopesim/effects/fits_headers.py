import yaml
import numpy as np
from astropy.io import fits
from . import Effect
from ..utils import check_keys


class ExtraFitsKeywords(Effect):
    """
    Extra FITS header keywords to be added to the pipeline FITS files

    These keywords are ONLY for keywords that should be MANUALLY ADDED to the
    headers after a simulation is read-out by the detector.

    Simulation parameters (Effect kwargs values, etc) will be added automatically
    by ScopeSim in a different function, but following this format.

    The dictionaries should be split into different HIERARCH lists:

    - HIERARCH ESO
     For ESO specific keywords
    - HIERARCH SIM
     For ScopeSim specific keywords, like simulation parameters
    - HIERARCH MIC
     For MICADO specific keywords, (unsure what these would be yet)

    More HIERARCH style keywords can also be added as needed for other use-cases

    Parameters
    ----------
    filename : str, optional
        Name of a .yaml nested dictionary file. See below for examples

    header_dict : nested dicts
        A series of nested python dictionaries following the format of the
        examples below. This keyword allows these dicts to be definied directly
        in the Effect yaml file, rather than in a seperate header keywords file.


    Yaml file format
    ----------------
    This document is a yaml document.
    Hence all new keywords should be specified in the form of yaml nested
    dictionaries.
    As each ``astropy.HDUList`` contains one or more extensions, the inital level is
    reserved for a list of keyword groups.
    For example::

    - ext_type: PrimaryHDU
      keywords:
        HIERARCH:
          ESO:
            ATM:
              TEMPERAT: -5

    - ext_number: [1, 2]
      keywords:
        HIERARCH:
          ESO:
            DET:
              DIT: [5, '[s] exposure length']   # example of adding a comment

    The keywords can be added to one or more extensions, based on one of the
    following ``ext_`` qualifiers: ``ext_name``, ``ext_number``, ``ext_type``

    Each of these ``ext_`` qualifiers can be a ``str`` or a ``list``.
    For a list, ScopeSim will add the keywords to all extensions matching the
    specified type/name/number

    The above example will result in the following keyword added to:

    - PrimaryHDU (ext 0):

      header['HIERARCH ESO ATM TEMPERAT'] = -5

    - Extensions 1 and 2 (regardless of type):

      header['HIERARCH ESO DET DIT'] = (5, '[s] exposure length')

    Resolved and un-resolved keywords
    ---------------------------------
    ScopeSim uses bang-strings to resolve global parameters.
    E.g: ``from_currsys('!ATMO.temperature')`` will resolve to a float
    These bang-strings will be resolved automatically in the ``keywords`` dictionary
    section.

    If the keywords bang-string should instead remain unresolved and the string
    added verbatim to the header, we use the ``unresolved_keywords`` dictionary
    section.

    Additionally, new functionality will be added to ScopeSim to resolve the
    kwargs/meta parameters of Effect objects.
    The format for this will be to use a new type: the hash-string.
    This will have this format::

      #<optical_element_name>.<effect_name>.<kwarg_name>

    For example, the temperature of the MICADO detector array can be accessed by::

      '#MICADO_DET.full_detector_array.temperature'

    In the context of the yaml file this would look like

    - ext_type: PrimaryHDU
      keywords:
        HIERARCH:
          ESO:
            DET
              TEMPERAT: '#MICADO_DET.full_detector_array.temperature'

    Obviously some though needs to be put into how exactly we list the simulation
    parameters in a coherent manner.
    But this is 'Zukunftsmusik'.
    For now we really just want an interface that can add the ESO header keywords,
    which can also be expanded in the future for our own purposes.

    Below is an example of some extra keywords for MICADO headers::

        - ext_type: PrimaryHDU
          keywords:
            HIERARCH:
              ESO:
                ATM:
                  TEMPERAT: '!ATMO.temperature'   # will be resolved via from_currsys
                  PWV: '!ATMO.pwv'
                  SEEING: 1.2
                DAR:
                  VALUE: '#<effect_name>.<kwarg_name>'   # will be resolved via effects
                DPR:
                  TYPE: 'some_type'
              SIM:
                random_simulation_keyword: some_value
              MIC:
                micado_specific: ['keyword', 'keyword comment']

          unresolved_keywords:
            HIERARCH:
              ESO:
                ATM:
                  TEMPERAT: '!ATMO.temperature'   # will be left as a string

        - ext_type: ImageHDU
          keywords:
            HIERARCH:
              SIM:
                hello: world
                hallo: welt
                grias_di: woed
                zdrasviute: mir
                salud: el mundo

    """
    def __init__(self, **kwargs):
        params = {"header_dict": None,
                  "filename": None}
        self.meta["z_order"] = [999]
        self.meta.update(params)
        self.meta.update(kwargs)

        self.dict_list = []
        filename = self.meta["filename"]
        if filename is not None:
            with open(filename) as f:
                # possible multiple yaml docs in a file
                # --> returns list even for a single doc
                dics = [dic for dic in yaml.full_load_all(f)]
                for dic in dics:
                    # format says yaml file contains list of dicts
                    if isinstance(dic, list):
                        self.dict_list += dic
                    # catch case where user forgets the list
                    elif isinstance(dic, dict):
                        self.dict_list += [dic]

        if self.meta["header_dict"] is not None:
            self.dict_list += [self.meta["header_dict"]]

    def apply_to(self, hdul, **kwargs):
        if isinstance(hdul, fits.HDUList):
            for dic in self.dict_list:
                resolved = flatten_dict(dic.get("keywords", {}), resolve=Ture)
                unresolved = flatten_dict(dic.get("unresolved_keywords", {}))
                exts = get_relevant_extensions(dic, hdul)
                for i in exts:
                    hdul[i].header.update(resolved)
                    hdul[i].header.update(unresolved)

        return hdlu


def get_relevant_extensions(dic, hdul):
    exts = []
    if dic.get("ext_name") is not None:
        exts += [i for i, hdu in enumerate(hdul)
                 if hdu.header["EXTNAME"] == dic["ext_name"]]
    elif dic.get("ext_number") is not None:
        ext_n = np.array(dic["ext_number"])
        exts += list(ext_n[ext_n<len(hdul)])
    elif dic.get("ext_type") is not None:
        cls = tuple([getattr(fits, cls_str) for cls_str in dic["ext_type"]])
        exts += [i for i, hdu in enumerate(hdul) if isinstance(hdu, cls)]

    return exts


def flatten_dict(dic, base_str="", flat_dict={}, resolve=False):
    for key, val in dic.items():
        new_str = base_str + f"{key} "
        if isinstance(val, dict):
            flatten_dict(val, new_str, flat_dict, resolve)
        else:
            flat_dict[new_str[:-1]] = val

    return flat_dict
