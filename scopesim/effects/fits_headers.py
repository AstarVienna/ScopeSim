import yaml
import numpy as np
from astropy.io import fits
from . import Effect
from ..utils import check_keys


class ExtraFitsKeywords(Effect):
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
