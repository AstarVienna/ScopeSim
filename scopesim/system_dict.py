# -*- coding: utf-8 -*-

import logging
from typing import TextIO
from io import StringIO
from collections.abc import Iterable, Mapping, MutableMapping

from more_itertools import ilen


class SystemDict(MutableMapping):
    def __init__(self, new_dict=None):
        self.dic = {}
        if isinstance(new_dict, Mapping):
            self.update(new_dict)
        elif isinstance(new_dict, Iterable):
            for entry in new_dict:
                self.update(entry)

    def update(self, new_dict: MutableMapping) -> None:
        # TODO: why do we check for dict here but not in the else?
        if isinstance(new_dict, Mapping) \
                and "alias" in new_dict \
                and "properties" in new_dict:
            alias = new_dict["alias"]
            if alias in self.dic:
                self.dic[alias] = recursive_update(self.dic[alias],
                                                   new_dict["properties"])
            else:
                self.dic[alias] = new_dict["properties"]
        else:
            # Catch any bang-string properties keys
            to_pop = []
            for key in new_dict:
                if key.startswith("!"):
                    self[key] = new_dict[key]
                    to_pop.append(key)
            for key in to_pop:
                new_dict.pop(key)

            if len(new_dict) > 0:
                self.dic = recursive_update(self.dic, new_dict)

    def __getitem__(self, key):
        if isinstance(key, str) and key.startswith("!"):
            # TODO: these should be replaced with key.removeprefix("!")
            #       once we can finally drop support for Python 3.8 UwU
            key_chunks = key[1:].split(".")
            entry = self.dic
            for key in key_chunks:
                if not isinstance(entry, Mapping):
                    raise KeyError(key)
                entry = entry[key]
            return entry
        return self.dic[key]

    def __setitem__(self, key, value):
        if isinstance(key, str) and key.startswith("!"):
            # TODO: these should be replaced with item.removeprefix("!")
            #       once we can finally drop support for Python 3.8 UwU
            *key_chunks, final_key = key[1:].split(".")
            entry = self.dic
            for key in key_chunks:
                if key not in entry:
                    entry[key] = {}
                entry = entry[key]
            entry[final_key] = value
        else:
            self.dic[key] = value

    def __delitem__(self, key):
        raise NotImplementedError("item deletion is not yet implemented for "
                                  f"{self.__class__.__name__}")

    def _yield_subkeys(self, key, value):
        for subkey, subvalue in value.items():
            if isinstance(subvalue, Mapping):
                yield from self._yield_subkeys(f"{key}.{subkey}", subvalue)
            else:
                yield f"!{key}.{subkey}"

    def __iter__(self):
        for key, value in self.dic.items():
            if isinstance(value, Mapping):
                yield from self._yield_subkeys(key, value)
            else:
                yield key

    def __len__(self) -> int:
        return ilen(iter(self))

    def _write_subdict(self, subdict: Mapping, stream: TextIO,
                       pad: str = "") -> None:
        pre = pad.replace("├─", "│ ").replace("└─", "  ")
        n_sub = len(subdict)
        for i_sub, (key, val) in enumerate(subdict.items()):
            subpre = "└─" if i_sub == n_sub - 1 else "├─"
            stream.write(f"{pre}{subpre}{key}: ")
            if isinstance(val, Mapping):
                self._write_subdict(val, stream, pre + subpre)
            else:
                stream.write(f"{val}")

    def write_string(self, stream: TextIO) -> None:
        """Write formatted string representation to I/O stream"""
        stream.write("SystemDict contents:")
        self._write_subdict(self.dic, stream, "\n")

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.dic!r})"

    def __str__(self) -> str:
        with StringIO() as str_stream:
            self.write_string(str_stream)
            output = str_stream.getvalue()
        return output


def recursive_update(old_dict: MutableMapping, new_dict: Mapping) -> MutableMapping:
    if new_dict is not None:
        for key in new_dict:
            if old_dict is not None and key in old_dict:
                if isinstance(old_dict[key], Mapping):
                    if isinstance(new_dict[key], Mapping):
                        old_dict[key] = recursive_update(old_dict[key],
                                                         new_dict[key])
                    else:
                        logging.warning("Overwriting dict: %s with non-dict: %s",
                                        old_dict[key], new_dict[key])
                        old_dict[key] = new_dict[key]
                else:
                    if isinstance(new_dict[key], Mapping):
                        logging.warning("Overwriting non-dict: %s with dict: %s",
                                        old_dict[key], new_dict[key])
                    old_dict[key] = new_dict[key]
            else:
                old_dict[key] = new_dict[key]

    return old_dict
