import logging


class SystemDict():
    def __init__(self, new_dict=None):
        self.dic = {}
        if isinstance(new_dict, dict):
            self.update(new_dict)
        elif isinstance(new_dict, list):
            for entry in new_dict:
                self.update(entry)

    def update(self, new_dict):
        if isinstance(new_dict, dict) \
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

    def __getitem__(self, item):
        if isinstance(item, str) and item.startswith("!"):
            item_chunks = item[1:].split(".")
            entry = self.dic
            for item in item_chunks:
                entry = entry[item]
            return entry
        return self.dic[item]

    def __setitem__(self, key, value):
        if isinstance(key, str) and key.startswith("!"):
            key_chunks = key[1:].split(".")
            entry = self.dic
            for key in key_chunks[:-1]:
                if key not in entry:
                    entry[key] = {}
                entry = entry[key]
            entry[key_chunks[-1]] = value
        else:
            self.dic[key] = value

    def __contains__(self, item):
        if isinstance(item, str) and item.startswith("!"):
            item_chunks = item[1:].split(".")
            entry = self.dic
            for item in item_chunks:
                if not isinstance(entry, dict) or item not in entry:
                    return False
                entry = entry[item]
            return True
        return item in self.dic

    def __repr__(self):
        msg = "<SystemDict> contents:"
        for key, val in self.dic.items():
            msg += f"\n{key}: "
            if isinstance(val, dict):
                for subkey in val.keys():
                    msg += f"\n  {subkey}: {val[subkey]}"
            else:
                msg += f"{val}\n"
        return msg


def recursive_update(old_dict, new_dict):
    if new_dict is not None:
        for key in new_dict:
            if old_dict is not None and key in old_dict:
                if isinstance(old_dict[key], dict):
                    if isinstance(new_dict[key], dict):
                        old_dict[key] = recursive_update(old_dict[key],
                                                         new_dict[key])
                    else:
                        logging.warning("Overwriting dict: %s with non-dict: %s",
                                        old_dict[key], new_dict[key])
                        old_dict[key] = new_dict[key]
                else:
                    if isinstance(new_dict[key], dict):
                        logging.warning("Overwriting non-dict: %s with dict: %s",
                                        old_dict[key], new_dict[key])
                    old_dict[key] = new_dict[key]
            else:
                old_dict[key] = new_dict[key]

    return old_dict
