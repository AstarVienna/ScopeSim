import os
import numpy as np
from scopesim.effects import effects_utils as eu


def setup(app):
    output_dir = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                 "../effects_docstrings"))
    os.makedirs(output_dir, exist_ok=True)

    efs_dict = eu.scopesim_effect_classes()
    eff_types = np.unique([eff.split(".")[0] for eff in efs_dict])
    eff_type_strs = {eff_type: f"{eff_type}\n{'='*len(eff_type)}\n\n"
                     for eff_type in eff_types}

    for eff, cls in efs_dict.items():
        new_str = f"{cls.__name__}\n{'+'*len(cls.__name__)}\n\n"
        if cls.__doc__ is not None:
            new_str += cls.__doc__.replace("\n    ", "\n") + "\n\n"
        else:
            new_str += ".. warning: Empty docstring\n\n"

        base_name = eff.split(".")[0]
        eff_type_strs[base_name] += new_str

    for eff_type, strs in eff_type_strs.items():
        with open(os.path.join(output_dir, f"{eff_type}.rst"), "w") as f:
            f.write(strs)


    toc_str = """Effect Descriptions
===================

.. toctree::
    :maxdepth: 2
    :caption: Contents:
    :glob:

    *
    
"""

    with open(os.path.join(output_dir, "index.rst"), "w") as f:
        f.write(toc_str)
