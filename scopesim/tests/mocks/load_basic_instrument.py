import os
from copy import deepcopy
from ... import UserCommands, OpticalTrain
from ... import rc


def load_example_optical_train(**kwargs):
    """
    Returns an basic example ``OpticalTrain`` object with IMG and SPEC modes

    Parameters
    ----------
    Any of the additional parameters taken by UserCommands.
    E.g. ``properties``, ``set_modes``

    Examples
    --------
    ::

        import scopesim as sim
        from scopesim.source import source_templates as st

        src = st.star()
        opt = sim.load_example_optical_train(set_modes=["spectroscopy"])
        opt["filter_wheel"].change_filter("J")

        opt.observe(src)
        hdul = opt.readout()[0]


    Returns
    -------
    opt: OpticalTrain

    """
    # old_local_path = deepcopy(rc.__config__["!SIM.file.local_packages_path"])
    # old_search_path = deepcopy(rc.__search_path__)

    inst_pkgs = os.path.dirname(__file__)
    if inst_pkgs not in rc.__search_path__:
        rc.__search_path__ = rc.__search_path__ + [inst_pkgs]
        rc.__config__["!SIM.file.local_packages_path"] = inst_pkgs

    cmd = UserCommands(use_instrument="basic_instrument", **kwargs)
    opt = OpticalTrain(cmd)

    # rc.__config__["!SIM.file.local_packages_path"] = old_local_path
    # rc.__search_path__ = old_search_path

    return opt
