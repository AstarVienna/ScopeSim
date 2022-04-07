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
    old_lpp = deepcopy(rc.__currsys__["!SIM.file.local_packages_path"])
    inst_pkgs = os.path.dirname(__file__)
    rc.__currsys__["!SIM.file.local_packages_path"] = inst_pkgs

    cmd = UserCommands(use_instrument="basic_instrument", **kwargs)
    opt = OpticalTrain(cmd)

    rc.__currsys__["!SIM.file.local_packages_path"] = old_lpp

    return opt