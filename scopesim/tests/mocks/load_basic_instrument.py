from unittest.mock import patch

from pathlib import Path
from scopesim import UserCommands, OpticalTrain, Simulation


def load_example_optical_train(**kwargs):
    """
    Return an basic example ``OpticalTrain`` object with IMG and SPEC modes.

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
    patched = {"!SIM.file.local_packages_path": str(Path(__file__).parent)}
    with patch.dict("scopesim.rc.__config__", patched):
        cmd = UserCommands(use_instrument="basic_instrument", **kwargs)
        opt = OpticalTrain(cmd)

    return opt


def example_simulation(*args, **kwargs):
    """Return basic example ``Simulation`` object with IMG and SPEC modes."""
    patched = {"!SIM.file.local_packages_path": str(Path(__file__).parent)}
    with patch.dict("scopesim.rc.__config__", patched):
        sim = Simulation("basic_instrument", *args, **kwargs)

    return sim
