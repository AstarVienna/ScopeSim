import os
import scopesim as sim


def test_simplecado():

    YAMLS = os.path.abspath(os.path.join(__file__, "../../mocks/yamls/"))

    src = sim.source.templates.empty_sky()

    cmd = sim.commands.UserCommands(sim_data_dir=YAMLS)
    cmd["SIM_DETECTOR_YAML"] = "SimpleCADO.yaml"
    cmd["SIM_PIXEL_SCALE"] = 0.004

    opt = sim.optics.optical_train.OpticalTrain(cmd)
    opt.observe(src)

    print(opt.image_plane.image)

