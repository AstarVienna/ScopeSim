import telescopy as tsp

def test_optics_spectrograph():
    assert tsp.optics.spectrograph.optics_spectrograph() == 42