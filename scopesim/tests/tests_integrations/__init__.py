import pytest
from scopesim import rc

if rc.__config__["!SIM.tests.run_integration_tests"] is False:
    pytestmark = pytest.mark.skip("Ignoring MICADO integration tests")
