import pytest
from scopesim import rc

run_test = rc.__config__["!SIM.tests.run_integration_tests"]
pytestmark = pytest.mark.skipif(run_test is False,
                                reason="Ignoring MICADO integration tests")