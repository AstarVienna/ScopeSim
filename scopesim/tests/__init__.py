import os
from scopesim import rc

rc.__config__["!SIM.tests.run_integration_tests"] = True
rc.__config__["!SIM.tests.run_skycalc_ter_tests"] = True
rc.__config__["!SIM.file.use_cached_downloads"] = False
rc.__config__["!SIM.reports.ip_tracking"] = False

if "TRAVIS" in os.environ:
    rc.__config__["!SIM.tests.run_skycalc_ter_tests"] = True
    rc.__config__["!SIM.tests.run_integration_tests"] = True
