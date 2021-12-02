# takes a package name as input
# reads in all yaml files within the package
# creates rst descriptions for each optical element
# shows the system throughput
# list all the available modes
# list the dependencies
# produces all the plots
# produces all the rst files
# produces a top-level rst / latex file that references all the generated ones
# can read in preambles / other package specific descriptions


class ReportGenerator:
    def __init__(self, pkg_name):
        self.pkg_name = pkg_name

    def rst_default_yaml(self):
        pass



