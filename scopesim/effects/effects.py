import os
from astropy.table import Table

from ..effects.data_container import DataContainer
from .. import base_classes as bc
from ..utils import from_currsys, write_report
from ..reports.rst_utils import table_to_rst
from .. import rc


class Effect(DataContainer):
    """
    The base class for representing the effects (artifacts) in an optical system

    The ``Effect`` class is conceived to independently apply the changes that
    an optical component (or series thereof) has on an incoming 3D description
    of an on-sky object. In other words, **an Effect object should receive a
    derivative of a ``Source`` object, alter it somehow, and return it**.

    The interface for the Effect base-class has been kept very general so that
    it can easily be sub-classed as data for new effects becomes available.
    Essentially, a sub-classed Effects object must only contain the following
    attributes:

    * ``self.meta`` - a dictionary to contain meta data.
    * ``self.apply_to(obj, **kwargs)`` - a method which accepts a
      Source-derivative and returns an instance of the same class as ``obj``
    * ``self.fov_grid(which="", **kwargs)``


    Parameters
    ----------
    See :class:`DataContainer` for input parameters

    Methods


    """

    def __init__(self, **kwargs):
        super(Effect, self).__init__(**kwargs)
        self.meta["z_order"] = []
        self.meta["include"] = True
        self.meta.update(kwargs)

    def apply_to(self, obj, **kwargs):
        if not isinstance(obj, (bc.FOVSetupBase, bc.SourceBase,
                                bc.FieldOfViewBase, bc.ImagePlaneBase,
                                bc.DetectorBase)):
            raise ValueError(f"object must one of the following: FOVSetupBase, "
                             f"Source, FieldOfView, ImagePlane, Detector: "
                             f"{type(obj)}")

        return obj

    def fov_grid(self, which="", **kwargs):
        """
        Returns the edges needed to generate FieldOfViews for an observation

        Parameters
        ----------
        which : str
            ["waveset", "edges", "shifts"] where:
            * waveset - wavelength bin extremes
            * edges - on sky coordinate edges for each FOV box
            * shifts - wavelength dependent FOV position offsets

        kwargs
        ------
        wave_min, wave_max : float
            [um] list of wavelength

        wave_mid : float
            [um] wavelength what will be centred on optical axis


        Returns
        -------
        waveset : list
            [um] N+1 wavelengths that set edges of N spectral bins

        edges : list of lists
            [arcsec] Contains a list of footprint lists

        shifts : list of 3 lists
            [wave, dx, dy] Contains lists corresponding to the (dx, dy) offset
            from the optical axis (0, 0) induced for each wavelength in (wave)
            [um, arcsec, arcsec]

        """
        self.update(**kwargs)
        return []

    def update(self, **kwargs):
        self.meta.update(kwargs)
        # self.update_bang_keywords()

    # def update_bang_keywords(self):
    #     for key in self.meta:
    #         if isinstance(self.meta[key], str) and self.meta[key][0] == "!":
    #             bang_key = self.meta[key]
    #             self.meta[key] = rc.__currsys__[bang_key]

    @property
    def include(self):
        return from_currsys(self.meta["include"])

    @include.setter
    def include(self, item):
        self.meta["include"] = item

    @property
    def display_name(self):
        return self.meta["name"]

    @property
    def meta_string(self):
        meta_str = ""
        max_key_len = max([len(key) for key in self.meta.keys()])
        for key in self.meta:
            if key not in ["comments", "changes", "description", "history",
                           "report_table_caption", "report_plot_caption",
                           "table"]:
                meta_str += "    {} : {}\n".format(key.rjust(max_key_len),
                                                   self.meta[key])
    
        return meta_str

    def report(self, filename=None, output="rst", rst_title_chars="*+",
               **kwargs):
        """
        For Effect objects, generates a report based on the data and meta-data
    
        This is to aid in the automation of the documentation process of the
        instrument packages in the IRDB.
    
        .. note:: If the Effect can generate a plot, this will be saved to disc
    
        Parameters
        ----------
        filename : str, optional
            Where to save the RST file
        output : str, optional
            ["rst", "latex"] Output file format
        rst_title_chars : 2-str, optional
            Two unique characters used to denote rst subsection headings.
            Options: = - ` : ' " ~ ^ _ * + # < >

        Additional parameters
        ---------------------
        Either from the ``self.meta["report"]`` dictionary or via ``**kwargs``

        "report_table_include": False
        "report_table_caption": ""
        "report_plot_caption": ""
        "report_plot_include": False
        "report_plot_file_formats": ["png"]
            Multiple formats can be saved. The last entry is used for the RST.
        "report_plot_filename": None
            If None, uses self.meta["name"] as the filename
        "file_description": str
            Taken from the header of a file, if available
        "class_description": str
            Taken from the docstring of the subclass
        "changes_str": list of str
            Take from the header of a file, if available

        Returns
        -------
        rst_str : str
            The full reStructureText string
    
        Notes
        -----
    
        The format of the RST output is as follows::
    
            <ClassType>: <effect name>
            **************************
            File Description: <description for file meta data>
            Class Description: <description from class docstring>
            Changes: <list of changes from file meta data> 
    
            Data
            ++++
            .. figure:: <Figure_name>.png
                If the <Effect> object contains a ``.plot()`` function, add
                plot and write it to disc
            Figure caption
    
            Table caption
            Table
                If the <Effect> object contains a ``.table()`` function, add
                a pprint version of the table
    
            Meta-data
            +++++++++
            ::
                A code block print out of the ``.meta`` dictionary
    
    
        """
        changes = self.meta.get("changes", [])
        changes_str = "- " + "\n- ".join([str(entry) for entry in changes])
        cls_doc = self.__doc__ if self.__doc__ is not None else "<no docstring>"
        cls_descr = cls_doc.lstrip().splitlines()[0]

        params = {"report_plot_filename": None,
                  "report_plot_file_formats": ["png"],
                  "report_plot_caption": "",
                  "report_plot_include": False,
                  "report_table_include": False,
                  "report_table_caption": "",
                  "report_table_rounding": None,
                  "report_image_path": "!SIM.reports.image_path",
                  "report_rst_path": "!SIM.reports.rst_path",
                  "report_latex_path": "!SIM.reports.latex_path",
                  "file_description": self.meta.get("description",
                                                    "<no description>"),
                  "class_description": cls_descr,
                  "changes_str": changes_str}
        params.update(self.meta)
        params.update(kwargs)
        params = from_currsys(params)

        rst_str = """
{}
{}
**Included by default**: ``{}``

**File Description**: {}

**Class Description**: {}

**Changes**:

{}

Data
{}
""".format(str(self),
           rst_title_chars[0] * len(str(self)),
           params["include"],
           params["file_description"],
           params["class_description"],
           params["changes_str"],
           rst_title_chars[1] * 4)

        if params["report_plot_include"] and hasattr(self, "plot"):
            fig = self.plot()

            if fig is not None:
                path = params["report_image_path"]
                fname = params["report_plot_filename"]
                if fname is None:
                    fname = self.meta["name"].lower().replace(" ", "_")

                for fmt in params["report_plot_file_formats"]:
                    fname = ".".join((fname.split(".")[0], fmt))
                    file_path = os.path.join(path, fname)

                    fig.savefig(fname=file_path)

                # rel_path = os.path.relpath(params["report_image_path"],
                #                            params["report_rst_path"])
                # rel_file_path = os.path.join(rel_path, fname)

                rst_str += """
.. figure:: {}
    :name: {}

    {}
""".format(fname,
           "fig:" + params.get("name", "<unknown Effect>"),
           params["report_plot_caption"])
    
        if params["report_table_include"]:
            rst_str += """
.. table::
    :name: {}

{}

{}
""".format("tbl:" + params.get("name"),
           table_to_rst(self.table, indent=4,
                        rounding=params["report_table_rounding"]),
           params["report_table_caption"])

        rst_str += """
Meta-data
{}
::

{}
""".format(rst_title_chars[1] * 9,
           self.meta_string)
    
        write_report(rst_str, filename, output)

        return rst_str

    def __repr__(self):
        name = self.meta.get("name", self.meta.get("filename", "<empty>"))
        return '{}: "{}"'.format(type(self).__name__, name)

    def __str__(self):
        return self.__repr__()
