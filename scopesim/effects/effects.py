# -*- coding: utf-8 -*-
"""Contains base class for effects."""

from pathlib import Path
from collections.abc import Mapping, MutableMapping
from dataclasses import dataclass, field, InitVar, fields
from typing import NewType, ClassVar

from .data_container import DataContainer
from .. import base_classes as bc
from ..utils import from_currsys, write_report
from ..reports.rst_utils import table_to_rst


# FIXME: This docstring is out-of-date for several reasons:
#   - Effects can act on objects other than Source (eg FOV, IMP, DET)
#   - fov_grid is outdated


# @dataclass(kw_only=True, eq=False)
class Effect:
    """
    Base class for representing the effects (artifacts) in an optical system.

    The ``Effect`` class is conceived to independently apply the changes that
    an optical component (or series thereof) has on an incoming 3D description
    of an on-sky object. In other words, **an Effect object should receive a
    derivative of a ``Source`` object, alter it somehow, and return it**.

    The interface for the Effect base-class has been kept very general so that
    it can easily be sub-classed as data for new effects becomes available.
    Essentially, a sub-classed Effects object must only contain the following
    attributes:

    * ``self.meta`` - a dictionary to contain metadata.
    * ``self.apply_to(obj, **kwargs)`` - a method which accepts a
      Source-derivative and returns an instance of the same class as ``obj``
    * ``self.fov_grid(which="", **kwargs)``


    Parameters
    ----------
    See :class:`DataContainer` for input parameters

    """

    z_order: ClassVar[tuple[int, ...]] = tuple()
    required_keys = set()

    def __init__(self, filename=None, **kwargs):
        self.data_container = DataContainer(filename=filename, **kwargs)
        self.meta = kwargs.get("meta", {})
        self.cmds = kwargs.get("cmds")

        self.meta.update(self.data_container.meta)
        self.meta["include"] = True
        self.meta.update(kwargs)

    def apply_to(self, obj, **kwargs):
        """TBA."""
        if not isinstance(obj, (bc.FOVSetupBase, bc.SourceBase,
                                bc.FieldOfViewBase, bc.ImagePlaneBase,
                                bc.DetectorBase)):
            raise ValueError("object must one of the following: FOVSetupBase, "
                             "Source, FieldOfView, ImagePlane, Detector: "
                             f"{type(obj)}")

        return obj

    def fov_grid(self, which="", **kwargs):
        """
        Return the edges needed to generate FieldOfViews for an observation.

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
        # self.update(**kwargs)
        return []

    # *******************************************************************
    # ported from DataContainer, previous base class
    # *******************************************************************

    @property
    def table(self):
        return self.data_container.table

    @table.setter
    def table(self, value):
        self.data_container.table = value

    @property
    def data(self):
        return self.data_container.data

    @property
    def _file(self):
        return self.data_container._file

    @_file.setter
    def _file(self, value):
        self.data_container._file = value

    # *******************************************************************

    @property
    def include(self) -> bool:
        return from_currsys(self.meta["include"], self.cmds)

    @include.setter
    def include(self, value: bool):
        self.meta["include"] = value

    @property
    def display_name(self) -> str:
        name = self.meta.get("name", self.meta.get("filename", "<untitled>"))
        if not hasattr(self, "_current_str"):
            return name
        current_str = from_currsys(self.meta[self._current_str], self.cmds)
        return f"{name} : [{current_str}]"

    @property
    def meta_string(self) -> str:
        padlen = 4 + len(max(self.meta, key=len))
        exclude = {"comments", "changes", "description", "history",
                   "report_table_caption", "report_plot_caption", "table"}
        meta_str = "\n".join(f"{key:>{padlen}} : {value}"
                             for key, value in self.meta.items()
                             if key not in exclude)
        return meta_str

    def report(self, filename=None, output="rst", rst_title_chars="*+",
               **kwargs):
        """
        For Effect objects, generates a report based on the data and meta-data.

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
        changes_str = "- " + "\n- ".join(str(entry) for entry in changes)
        cls_doc = self.__doc__ if self.__doc__ is not None else "<no docstring>"
        cls_descr = cls_doc.lstrip().splitlines()[0]

        params = {
            "report_plot_filename": None,
            "report_plot_file_formats": ["png"],
            "report_plot_caption": "",
            "report_plot_include": getattr(self, "report_plot_include", False),
            "report_table_include": getattr(self, "report_table_include", False),
            "report_table_caption": "",
            "report_table_rounding": getattr(self, "report_table_rounding", None),
            "report_image_path": "!SIM.reports.image_path",
            "report_rst_path": "!SIM.reports.rst_path",
            "report_latex_path": "!SIM.reports.latex_path",
            "file_description": self.meta.get("description",
                                              "<no description>"),
            "class_description": cls_descr,
            "changes_str": changes_str,
        }
        params.update(self.meta)
        params.update(kwargs)
        params = from_currsys(params, self.cmds)

        rst_str = f"""
{str(self)}
{rst_title_chars[0] * len(str(self))}
**Included by default**: ``{params["include"]}``

**File Description**: {params["file_description"]}

**Class Description**: {params["class_description"]}

**Changes**:

{params["changes_str"]}

Data
{rst_title_chars[1] * 4}
"""

        if params["report_plot_include"] and hasattr(self, "plot"):
            from matplotlib.figure import Figure
            fig = self.plot()
            # HACK: plot methods should always return the same, while this is
            #       not sorted out, deal with both fig and ax
            if not isinstance(fig, Figure):
                fig = fig.figure

            if fig is not None:
                path = params["report_image_path"]
                fname = params["report_plot_filename"]
                if fname is None:
                    fname = self.meta["name"].lower().replace(" ", "_")

                for fmt in params["report_plot_file_formats"]:
                    fname = ".".join((fname.split(".")[0], fmt))
                    file_path = Path(path, fname)

                    fig.savefig(fname=file_path)

                # rel_path = os.path.relpath(params["report_image_path"],
                #                            params["report_rst_path"])
                # rel_file_path = os.path.join(rel_path, fname)

                # TODO: fname is set in a loop above, so using it here in the
                #       fstring will only access the last value from the loop,
                #       is that intended?
                rst_str += f"""
.. figure:: {fname}
    :name: {"fig:" + params.get("name", "<unknown Effect>")}

    {params["report_plot_caption"]}
"""

        if params["report_table_include"]:
            rst_str += f"""
.. table::
    :name: {"tbl:" + params.get("name")}

{table_to_rst(self.table, indent=4, rounding=params["report_table_rounding"])}

{params["report_table_caption"]}
"""

        rst_str += f"""
Meta-data
{rst_title_chars[1] * 9}
::

{self.meta_string}
"""

        write_report(rst_str, filename, output)

        return rst_str

    def info(self) -> None:
        """Print basic information on the effect, notably the description."""
        if (desc := self.meta.get("description")) is not None:
            print(f"{self}\nDescription: {desc}")
        else:
            print(self)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(**{self.meta!r})"

    def __str__(self) -> str:
        return f"{self.__class__.__name__}: \"{self.display_name}\""

    def __getitem__(self, item):
        if isinstance(item, str) and item.startswith("#"):
            if len(item) > 1:
                if item.endswith("!"):
                    key = item.removeprefix("#").removesuffix("!")
                    if len(key) > 0:
                        value = from_currsys(self.meta[key], self.cmds)
                    else:
                        value = from_currsys(self.meta, self.cmds)
                else:
                    value = self.meta[item.removeprefix("#")]
            else:
                value = self.meta
        else:
            raise ValueError(f"__getitem__ calls must start with '#': {item}")

        return value

    def _get_path(self):
        if any(key not in self.meta for key in ("path", "filename_format")):
            return None
        return Path(self.meta["path"],
                    from_currsys(self.meta["filename_format"], self.cmds))
