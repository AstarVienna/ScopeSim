from warnings import warn

from astropy.table import Table
from astropy.io import ascii as ioascii
from astropy.io import fits

from .. import utils


class DataContainer:
    """
    A class to hold data files needed by all Effects objects.

    Parameters
    ----------
    filename : str
        Path to file containing data.
        Accepted formats: ASCII table, FITS table, FITS image

    table : astropy.Table
        An astropy Table containing data

    array_dict : dict
        A dictionary out of which an astropy.Table object can be constructed.

    kwargs :
        addition meta data


    Notes
    -----
    If a table is to be generated from an ``array_dict`` parameter, column
    units can be passed as keyword arguments (kwargs) using the following
    format::

        Datacontainer(... , <column name>_unit="<unit string>")

    where unit string is a string recognised by ``astropy.units``.
    Any additional table meta-data can also be passed using this format.


    Attributes
    ----------
    data : astropy.Table, fits.HDUList
        A generic property method which returns the data from the file. Any
        function calling this should be prepared to handle both data formats

    meta : dict
        Contains all meta data read in from the file's header, and/or passed
        via kwargs.

    table : astropy.Table
        If the file has a table format (ASCII of FITS) it is read in
        immediately and stored in ``.table``

    _file : HDUList pointer
        If the file is a FITS image or cube, the data is only read in when
        needed in order to save on memory usage. ``._file`` contains a pointer
        to the data open FITS file.

    """

    meta = None

    def __init__(self, filename=None, table=None, array_dict=None, cmds=None,
                 **kwargs):
        self.cmds = cmds
        # Setting a default for cmds cannot be done here, because from_currsys
        # checks whether cmds is None. TODO: make this possible.
        # if self.cmds is None:
        #     from scopesim import UserCommands
        #     self.cmds = UserCommands()

        if filename is None and "file_name" in kwargs:
            warn("The 'file_name' kwarg is deprecated and will raise an error "
                 "in the future, please use 'filename' instead!",
                 DeprecationWarning, stacklevel=2)
            filename = kwargs["file_name"]

        # A !-string filename is immediately resolved, but the file is not yet
        # searched for. This makes sense, as the filename should end up in the
        # FITS headers, and should therefor be independent of any particular
        # system.
        #
        # It would perhaps be even better to not resolve a !-string filename
        # here, but that currently does not make sense because the file is
        # directly read in. That is, the file would not be read again if
        # someone changes the value that the !-string points to. So changing
        # the value a !-string points to would lead to an inconsistent state.
        # Immediately resolving the !-string prevents such an inconsistency.
        if isinstance(filename, str) and filename.startswith("!"):
            filename = utils.from_currsys(filename, self.cmds)

        # A derived clas might have set .meta before calling super().__init__()
        if self.meta is None:
            self.meta = {}
        self.meta.update({
            "filename": filename,
            "description": "",
            "history": [],
            "name": "<empty>",
        })
        self.meta.update(kwargs)

        self.headers = []
        self.table = None
        self._file = None

        # Need to check whether the file exists before trying to load it.
        filename_full = utils.find_file(self.meta["filename"])
        if filename_full is not None:
            if self.is_fits:
                self._load_fits()
            else:
                self._load_ascii()

        if table is not None:
            self._from_table(table)

        if array_dict is not None:
            self._from_arrays(array_dict)

        # Get units as they are given in the kworgs.
        unit_dict_kwargs = {
            k: v
            for k, v in kwargs.items()
            if k.endswith("_unit")
        }

        # Collect the units from the table. self.meta should contain both
        # kwargs and headers from the file at this point.
        unit_dict = {
            k: v
            for k, v in self.meta.items()
            if k.endswith("_unit")
        }

        # Verify that any units given in the kwargs
        # are the same as those given in the table.
        for k, v_kwargs in unit_dict_kwargs.items():
            v_table = unit_dict.get(k, v_kwargs)
            if v_kwargs != v_table:
                raise ValueError(f"Column {k} has unit {v_table} in table, but {v_kwargs} in kwargs")

        # Add units from kwargs if they are given.
        if self.table:
            for column in self.table.columns.values():
                key_unit = f"{column.name}_unit"
                if key_unit not in unit_dict:
                    continue
                unit_kwargs = unit_dict[key_unit]
                if column.unit is None:
                    column.unit = unit_kwargs

    def _from_table(self, table):
        self.table = table
        self.headers += [table.meta]
        self.meta.update(table.meta)
        self.meta["history"] += ["Table added directly"]

    def _from_arrays(self, array_dict):
        data = []
        colnames = []
        for key, val in array_dict.items():
            data += [val]
            colnames += [key]
        self.table = Table(names=colnames, data=data)
        self.headers += [None]
        self.meta["history"] += ["Table generated from arrays"]
        self.table.meta.update(self.meta)

    def _load_ascii(self):
        filename_full = utils.find_file(self.meta["filename"])
        self.table = ioascii.read(filename_full,
                                  format="basic", guess=False)
        hdr_dict = utils.convert_table_comments_to_dict(self.table)
        if isinstance(hdr_dict, dict):
            self.headers += [hdr_dict]
        else:
            self.headers += [None]

        self.meta.update(self.table.meta)
        self.meta.update(hdr_dict)
        # self.table.meta.update(hdr_dict)
        self.table.meta.update(self.meta)
        if self.table.meta.get("cmds"):
            self.table.meta.pop("cmds")

        self.meta["history"] += ["ASCII table read from "
                                 f"{self.meta['filename']}"]

    def _load_fits(self):
        filename_full = utils.find_file(self.meta["filename"])
        self._file = fits.open(filename_full)
        for ext in self._file:
            self.headers += [ext.header]

        self.meta.update(dict(self._file[0].header))
        self.meta["history"] += ["Opened handle to FITS file "
                                 f"{self.meta['filename']}"]

    def get_data(self, ext=0, layer=None):
        """
        Return either a table or a ImageHDU object.

        .. note:: Use this call for reading in individual FITS extensions.

           The ``.data`` handle will read in **all** extensions and return an
           HDUList object

        Parameters
        ----------
        ext : int
        layer : int
            If the FITS extension is a data cube, layer corresponds to a slice
            from this cube of ``<ImageHDU>.data[layer, :, :]``

        Returns
        -------
        data_set : astropy.Table, fits.ImageHDU

        """
        data_set = None
        if self.is_fits:
            if isinstance(self._file[ext], fits.BinTableHDU):
                data_set = Table.read(self._file[ext], format="fits")
            else:
                if self._file[ext].data is not None:
                    data_dims = len(self._file[ext].data.shape)
                    if data_dims == 3 and layer is not None:
                        data_set = self._file[ext].data[layer]
                    else:
                        data_set = self._file[ext].data
        else:
            data_set = self.table

        return data_set

    @property
    def is_fits(self):
        flag = False
        if isinstance(self._file, fits.HDUList):
            flag = True
        elif isinstance(self.meta["filename"], str):
            filename_full = utils.find_file(self.meta["filename"])
            flag = utils.is_fits(filename_full)

        return flag

    @property
    def data(self):
        data_set = None
        if self.is_fits:
            for ii in range(len(self._file)):
                data_set = self.get_data(ii)
                if data_set is not None:
                    break
        else:
            data_set = self.table

        return data_set

    def validate(self, etype):
        etype_colname = utils.real_colname("ETYPE", self.meta.colnames)
        return self.meta[etype_colname] == etype
