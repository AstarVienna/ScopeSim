from astropy.table import Table
from astropy.io import ascii as ioascii
from astropy.io import fits

from .. import utils


class DataContainer:
    """
    A class to hold data files needed by all Effects objects

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
    If a table is to be generated from an ``array_dict`` parameter, column units
    can be passed as keyword arguments (kwargs) using the following format::

        Datacontainer(... , <column name>_unit="<unit string>")

    where unit string is a string recognised by ``astropy.units``.
    Any additional table meta-data can also be passed using this format.


    Attributes
    ----------
    data : astropy.Table, fits.HDUList
        A generic property method which returns the data from the file. Any
        function calling this should be prepared to handle both data formats

    meta : dict
        Contains all meta data read in from the file's header, and/or passed via
        kwargs

    table : astropy.Table
        If the file has a table format (ASCII of FITS) it is read in
        immediately and stored in ``.table``

    _file : HDUList pointer
        If the file is a FITS image or cube, the data is only read in when
        needed in order to save on memory usage. ``._file`` contains a pointer
        to the data open FITS file.

    """
    def __init__(self, filename=None, table=None, array_dict=None, **kwargs):

        if filename is None and "file_name" in kwargs:
            filename = kwargs["file_name"]

        filename = utils.find_file(filename)
        self.meta = {"filename": filename,
                     "description": "",
                     "history": [],
                     "name": "<empty>"}
        self.meta.update(kwargs)

        self.headers = []
        self.table = None
        self._file = None

        if filename is not None:
            if self.is_fits:
                self._load_fits()
            else:
                self._load_ascii()

        if table is not None:
            self._from_table(table)

        if array_dict is not None:
            self._from_arrays(array_dict)

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
        self.table = ioascii.read(self.meta["filename"],
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
        self.meta["history"] += ["ASCII table read from {}"
                                 "".format(self.meta["filename"])]

    def _load_fits(self):
        self._file = fits.open(self.meta["filename"])
        for ext in self._file:
            self.headers += [ext.header]

        self.meta.update(dict(self._file[0].header))
        self.meta["history"] += ["Opened handle to FITS file {}"
                                 "".format(self.meta["filename"])]

    def get_data(self, ext=0, layer=None):
        """
        Returns either a table or a ImageHDU object

        .. note:: Use this call for reading in individual FITS extensions.

           The ``.data'' handle will read in **all** extensions and return an
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
            flag = utils.is_fits(self.meta["filename"])

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
