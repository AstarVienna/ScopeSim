import os

from astropy.table import TableFormatter
from docutils.core import publish_doctree, publish_parts
from docutils.nodes import comment, literal_block
import yaml

from .. import rc
from ..utils import from_currsys


def walk(node, context_code=None):
    """
    Recursively walk through a docutils doctree and run/plot code blocks

    Parameters
    ----------
    node : docutils.node.Node
    context_code : str, optional
        A code string inherited from previous code/comment nodes

    Returns
    -------
    context_code : str
        Code to be inherited by subsequent code/comment nodes

    """

    if context_code is None:
        context_code = """
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
"""

    if isinstance(node, (comment, literal_block)):
        if isinstance(node, comment):
            context_code = process_comment_code(node, context_code)
        elif isinstance(node, literal_block):
            context_code = process_literal_code(node, context_code)

        exec(context_code)

    elif hasattr(node, "children"):
        for child in node.children:
            context_code = walk(child, context_code)

    return context_code


def process_comment_code(node, context_code):
    """Add code from a ``comment`` node to the context_code string"""
    result = node.rawsource.split("---")
    options = yaml.full_load(result[0]) if len(result) > 1 else {}
    new_code = result[-1]
    context_code = process_code(context_code, new_code, options)

    return context_code


def process_literal_code(node, context_code):
    """Add code from a ``literal_block`` node to the context_code string"""
    new_code = node.rawsource
    attribs = node.attributes["classes"]
    action = [att for att in attribs if "format-" not in att]
    format = [att.replace("format-", "") for att in attribs if "format-" in att]
    format = format if len(format) > 0 else ["png"]
    options = {"action": action,
               "format": format}
    if len(node.attributes["names"]) > 0:
        options["name"] = node.attributes["names"][0]
    context_code = process_code(context_code, new_code, options)

    return context_code


def process_code(context_code, code, options):
    """
    Extracts and adds code from the node text to the context_code string

    Code can be passed in either a ``literal_block`` or ``comment`` docutils
    node. The options regarding what to do with the code are included in either
    the :class: and :name: tag of a ``literal_block`` node, or in a yaml header
    section in a ``comment`` node.

    See below for examples of code blocks.

    Options for controlling what happens to the code block are as follows:

    - name: The filename for the plot

    - format: Any matplotlib-accepted file format, e.g: png, pdf, svg, etc.
        Multiple file formats can be

    - action: [reset, clear-figure, plot]
        - ``reset`` clears the context_code string. By default code is saved
          from previous code blocks.
        - ``clear-figure`` adds ``plt.clf()`` to the context string before the
          current code block is added
        - ``plot`` adds ``plt.savefig({name}.{format}) to the context string.
          If multiple formats are given, these will be iterated over.

    For ``comment`` blocks, options should be given in a yaml style header,
    separated form the code block by exactly three (3) hyphens ('---')

    For ``literal_block`` blocks, we hijack the ``:class:`` and ``:name:``
    attributes. See the examples below. All action keywords are passed to
    ``:class:``. Format keys are passed as ``format-<key>``.

    Parameters
    ----------
    context_code : str
        code from previous Nodes
    code : str
        code from the current Node
    options : dict
        options dictionary derived from Node attributes or RST header blocks

    Returns
    -------
    context_code : str

    Examples
    --------

    Example of ``literal_block`` code block::

        .. code::
            :name: my_fug
            :class: reset, clear-figure, plot, format-png

            plt.plot([0,1], [1,1])

        .. figure:: my_fug.png
            :name: fig:my_fug

    Example of a ``comment`` code block::

        ..
            name: my_fug2
            format: [jpg, svg]
            action: [reset, clear-figure, plot]
            ---
            plt.plot([0,1], [1,0])

        .. figure:: my_fug2.jpg
            :name: fig:my_fug2

    """

    if "reset" in options.get("action", []):
        context_code = ""

    if "clear-figure" in options.get("action", []):
        context_code += "\nimport matplotlib.pyplot as plt\nplt.clf()\n"

    if "execute" in options.get("action", []):
        context_code += code

    if "plot" in options.get("action", []):
        context_code += code

        img_path = options.get("path", rc.__config__["!SIM.reports.image_path"])
        formats = options.get("format", ["png"])
        formats = [formats] if isinstance(formats, str) else formats
        for fmt in formats:

            fname = options.get("name", "untitled").split(".")[0]
            fname = ".".join([fname, fmt])
            fname = os.path.join(img_path, fname)
            context_code += '\nplt.savefig("{}")'.format(fname)

    return context_code


def plotify_rst_text(rst_text):
    """
    Generates and saves plots from code blocks in an RST string

    The save directory for the plots defaults to
    ``scopesim.rc.__config__["!SIM.reports.image_path"]``.
    This can be overridden ONLY inside a COMMENT block using the path keyword.

    Parameters
    ----------
    rst_text : str
        Any RST text string

    Examples
    --------
    The following rst text will generate a plot and save it in three formats::

        A basic plot
        ============
        Let's make a basic plot using the comment style

        ..
            action: plot
            format: [pdf, png]
            name: my_fug
            path: "./images/"
            ---
            plt.plot([0,1], [0,1])

        And now a plot using the code block style

        .. code::
            :name: my_fug3
            :class: reset, plot, format-jpg, format-svg

            import matplotlib.pyplot as plt
            plt.plot([0,1], [1,1])

    The plot are created be calling::

        plotify_rst_text(rst_text)

    where ``rst_text`` is a string holding the full RST text from above.

    Notes
    -----
    * Possible actions are: [reset, clear-figure, plot]
    * Code is retained between code blocks in the same string, so we do not need
      to re-write large sections of code
    * By default ``import numpy as np`` and ``import matplotlib.pyplot as plt``
      are loaded automatically, so these do not need to be explicitly specified
      in each code block.
    * THE EXCEPTION is when the action ``reset`` is specified. This clears the
      code_context variable.

    """

    walk(publish_doctree(rst_text))


def latexify_rst_text(rst_text, filename=None, path=None, title_char="=",
                      float_figures=True, use_code_box=True):
    """
    Converts an RST string (block of text) into a LaTeX string

    NOTE: plots will NOT be generated with this command. For that we must invoke
    the ``plotfiy`` command.

    Parameters
    ----------
    rst_text : str
    filename : str, optional
        What to name the latex file. If None, the filename is derived from the
        text title
    path : str, optional
        Where to save the latex file
    title_char : str, optional
        The character used to underline the rst text title. Usually "=".
    float_figures : bool, optional
        Set to False if figures should not be placed by LaTeX.
        Replaces all ``\begin{figure}`` with ``\begin{figure}[H]``
    use_code_box : bool, optional
        Adds a box around quote blocks

    Returns
    -------
    tex_str : str
        The same string or block of text in LaTeX format

    Examples
    --------
    ::
        rst_text = '''
        Meaning of life
        ===============

        Apparently it's 42. Lets plot

        ..
            action: plot
            name: my_plot
            path: ./images/
            ---
            plt.plot([0,1], [0,1])

        .. figure:: images/my_plot.png
            :name: label-my-plot

            This is an included figure caption

        '''

        plotify_rst_text(rst_text)
        latexify_rst_text(rst_text, filename="my_latex_file", path="./")

    """

    if path is None:
        path = from_currsys(rc.__config__["!SIM.reports.latex_path"])

    if filename is None:
        filename = rst_text.split(title_char)[0].strip().replace(" ", "_")

    text = "Title\n<<<<<\nSubtitle\n>>>>>>>>\n\n"
    parts = publish_parts(text + rst_text, writer_name="latex")

    if not float_figures:
        parts["body"] = parts["body"].replace('begin{figure}',
                                              'begin{figure}[H]')

    if use_code_box:
        parts["body"] = parts["body"].replace('begin{alltt}',
                                              'begin{alltt}\n\\begin{lstlisting}[frame=single]')
        parts["body"] = parts["body"].replace('end{alltt}',
                                              'end{lstlisting}\n\\end{alltt}')

    filename = filename.split(".")[0] + ".tex"
    file_path = os.path.join(path, filename)
    with open(file_path, "w") as f:
        f.write(parts["body"])

    tex_str = parts["body"]

    return tex_str


def rstify_rst_text(rst_text, filename=None, path=None, title_char="="):
    """ The same as ``latexify_rst_text```, but the output is in RST format """
    if path is None:
        path = from_currsys(rc.__config__["!SIM.reports.rst_path"])

    if filename is None:
        filename = rst_text.split(title_char)[0].strip().replace(" ", "_")

    filename = filename.split(".")[0] + ".rst"
    file_path = os.path.join(path, filename)
    with open(file_path, "w") as f:
        f.write(rst_text)

    return rst_text


def table_to_rst(tbl, indent=0, rounding=None):
    if isinstance(rounding, int):
        for col in tbl.itercols():
            if col.info.dtype.kind == 'f':
                col.info.format = '.{}f'.format(rounding)
    
    tbl_fmtr = TableFormatter()
    lines, outs = tbl_fmtr._pformat_table(tbl, max_width=-1, max_lines=-1,
                                          show_unit=False)
    i = outs["n_header"] - 1
    lines[i] = lines[i].replace("-", "=")
    lines = [lines[i]] + lines + [lines[i]]

    indent = " " * indent
    rst_str = indent + ("\n" + indent).join(lines)

    return rst_str
