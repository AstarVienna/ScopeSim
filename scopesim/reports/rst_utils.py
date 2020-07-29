import os
from docutils.core import publish_doctree, publish_parts
from docutils.nodes import comment, literal_block
import yaml

from .. import rc


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
    new_code = result[-1]
    options = yaml.load(result[0]) if len(result) > 1 else {}
    context_code = process_code(context_code, new_code, options)

    return context_code


def process_literal_code(node, context_code):
    """Add code from a ``literal_block`` node to the context_code string"""
    new_code = node.rawsource
    attribs = node.attributes["classes"]
    action = [att for att in attribs if "format-" not in att]
    format = [att.replace("format-", "") for att in attribs if "format-" in att]
    options = {"action": action,
               "format": format,
               "name": node.attributes["names"][0]}
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
          If more than one formats are given, these will be iterated over.

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
        context_code += "\nplt.clf()\n"

    if "plot" in options.get("action", []):
        context_code += code

        formats = options.get("format", ["png"])
        formats = [formats] if isinstance(formats, str) else formats
        for fmt in formats:
            img_path = rc.__config__["!SIM.reports.image_path"]
            fname = options.get("name", "untitled").split(".")[0]
            fname = ".".join([fname, fmt])
            fname = os.path.join(img_path, fname)
            context_code += '\nplt.savefig("{}")'.format(fname)

    return context_code


def plotify_rst_text(rst_text):
    doctree = publish_doctree(rst_text)
    walk(doctree)


def latexify_rst_text(rst_text, filename=None, path=None):
    if path is None:
        path = rc.__config__["!SIM.reports.latex_path"]

    if filename is None:
        filename = rst_text.split("===")[0].strip().replace(" ", "_") + ".tex"

    text = "Title\n<<<<<\nSubtitle\n>>>>>>>>\n\n"
    parts = publish_parts(text + rst_text, writer_name="latex")

    with open(os.path.join(path, filename), "w") as f:
        f.write(parts["body"])


def rstify_rst_text(rst_text, filename=None, path=None):
    if path is None:
        path = rc.__config__["!SIM.reports.rst_path"]

    if filename is None:
        filename = rst_text.split("===")[0].strip().replace(" ", "_") + ".rst"

    with open(os.path.join(path, filename), "w") as f:
        f.write(rst_text)
