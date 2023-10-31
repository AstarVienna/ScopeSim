import pytest
from unittest.mock import patch

import os
from pathlib import Path
import shutil
from docutils.core import publish_doctree

from scopesim.reports import rst_utils as ru
from scopesim.tests.mocks.py_objects import report_objects as ro
from scopesim.tests.mocks.py_objects import effects_objects as eo


CLEAN_UP = True
PLOTS = False

IMG_PATH = Path(__file__).parent / "images_temp"
LATEX_PATH = Path(__file__).parent / "latex_temp"
RST_PATH = Path(__file__).parent / "rst_temp"


@pytest.fixture(scope="module", autouse=True)
def setup_and_teardown():
    paths = [IMG_PATH, LATEX_PATH, RST_PATH]
    for path in paths:
        path.mkdir(parents=True, exist_ok=True)

    patched = {"!SIM.reports.image_path": str(IMG_PATH.absolute()),
               "!SIM.reports.latex_path": str(LATEX_PATH.absolute()),
               "!SIM.reports.rst_path": str(RST_PATH.absolute()),
               }
    with patch.dict("scopesim.rc.__config__", patched):
        with patch.dict("scopesim.rc.__currsys__", patched):
            yield

    if CLEAN_UP:
        for path in paths:
            shutil.rmtree(path)


class TestWalk:
    def test_context_code_is_empty_with_reset_flag(self):
        doctree = publish_doctree(ro.reset_comment_snippet)
        code = ru.walk(doctree)
        assert code == ""


class TestPlotRstText:
    # @pytest.mark.skip(reason=("This produces a DeprecationWarning about a "
    #                           "module called py23. Find out what that is and "
    #                           "remove/replace it."))
    def test_image_file_exists_for_comment_node(self):
        assert IMG_PATH.exists()
        ru.plotify_rst_text(ro.comment_plot_snippet)
        assert (IMG_PATH / "my_fug.png").exists()
        assert (IMG_PATH / "my_fug.pdf").exists()

    # @pytest.mark.skip(reason=("This produces a DeprecationWarning about a "
    #                           "module called py23. Find out what that is and "
    #                           "remove/replace it."))
    def test_image_file_exists_for_comment_node_with_escapable_name(self):
        """Test whether images are created with escapable names.

        That is, on windows, plotify_rst_text should not create
        images_temp\ty_fug.pdf, because that has a tab character in it.
        """
        assert IMG_PATH.exists()
        ru.plotify_rst_text(ro.comment_plot_snippet_with_escapable_name)
        assert (IMG_PATH / "ty_fug.png").exists()
        assert (IMG_PATH / "ty_fug.pdf").exists()

    def test_image_file_exists_for_literal_node(self):
        ru.plotify_rst_text(ro.literal_plot_snippet)
        assert (IMG_PATH / "my_fug3.svg").exists()
        assert (IMG_PATH / "my_fug3.png").exists()


class TestLatexifyRstText:
    def test_stuff(self):
        ru.latexify_rst_text(ro.big_rst_text)
        assert (LATEX_PATH / "This_parrot_goes_vrooom.tex").exists()


class TestRstifyRstText:
    def test_stuff(self):
        ru.rstify_rst_text(ro.big_rst_text)
        assert (RST_PATH / "This_parrot_goes_vrooom.rst").exists()


@pytest.mark.skip(reason="Ignoring for Github Actions")
# @pytest.mark.skip(reason=("This produces a DeprecationWarning about a "
#                           "module called py23. Find out what that is and "
#                           "remove/replace it."))
class TestPlotifyRstText:
    @pytest.mark.parametrize("fname", ["my_fug_A.pdf",
                                       "my_fug_B.svg",
                                       "my_fug_C.png"])
    def test_stuff(self, fname):
        # FIXME: find a better solution (monkeypatch?) than this chdir stuff
        os.chdir(IMG_PATH.parent)
        ru.plotify_rst_text(ro.big_rst_text)
        assert (IMG_PATH / fname).exists()


class TestEffectReport:
    def test_all_parts_are_created_in_rc_folders(self):
        # patched = {"!SIM.reports.image_path": str(IMG_PATH.absolute()),
        #            "!SIM.reports.latex_path": str(LATEX_PATH.absolute()),
        #            "!SIM.reports.rst_path": str(RST_PATH.absolute()),
        #            }
        # with patch.dict("scopesim.rc.__currsys__", patched):
        det_list = eo._detector_list()
        rst_text = det_list.report()
        ru.rstify_rst_text(rst_text, title_char="*", filename=det_list.meta["name"])
        ru.latexify_rst_text(rst_text, title_char="*", filename=det_list.meta["name"])
