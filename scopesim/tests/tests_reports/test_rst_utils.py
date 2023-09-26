import pytest
from pathlib import Path
import shutil
from docutils.core import publish_doctree

from scopesim.reports import rst_utils as ru
from scopesim.tests.mocks.py_objects import report_objects as ro
from scopesim import rc
rc.__config__["!SIM.reports.image_path"] = "./images_temp/"
rc.__config__["!SIM.reports.latex_path"] = "./latex_temp/"
rc.__config__["!SIM.reports.rst_path"] = "./rst_temp/"

CLEAN_UP = True
PLOTS = False

IMG_PATH = Path(rc.__config__["!SIM.reports.image_path"])
LATEX_PATH = Path(rc.__config__["!SIM.reports.latex_path"])
RST_PATH = Path(rc.__config__["!SIM.reports.rst_path"])


def setup_module():
    for path in [IMG_PATH, LATEX_PATH, RST_PATH]:
        if not path.exists():
            path.mkdir()
    rc.__config__["!SIM.reports.image_path"] = str(IMG_PATH)
    rc.__config__["!SIM.reports.latex_path"] = str(LATEX_PATH)
    rc.__config__["!SIM.reports.rst_path"] = str(RST_PATH)


def teardown_module():
    for path in [IMG_PATH, LATEX_PATH, RST_PATH]:
        if CLEAN_UP:
            shutil.rmtree(path)


class TestWalk:
    def test_context_code_is_empty_with_reset_flag(self):
        doctree = publish_doctree(ro.reset_comment_snippet)
        code = ru.walk(doctree)
        assert code == ""


class TestPlotRstText:
    @pytest.mark.skip(reason=("This produces a DeprecationWarning about a "
                              "module called py23. Find out what that is and "
                              "remove/replace it."))
    def test_image_file_exists_for_comment_node(self):
        assert IMG_PATH.exists()
        ru.plotify_rst_text(ro.comment_plot_snippet)
        assert (IMG_PATH / "my_fug.png").exists()
        assert (IMG_PATH / "my_fug.pdf").exists()

    @pytest.mark.skip(reason=("This produces a DeprecationWarning about a "
                              "module called py23. Find out what that is and "
                              "remove/replace it."))
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
        print(IMG_PATH)
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
class TestPlotifyRstText:
    def test_stuff(self):
        print(IMG_PATH)
        ru.plotify_rst_text(ro.big_rst_text)
        fnames = ["my_fug_A.pdf", "my_fug_B.svg", "my_fug_C.png"]
        for fname in fnames:
            assert (IMG_PATH / fname).exists()


class TestEffectReport:
    def test_all_parts_are_created_in_rc_folders(self):
        from scopesim.tests.mocks.py_objects import effects_objects as eo
        det_list = eo._detector_list()
        rst_text = det_list.report()
        ru.rstify_rst_text(rst_text, title_char="*", filename=det_list.meta["name"])
        ru.latexify_rst_text(rst_text, title_char="*", filename=det_list.meta["name"])
