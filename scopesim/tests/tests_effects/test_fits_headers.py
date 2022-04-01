import os
import pytest
from pytest import raises
from astropy.io import fits
import numpy as np

from scopesim.effects import ExtraFitsKeywords, EffectsMetaKeywords
from scopesim.effects import fits_headers as fh
import scopesim as sim

YAMLS_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../mocks/yamls/"))

@pytest.fixture(scope="function")
def simplecado_opt():
    simplecado_yaml = os.path.join(YAMLS_PATH, "SimpleCADO.yaml")
    cmd = sim.UserCommands(yamls=[simplecado_yaml])
    return sim.OpticalTrain(cmd)


@pytest.fixture(scope="function")
def comb_hdul():
    pri = fits.PrimaryHDU(header=fits.Header({"EXTNAME": "PriHDU"}))
    im = fits.ImageHDU(header=fits.Header({"EXTNAME": "ImHDU"}))
    tbl = fits.BinTableHDU(header=fits.Header({"EXTNAME": "BinTblHDU"}))
    hdul = fits.HDUList([pri, im, tbl])

    return hdul


# taken from the example yaml in docstring of ExtraFitsKeywords
@pytest.fixture(scope="function")
def yaml_string():
    return """
- ext_type: PrimaryHDU
  keywords:
    HIERARCH:
      ESO:
        TEL:
          area: ["!TEL.area", "default = 0"]
          pixel_scale: "!INST.pixel_scale"
          its_over: 9000
        DAR:
          VALUE: '#dark_current.value'   # will be resolved via effects
        DPR:
          TYPE: 'some_type'
      SIM:
        random_simulation_keyword: some_value
      MIC:
        micado_specific: ['keyword', 'keyword comment']

  unresolved_keywords:
    HIERARCH:
      ESO:
        ATM:
          TEMPERAT: '!ATMO.temperature'   # will be left as a string

- ext_type: ImageHDU
  keywords:
    HIERARCH:
      SIM:
        hello: world
        hallo: welt
        grias_di: woed
        zdrasviute: mir
        salud: el mundo
"""


class TestExtraFitsKeywordsInit:
    def test_initialises_with_nothing(self):
        assert isinstance(ExtraFitsKeywords(), ExtraFitsKeywords)

    @pytest.mark.usefixtures("yaml_string")
    def test_initialies_with_yaml_string(self, yaml_string):
        eff = ExtraFitsKeywords(yaml_string=yaml_string)
        assert isinstance(eff, ExtraFitsKeywords)


# @pytest.mark.usefixtures("simplecado_opt")
@pytest.mark.usefixtures("comb_hdul")
class TestExtraFitsKeywordsApplyTo:
    def test_works_if_no_resolve_or_opticaltrain(self, comb_hdul):
        header_dict = {"ext_type": "PrimaryHDU",
                       "unresolved_keywords":
                           {"SIM":
                                {"dark_current": "#dark_current.value"}
                            }
                       }
        eff = ExtraFitsKeywords(header_dict=header_dict)
        hdul = eff.apply_to(comb_hdul)
        hdr = hdul[0].header

        assert hdr["HIERARCH SIM dark_current"] == '#dark_current.value'

    def test_errors_if_resolve_and_no_opticaltrain(self, comb_hdul):
        header_dict = {"ext_type": "PrimaryHDU",
                       "keywords":
                           {"SIM":
                                {"dark_current": "#dark_current.value"}
                            }
                       }
        eff = ExtraFitsKeywords(header_dict=header_dict)
        with pytest.raises(ValueError):
            hdul = eff.apply_to(comb_hdul)

    @pytest.mark.usefixtures("simplecado_opt")
    def test_resolves_hash_strings_with_opticaltrain(self, simplecado_opt,
                                                      comb_hdul):
        header_dict = {"ext_name": "PriHDU",
                       "keywords":
                           {"SIM":
                                {"dark_current": "#dark_current.value",
                                 "telescope_area": "!TEL.area"}
                            }
                       }
        eff = ExtraFitsKeywords(header_dict=header_dict)
        hdul = eff.apply_to(comb_hdul, resolve=True,
                            optical_train=simplecado_opt)

        assert hdul[0].header["HIERARCH SIM dark_current"] == 0.1
        assert hdul[0].header["HIERARCH SIM telescope_area"] == 0.

    @pytest.mark.usefixtures("yaml_string")
    @pytest.mark.usefixtures("simplecado_opt")
    def test_full_yaml_string(self, yaml_string, simplecado_opt, comb_hdul):
        eff = ExtraFitsKeywords(yaml_string=yaml_string)
        hdul = eff.apply_to(comb_hdul, optical_train=simplecado_opt)
        pri_hdr = hdul[0].header

        # resolved keywords
        assert pri_hdr["HIERARCH ESO TEL area"] == 0             # !-str
        assert pri_hdr["HIERARCH ESO TEL pixel_scale"] == 0.004  # !-str
        assert pri_hdr["HIERARCH ESO DAR VALUE"] == 0.1          # #-str
        assert pri_hdr["HIERARCH ESO TEL its_over"] == 9000      # normal
        # unresolved keywords
        assert pri_hdr["HIERARCH ESO ATM TEMPERAT"] == '!ATMO.temperature'
        # comments
        assert pri_hdr.comments["HIERARCH ESO TEL area"] == "[m2] default = 0"
        # ImageHDU header
        assert hdul[1].header["HIERARCH SIM grias_di"] == "woed"


@pytest.mark.usefixtures("comb_hdul")
class TestGetRelevantExtensions:
    def test_works_for_ext_name(self, comb_hdul):
        dic = {"ext_name": "PriHDU"}
        exts = fh.get_relevant_extensions(dic, comb_hdul)
        answer = [0]

        assert np.all([ans in exts for ans in answer])
        assert len(exts) == len(answer)

    def test_works_for_ext_number(self, comb_hdul):
        dic = {"ext_number": [1, 2, 3]}
        exts = fh.get_relevant_extensions(dic, comb_hdul)
        answer = [1, 2]

        assert np.all([ans in exts for ans in answer])
        assert len(exts) == len(answer)

    @pytest.mark.parametrize("ext_type, answer",
                             [("PrimaryHDU", [0]),
                              (["ImageHDU", "PrimaryHDU"], [0, 1])])
    def test_works_for_ext_type(self, comb_hdul, ext_type, answer):
        dic = {"ext_type": ext_type}
        exts = fh.get_relevant_extensions(dic, comb_hdul)

        assert np.all([ans in exts for ans in answer])
        assert len(exts) == len(answer)


class TestFlattenDict:
    def test_works(self):
        dic = {"HIERARCH":
                   {"ESO":
                        {"ATM":
                             {"PWV": 1.0, "AIRMASS": 2.0},
                         "DPR": {"TYPE": "DARK"}},
                    "SIM": {
                        "area": ("!TEL.area", "area")}
                    }
               }
        flat_dict = fh.flatten_dict(dic)
        assert flat_dict["HIERARCH ESO ATM PWV"] == 1.0
        assert flat_dict["HIERARCH ESO ATM AIRMASS"] == 2.0
        assert flat_dict["HIERARCH ESO DPR TYPE"] == "DARK"
        assert flat_dict["HIERARCH SIM area"][0] == "!TEL.area"
        assert flat_dict["HIERARCH SIM area"][1] == "area"

    def test_resolves_bang_strings(self):
        dic = {"SIM": {"random_seed": "!SIM.random.seed"}}
        flat_dict = fh.flatten_dict(dic, resolve=True)
        assert flat_dict["SIM random_seed"] == None

    @pytest.mark.usefixtures("simplecado_opt")
    def test_resolves_hash_strings(self, simplecado_opt):
        dic = {"SIM": {"dark_current": "#dark_current.value"}}
        flat_dict = fh.flatten_dict(dic, resolve=True,
                                    optics_manager=simplecado_opt.optics_manager)
        assert flat_dict["SIM dark_current"] == 0.1


@pytest.mark.usefixtures("yaml_string")
@pytest.mark.usefixtures("simplecado_opt")
class TestEffectsMetaKeywordsApplyTo:
    def test_effect_meta_in_header(self, yaml_string, simplecado_opt, comb_hdul):
        eff = EffectsMetaKeywords()
        hdul = eff.apply_to(comb_hdul, optical_train=simplecado_opt)
        pri_hdr = hdul[0].header

        assert pri_hdr["SIM EFF0 class"] == "DetectorList"
        assert pri_hdr["SIM EFF0 array_dict pixsize"] == "list:[0.015]"

    def test_effect_meta_in_secondary_header(self, yaml_string, simplecado_opt,
                                             comb_hdul):
        eff = EffectsMetaKeywords(ext_number=1, keyword_prefix="HIERARCH GOKU")
        hdul = eff.apply_to(comb_hdul, optical_train=simplecado_opt)
        pri_hdr = hdul[0].header
        sec_hdr = hdul[1].header

        assert "GOKU EFF0 class" not in pri_hdr
        assert sec_hdr["GOKU EFF0 class"] == "DetectorList"
        assert sec_hdr["GOKU EFF0 array_dict pixsize"] == "list:[0.015]"


@pytest.mark.usefixtures("simplecado_opt")
class TestSourceDescriptionFitsKeywordsApplyTo: