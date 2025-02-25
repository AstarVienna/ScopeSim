# pylint: disable=missing-class-docstring
# pylint: disable=missing-function-docstring

import pytest
from unittest.mock import patch

from pathlib import Path
from tempfile import TemporaryDirectory
from astropy.io import fits
from astropy import units as u
import numpy as np

from scopesim.effects import fits_headers as fh
from scopesim.source.source_templates import star
import scopesim as sim


@pytest.fixture(name="simplecado_opt", scope="function")
def fixture_simplecado_opt(mock_path_yamls):
    with patch("scopesim.rc.__currsys__"):
        simplecado_yaml = str(mock_path_yamls / "SimpleCADO.yaml")
        cmd = sim.UserCommands(yamls=[simplecado_yaml])
        yield sim.OpticalTrain(cmd)


@pytest.fixture(name="comb_hdul", scope="function")
def fixture_comb_hdul():
    pri = fits.PrimaryHDU(header=fits.Header({"EXTNAME": "PriHDU"}))
    im = fits.ImageHDU(header=fits.Header({"EXTNAME": "ImHDU"}))
    tbl = fits.BinTableHDU(header=fits.Header({"EXTNAME": "BinTblHDU"}))
    hdul = fits.HDUList([pri, im, tbl])

    return hdul


# taken from the example yaml in docstring of ExtraFitsKeywords
@pytest.fixture(name="yaml_string", scope="function")
def fixture_yaml_string():
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
    EXTNAME: "DET++.DATA"
"""


class TestExtraFitsKeywordsInit:
    def test_initialises_with_nothing(self):
        assert isinstance(fh.ExtraFitsKeywords(), fh.ExtraFitsKeywords)

    def test_initialies_with_yaml_string(self, yaml_string):
        eff = fh.ExtraFitsKeywords(yaml_string=yaml_string)
        assert isinstance(eff, fh.ExtraFitsKeywords)


class TestExtraFitsKeywordsApplyTo:
    def test_works_if_no_resolve_or_opticaltrain(self, comb_hdul):
        header_dict = {"ext_type": "PrimaryHDU",
                       "unresolved_keywords":
                           {"SIM":
                            {"dark_current": "#dark_current.value"}
                            }
                       }
        eff = fh.ExtraFitsKeywords(header_dict=header_dict)
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
        eff = fh.ExtraFitsKeywords(header_dict=header_dict)
        with pytest.raises(ValueError):
            _ = eff.apply_to(comb_hdul)

    def test_resolves_hash_strings_with_opticaltrain(self, simplecado_opt,
                                                     comb_hdul):
        header_dict = {"ext_name": "PriHDU",
                       "keywords":
                           {"SIM":
                            {"dark_current": "#dark_current.value",
                             "telescope_area": "!TEL.area"}
                            }
                       }
        eff = fh.ExtraFitsKeywords(header_dict=header_dict)
        hdul = eff.apply_to(comb_hdul, resolve=True,
                            optical_train=simplecado_opt)

        assert hdul[0].header["HIERARCH SIM dark_current"] == 0.1
        assert hdul[0].header["HIERARCH SIM telescope_area"] == 1.

    def test_full_yaml_string(self, yaml_string, simplecado_opt, comb_hdul):
        eff = fh.ExtraFitsKeywords(yaml_string=yaml_string)
        hdul = eff.apply_to(comb_hdul, optical_train=simplecado_opt)
        pri_hdr = hdul[0].header

        # resolved keywords
        assert pri_hdr["HIERARCH ESO TEL area"] == 1             # !-str
        assert pri_hdr["HIERARCH ESO TEL pixel_scale"] == 0.004  # !-str
        assert pri_hdr["HIERARCH ESO DAR VALUE"] == 0.1          # #-str
        assert pri_hdr["HIERARCH ESO TEL its_over"] == 9000      # normal
        # unresolved keywords
        assert pri_hdr["HIERARCH ESO ATM TEMPERAT"] == '!ATMO.temperature'
        # comments
        assert pri_hdr.comments["HIERARCH ESO TEL area"] == "[m2] default = 0"
        # ImageHDU header
        assert hdul[1].header["HIERARCH SIM grias_di"] == "woed"
        assert hdul[1].header["EXTNAME"] == "DET1.DATA"


class TestGetRelevantExtensions:
    def test_works_for_ext_name(self, comb_hdul):
        dic = {"ext_name": "PriHDU"}
        exts = fh.get_relevant_extensions(dic, comb_hdul)
        answer = [0]

        assert all(ans in exts for ans in answer)
        assert len(exts) == len(answer)

    def test_works_for_ext_number(self, comb_hdul):
        dic = {"ext_number": [1, 2, 3]}
        exts = fh.get_relevant_extensions(dic, comb_hdul)
        answer = [1, 2]

        assert all(ans in exts for ans in answer)
        assert len(exts) == len(answer)

    @pytest.mark.parametrize("ext_type, answer",
                             [("PrimaryHDU", [0]),
                              (["ImageHDU", "PrimaryHDU"], [0, 1])])
    def test_works_for_ext_type(self, comb_hdul, ext_type, answer):
        dic = {"ext_type": ext_type}
        exts = fh.get_relevant_extensions(dic, comb_hdul)

        assert all(ans in exts for ans in answer)
        assert len(exts) == len(answer)


class TestFlattenDict:
    def test_works(self):
        dic = {
            "HIERARCH": {
                "ESO": {
                    "ATM": {
                        "PWV": 1.0,
                        "AIRMASS": 2.0,
                    },
                    "DPR": {
                        "TYPE": "DARK",
                    }
                },
                "SIM": {
                    "area": ("!TEL.area", "area"),
                    "SRC0": {"scaling_unit": u.mag},
                },
            },
        }
        flat_dict = fh.flatten_dict(dic)
        assert flat_dict["HIERARCH ESO ATM PWV"] == 1.0
        assert flat_dict["HIERARCH ESO ATM AIRMASS"] == 2.0
        assert flat_dict["HIERARCH ESO DPR TYPE"] == "DARK"
        assert flat_dict["HIERARCH SIM area"][0] == "!TEL.area"
        assert flat_dict["HIERARCH SIM area"][1] == "area"
        assert flat_dict["HIERARCH SIM SRC0 scaling_unit"] == "mag"

    def test_resolves_bang_strings(self):
        # TODO: Use fixtures, because success depends on order of tests.
        dic = {"SIM": {"random_seed": "!SIM.random.seed"}}
        flat_dict = fh.flatten_dict(dic, resolve=True)
        assert flat_dict["SIM random_seed"] is None

    def test_resolves_hash_strings(self, simplecado_opt):
        dic = {"SIM": {"dark_current": "#dark_current.value"}}
        flat_dict = fh.flatten_dict(
            dic, resolve=True, optics_manager=simplecado_opt.optics_manager)
        assert flat_dict["SIM dark_current"] == 0.1


class TestEffectsMetaKeywordsApplyTo:
    def test_effect_meta_in_header(self, yaml_string, simplecado_opt,
                                   comb_hdul):
        eff = fh.EffectsMetaKeywords()
        hdul = eff.apply_to(comb_hdul, optical_train=simplecado_opt)
        pri_hdr = hdul[0].header

        assert pri_hdr["SIM EFF0 class"] == "DetectorList"
        assert pri_hdr["SIM EFF0 array_dict pixsize"] == "list:[0.015]"

    def test_effect_meta_in_secondary_header(self, yaml_string, simplecado_opt,
                                             comb_hdul):
        eff = fh.EffectsMetaKeywords(ext_number=1,
                                     keyword_prefix="HIERARCH GOKU")
        hdul = eff.apply_to(comb_hdul, optical_train=simplecado_opt)
        pri_hdr = hdul[0].header
        sec_hdr = hdul[1].header

        assert "GOKU EFF0 class" not in pri_hdr
        assert sec_hdr["GOKU EFF0 class"] == "DetectorList"
        assert sec_hdr["GOKU EFF0 array_dict pixsize"] == "list:[0.015]"


class TestSourceDescriptionFitsKeywordsApplyTo:
    def test_source_meta_keywords_in_header(self, simplecado_opt, comb_hdul):
        star1, star2 = star(flux=1*u.ABmag), star()
        star1.meta["hello"] = "world"
        star2.meta["servus"] = "oida"

        simplecado_opt._last_source = star1 + star2
        eff = fh.SourceDescriptionFitsKeywords()
        hdul = eff.apply_to(comb_hdul, optical_train=simplecado_opt)
        pri_hdr = hdul[0].header

        assert pri_hdr["SIM SRC0 hello"] == "world"
        assert pri_hdr["SIM SRC1 servus"] == "oida"

    def test_value_is_longer_than_80_characters(self, simplecado_opt,
                                                comb_hdul):
        star1, star2 = star(flux=1*u.ABmag), star()
        star1.meta["function_call"] *= 5

        simplecado_opt._last_source = star1 + star2
        simplecado_opt._last_source
        eff = fh.SourceDescriptionFitsKeywords()
        hdul = eff.apply_to(comb_hdul, optical_train=simplecado_opt)
        pri_hdr = hdul[0].header

        assert len(pri_hdr["FNSRC0"]) == 120

        # save to disk, what happens to cards longer than 80 characters
        with TemporaryDirectory() as tmpdir:
            fname = Path(tmpdir, "test.fits")
            hdul.writeto(fname)
            tmp_hdr = fits.getheader(fname)

        assert len(tmp_hdr["FNSRC0"]) == 120


class TestSimulationConfigFitsKeywordsApplyTo:
    def test_sys_dict_dicts_are_added_to_header(self, simplecado_opt,
                                                comb_hdul):
        eff = fh.SimulationConfigFitsKeywords()
        hdul = eff.apply_to(comb_hdul, optical_train=simplecado_opt)
        pri_hdr = hdul[0].header

        assert pri_hdr["SIM CONFIG DET ndit"] == 1
        assert pri_hdr["SIM CONFIG SIM random seed"] is None

    def test_bang_string_untouched_for_resolve_false(self, simplecado_opt,
                                                     comb_hdul):
        eff = fh.SimulationConfigFitsKeywords(resolve=False)
        hdul = eff.apply_to(comb_hdul, optical_train=simplecado_opt)
        pri_hdr = hdul[0].header

        assert pri_hdr["SIM CONFIG DET ndit"] == "!OBS.ndit"


class TestAllFitsKeywordEffects:
    def test_works(self, simplecado_opt):
        hdr_dic = {"ext_type": "PrimaryHDU",
                   "keywords": {
                       "HIERARCH": {
                           "ESO": {
                               "INS": {
                                   "pixel_scale": "!INST.pixel_scale",
                                   "dark_current": "#dark_current.value"
                               }
                           },
                           "SIM":
                               {"hello": "world"}
                       }
                   }
                   }
        cmds = simplecado_opt.cmds
        extra_keys = fh.ExtraFitsKeywords(header_dict=hdr_dic, cmds=cmds)
        opt_keys = fh.EffectsMetaKeywords(cmds=cmds)
        src_keys = fh.SourceDescriptionFitsKeywords(cmds=cmds)
        config_keys = fh.SimulationConfigFitsKeywords(cmds=cmds)

        for eff in [src_keys, opt_keys, config_keys, extra_keys]:
            simplecado_opt.optics_manager.add_effect(eff)

        star1 = star(flux=1 * u.ABmag)
        simplecado_opt.observe(star1)
        hdul = simplecado_opt.readout()[0]

        hdr = hdul[0].header
        assert hdr["ESO INS dark_current"] == 0.1
        assert hdr["ESO INS pixel_scale"] == 0.004
        assert hdr["SIM EFF0 class"] == "DetectorList"
        assert hdr["SIM SRC0 class"] == "Table"
        assert hdr["SIM SRC0 photometric_system"] == "ab"
        assert hdr["SIM CONFIG SIM random seed"] is None
