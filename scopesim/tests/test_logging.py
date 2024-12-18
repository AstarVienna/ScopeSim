# -*- coding: utf-8 -*-
"""Test to make sure the logging configuration is applied."""

import logging
from importlib import reload

import pytest

import scopesim as sim


@pytest.fixture(scope="function")
def reload_scopesim():
    """Temprarily disable the global configure_logging fixture."""
    base_logger = logging.getLogger("astar")
    handlers = base_logger.handlers
    prop = base_logger.propagate
    # Reload scopesim to apply logging configuration again
    reload(sim)
    yield
    # Restore
    base_logger.handlers = handlers
    base_logger.propagate = prop


@pytest.mark.usefixtures("reload_scopesim")
def test_loggers_are_configured():
    base_logger_dict = sim.rc.__logging_config__["loggers"]["astar"]

    base_logger = logging.getLogger("astar")
    sim_logger = base_logger.getChild("scopesim")

    base_logger_level = logging.getLevelName(base_logger.getEffectiveLevel())
    assert base_logger_level == base_logger_dict["level"]

    sim_logger_level = logging.getLevelName(sim_logger.getEffectiveLevel())
    assert sim_logger_level == "DEBUG"

    assert base_logger.propagate == base_logger_dict["propagate"]
    assert sim_logger.propagate

    for handler, name in zip(base_logger.handlers,
                             base_logger_dict["handlers"]):
        handler.name == name


def test_log_to_file():
    base_logger = logging.getLogger("astar")
    sim.log_to_file(enable=True)
    assert any(handler.name == "file" for handler in base_logger.handlers)
    sim.log_to_file(enable=False)
    assert not any(handler.name == "file" for handler in base_logger.handlers)


def test_set_console_log_level():
    base_logger = logging.getLogger("astar")
    sim.set_console_log_level("ERROR")
    assert base_logger.handlers[0].level == logging.ERROR
    sim.set_console_log_level()
    assert base_logger.handlers[0].level == logging.INFO
