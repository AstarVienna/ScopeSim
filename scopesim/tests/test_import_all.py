"""Test whether all Python files can be imported."""
import importlib.util
import pkgutil

import scopesim


def import_submodules(package):
    """ Import all submodules of a module, recursively, including subpackages

    :param package: package (name or actual module)
    :type package: str | module
    :rtype: dict[str, types.ModuleType]

    https://stackoverflow.com/a/25562415
    """
    if isinstance(package, str):
        package = importlib.import_module(package)
    results = {}
    for loader, name, is_pkg in pkgutil.walk_packages(package.__path__):
        full_name = package.__name__ + '.' + name
        results[full_name] = importlib.import_module(full_name)
        if is_pkg:
            results.update(import_submodules(full_name))
    return results


def test_import_all():
    """Test whether all Python files can be imported.

    This ensures that all files are included in the code coverage.
    """
    import_submodules(scopesim)
