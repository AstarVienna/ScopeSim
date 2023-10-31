# Guide for ScopeSim test

This document is intended as a brief summary of pytest's fixture functionality in the context of what is needed for tests in ScopeSim.
It should serve as a guideline on how to design unit- and integration tests to a common standard, so that similar tasks are performed in a standardized way throughout the ScopeSim test suite.
In particular, the use of setup and teardown actions, mock files and patching of globals will be discussed here.

## Fixtures and their use

Pytest provides the ability to use [fixtures](https://docs.pytest.org/en/6.2.x/fixture.html) for managing the environment around a test.
This can come in the form of providing objects used by a test (available to the test like a function parameter) or via any kind of setup and teardown actions.
Fixtures can have different scopes (function, class, module, etc.) which determines when they are freshly executed.
I.e. if a fixture with scope "class" produces an object needed by multiple test in a class, the same instance of that object will be shared by all test in the class.
Thus, if a test might modify that object, it can be better to use "function" scope (the default), to have a clean instance for each test, resulting in a more reliable and repeatable test execution.

There are two principle ways how fixtures can be used by tests: either an argument to a test function (which makes sense for fixtures that `return` or `yield` objects) or via the decorator `@pytest.mark.usefixtures("fixture_name", "another_fixture_name")` around a test function or class.
The latter case makes sense, if the fixture does not return an object that could be used by the test, but rather modifies the environment around the test in some other way.
Using this decorator will ensure the fixture is called anyway, even if its name is not included as an argument to the test function.
Note that using the `usefixtures` decorator on a test function that _does_ use the same fixture as an argument is generally redundant, as it doesn't do anything differently.

The following examples should serve to illustrate these two way of using fixtures:
```python
import some_global

@pytest.fixture(scope="class")
def setup_global():
    some_global.modify()

@pytest.fixture(scope="function")
def provide_object():
    my_object = MyClass()
    return my_object

@pytest.mark.usefixtures("setup_global")
class TestSomething:
    # don't need to decorate "provide_object", as it's used as an argument
    def test_something(self, provide_object):
        assert provide_object.something()
```

## Setup and Teardown

For some test cases it might be required to perform more or less complicated setup and teardown procedures.
This should generally be done with pytest's `yield`-fixture feature, which allows execution of further code in the fixture once the control flow is passed back out of the test using the fixture.

The following example will perform unspecified setup and teardown actions before and after each test class is used, that was decorated with this fixture.
```python
@pytest.fixture(scope="class")
def setup_and_teardown():
    # perform actions required for setup
    yield
    # perform actions required for teardown

@pytest.mark.usefixtures("setup_and_teardown")
class TestSomething:
    def test_something(self):
        assert something
```

This functionality can also be combined with a context manager for setup and teardown, e.g.:
```python
@pytest.fixture(scope="class")
def setup_and_teardown():
    with SetupContext():
        yield
```

The fixture can also `yield` and object if required.
The following example will yield an unspecified object to each function that uses the fixture.
Note that the use of `usefixtures` is not required if the returned/yielded value is used, as discussed above.
```python
@pytest.fixture(scope="function")
def setup_object():
    my_object = MyClass()
    yield my_object
    my_object.teardown()

def test_something_else(setup_object):
    assert something(setup_object)
```

In the case that the setup and teardown is required on the module level (such as a temporary directory), the keywords `scope="module"` and `autouse=True` can be used.
In that case, the fixture doesn't have to be "used" or mentioned anywhere else in that module. It will only be called once per module.
```python
@pytest.fixture(scope="module", autouse=True)
def setup_and_teardown():
    # perform actions required for setup
    yield
    # perform actions required for teardown
```

## Mock files

Whenever a mock data file, such as would ordinarily be located in an IRDB instrument package, is passed to a class or function as part of a test, the path of that file should be given as an absolute path, if possible using the mock path fixtures described below.
ScopeSim offers the possibility to refer to a file only by its name, if it is located in a directory present in `rc.__search_path__`.
While this is a feature intended to be used in normal applications, in the context of testing, it is preferred to make it more explicit where a file should be found.
This is mostly to avoid confusion if more than one file of the same name exists in different subfolders of the root mock directory, but also to reduce the risk of unexpected failures as well as "false positive" passes.
In some cases however, an internal process within ScopeSim might look for a file other than what can be specified in the test. In that case, `rc.__search_path__` needs to be patched to point to the correct location (and only there), see below for details.

### Predefined global mocks

The following mock paths (all as `pathlib.Path` objects) are available as global fixtures:

- `mock_dir`: root directory for all mock files
- `mock_path`: equivalent to `mock_dir / "files"`, location of most generic mock files.
- `mock_path_yamls`: equivalent to `mock_dir / "yamls"`, location of some mock YAML files.
- `mock_path_micado`: equivalent to `mock_dir / "MICADO_SCAO_WIDE"`, location of MICADO-specific mock files. Discouraged for new development, see below.

### On the use of package-specific mocks

Several tests currently use mock files specific to certain IRDB instrument packages, most notably MICADO.
Any new development should ideally use a generic non-instrument-specific test and mock configuration, making use of the `mocks.basic_instrument` as much as possible, and extending it's functionality where needed.

## Patching globals

Sometimes it is needed to simulate certain global configurations for a test case.
In those cases, it is crucial to not just modify that global, otherwise there is a significant risk of "polluting" these settings for other tests!
Instead, use e.g. `unittest.mock.patch` from the standard library to make sure those modifications stay encapsulated.
If the global to be patched implements the `MutableMapping` protocol, it is possible to use `patch.dict` for straightforward dict-style patching.

The following examples illustrate the use of this patching:

_Sidenote: These examples make use of `rc.__currsys__`, which is expected to see substantial change in future versions of ScopeSim. This guide will be updated when that occurs._

For all examples in this section, we will assume the following has been imported:
```python
from unittest.mock import patch
```

For simple cases, the easiest way can be to use `patch.dict` as a decorator:
```python
@patch.dict("scopesim.rc.__currsys__",
            {"!OBS.detector_readout_mode": "fast"})
```             

For more complex cases, or ones where another fixture (e.g. a mock path, see above) needs to be accessed, `patch` or `patch.dict` can be used as a context manager:
```python
patched = {"!DET.width": 4.2,
           "!DET.pixel_size": 0.1}
with patch.dict("scopesim.rc.__currsys__", patched):
    assert something

# use global fixtures mock_path and mock_path_yamls
with patch("scopesim.rc.__search_path__", [mock_path, mock_path_yamls]):
    assert something
```

Patching can also be used in `yield`-fixtures (see above). These fixtures can then be applied via e.g. `@pytest.mark.usefixtures("no_file_error")`.
This kind of fixture can (like any other fixture) also accept other fixtures, such as the default mock file paths.
```python
@pytest.fixture(scope="function")
def no_file_error():
    """Patch currsys to avoid missing file error."""
    patched = {"!SIM.file.error_on_missing_file": False}
    with patch.dict("scopesim.rc.__currsys__", patched):
        yield
```

### Predefined global patches

The following global `yield`-fixtures are available for patching, to be used via `@pytest.mark.usefixtures()`

- `patch_mock_path` and `patch_mock_path_micado`: patches `rc.__search_path__` to `mock_path` and `mock_path_micado`, respectively (see above).
Note that in these cases, _only_ that path is present in the patched `rc.__search_path__`.
- `patch_all_mock_paths`: like `patch_mock_path`, but also patches `"!SIM.file.local_packages_path"` with `mock_dir`, which is needed by some tests.
- `no_file_error`: sets `rc.__currsys__["!SIM.file.error_on_missing_file"] = False`.
Allowing files to not be present in a specific location is a feature of ScopeSim, used e.g. to determine if something needs to be downloaded or looked for at another location.
However, the tests are generally run with this set to `True`, to spot any cases of files missing unintentionally.
If a test needs the "silent missing" functionality, this patch needs to be applied.
- `protect_currsys`: creates a copy of `rc.__currsys__` for the scope of the test, to avoid polluting the global one.
Should be used around everything the uses `OpticalTrain`.

## Misc

Other points to keep in mind:

- Tests should be able to run as a part of the whole test suite (obviously), but also in a standalone way on a per-module basis.
When adding or modifiying tests, make sure to try them in both settings.
- Tests should be limited in their memory consumption to facilitate running in GitHub Actions.
To still include memory-heavy tests, decorate them with `pytest.mark.skip` (and a meaningful `reason`) before pushing to GitHub.
In the simplest way, comment out that line for local testing with more memory.
