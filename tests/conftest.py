

import os
import pytest

from tavi.meta.decorators.singleton import reset_Singletons


if not os.environ.get("env"):
    os.environ["env"] = "test"

@pytest.fixture(autouse=True)
def _reset_Singletons(request):
    if not "integration" in request.keywords:
        reset_Singletons()
    yield

@pytest.fixture(scope="class", autouse=True)
def _reset_class_scope_Singletons(request):
    if not "integration" in request.keywords:
        reset_Singletons()
    yield

@pytest.fixture(scope="module", autouse=True)
def _reset_module_scope_Singletons(request):
    if not "integration" in request.keywords:
        reset_Singletons()
    yield



from util.Config_helpers import Config_override_fixture