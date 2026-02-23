"""Application configuration management."""

import copy
import importlib.resources as resources
import logging
import os
import re
import shutil
import subprocess
import sys
from enum import StrEnum
from pathlib import Path
from typing import IO, Any, Dict, TypeVar

import yaml
from ruamel.yaml import YAML as rYAML  # noqa: N811

from tavi import __version__ as tavi_version
from tavi.meta.decorators.singleton import Singleton
from tavi.meta.time import isoFromTimestamp, timestamp

package_name = "tavi"


class DeployEnvEnum(StrEnum):
    """Deployment environment enumeration."""

    NEXT = f"{package_name}_next"
    QA = f"{package_name}_qa"
    PROD = f"{package_name}_prod"


def isTestEnv() -> bool:
    """
    Check if the current environment is a test environment.

    This is determined by the presence of the "env" environment variable
    and whether it contains the string "test".
    """
    env = os.environ.get("env")
    return env and "test" in env and "conftest" in sys.modules


def _find_root_dir() -> str:
    """Find the root directory of the TAVI module."""
    try:
        MODULE_ROOT = Path(sys.modules[package_name].__file__).parent

        # Using `"test" in env` here allows different versions of "[category]_test.yml" to be used for different
        #  test categories: e.g. unit tests use "test.yml" but integration tests use "integration_test.yml".
        if isTestEnv():
            # WARNING: there are now multiple "conftest.py" at various levels in the test hierarchy.
            MODULE_ROOT = MODULE_ROOT.parent.parent / "tests"
    except Exception as e:
        raise RuntimeError(f"Unable to determine {package_name} module-root directory") from e

    return str(MODULE_ROOT)


# update config with self._config, deep_update does not work with ruamel.yaml
def merge_dicts(d1: Dict[str, Any], d2: Dict[str, Any]) -> Dict[str, Any]:
    """Merge two dictionaries with d2 overwriting values in d1."""
    for key, value in d2.items():
        if isinstance(value, dict) and key in d1:
            merge_dicts(d1[key], value)
        else:
            d1[key] = value
    return d1


@Singleton
class _Resource:
    """Resource manager for loading application resources."""

    _package_mode: bool
    _resources_path: str
    _logger = logging.getLogger(__name__ + ".Resource")

    def __init__(self) -> None:
        """Initialize the resource manager."""
        # where the location of resources are depends on whether or not this is in package mode
        self._package_mode = not self._existsInPackage("application.yml")
        if self._package_mode:
            self._logger.debug("In package mode")
            self._resources_path = "/resources/"
        else:
            self._logger.debug("Not in package mode")
            self._resources_path = os.path.join(_find_root_dir(), "resources/")

    def _existsInPackage(self, sub_path: str) -> bool:
        """Check if a resource exists in the package."""
        with resources.path(f"{package_name}.resources", sub_path) as path:
            return os.path.exists(path)

    def exists(self, sub_path: str) -> bool:
        """Check if a resource exists."""
        if self._package_mode:
            return self._existsInPackage(sub_path)
        else:
            return os.path.exists(self.getPath(sub_path))

    def getPath(self, sub_path: str) -> str:
        """Get the path to a resource."""
        if sub_path.startswith("/"):
            return os.path.join(self._resources_path, sub_path[1:])
        else:
            return os.path.join(self._resources_path, sub_path)

    def read(self, sub_path: str) -> str:
        """Read a resource file."""
        with self.open(sub_path, "r") as file:
            return file.read()

    def open(self, sub_path: str, mode: str) -> IO[Any]:
        """Open a resource file."""
        if self._package_mode:
            with resources.path(f"{package_name}.resources", sub_path) as path:
                return open(path, mode)
        else:
            return open(self.getPath(sub_path), mode)


Resource = _Resource()

KeyType = TypeVar("KeyType")


def deep_update(mapping: Dict[KeyType, Any], *updating_mappings: Dict[KeyType, Any]) -> Dict[KeyType, Any]:
    """Crawl through initial dictionary while appending and updating entries in nested dicts."""
    updated_mapping = mapping.copy()
    for updating_mapping in updating_mappings:
        for k, v in updating_mapping.items():
            if k in updated_mapping and isinstance(updated_mapping[k], dict) and isinstance(v, dict):
                updated_mapping[k] = deep_update(updated_mapping[k], v)
            else:
                updated_mapping[k] = v
    return updated_mapping


@Singleton
class _Config:
    _config: Dict[str, Any] = {}
    _logger = logging.getLogger(__name__ + ".Config")
    _property_change_warnings: Dict[str, str] = {}
    _default_env: str = "application.yml"

    def __init__(self) -> None:
        self._property_change_warnings = {
            "version.start": (
                "It is NOT ADVISED to change the `version.start` property"
                " WITHOUT CAREFUL CONSIDERATION of the file indexes on disk."
            ),
        }

        self.reload()
        # if -user.yml exists, then load it by default
        if self.shouldSwapToUserYml():
            self._logger.info("Found user configuration file, swapping to it")
            self.swapToUserYml()

    def shouldSwapToUserYml(self) -> bool:
        """Check if the app should use a user config."""
        isExplicitEnv = "env" in os.environ
        isUserYmlExists = (self._userHome() / f"{package_name}-user.yml").exists()
        return isUserYmlExists and not isExplicitEnv and not isTestEnv()

    def _fix_directory_properties(self) -> None:
        """Expand ~/ to a specified path."""

        def expandhome(direc: str) -> str:
            if "~" in direc:
                return str(Path(direc).expanduser())
            else:
                return direc

        if "instrument" in self._config and "home" in self._config["instrument"]:
            self._config["instrument"]["home"] = expandhome(self._config["instrument"]["home"])
        if "samples" in self._config and "home" in self._config["samples"]:
            self._config["samples"]["home"] = expandhome(self._config["samples"]["home"])

    def configureForDeploy(self) -> None:
        version = self.packageVersion()
        if "dev" in version:
            self.mergeAndExport(DeployEnvEnum.NEXT)
        elif "rc" in version:
            self.mergeAndExport(DeployEnvEnum.QA)
        else:
            self.mergeAndExport(DeployEnvEnum.PROD)

    def mergeAndExport(self, envName: str) -> None:
        """Merge/export the current configuration with the specified environment configuration."""
        self._logger.debug(f"Merging/exporting configuration with environment: {envName}")
        self.refresh(envName, False)

        ryaml = rYAML()
        ryaml.default_flow_style = False
        ryaml.indent(mapping=2, sequence=4, offset=2)

        with Resource.open(self._default_env, "r") as file:
            config = ryaml.load(file)

        # overwrite default env application.yml with exported config
        merge_dicts(config, self._config)

        # Export the merged configuration to application.yml
        with Resource.open(self._default_env, "w") as file:
            ryaml.dump(config, file)

    def reload(self, env_name: str = None) -> None:
        # use refresh to do initial load, clearing shouldn't matter
        self.refresh(self._default_env, True)

        # ---------- internal values: --------------------------
        # allow "resources" relative paths to be entered into the "yml"
        #   using "${module.root}"
        self._config["module"] = {}
        self._config["module"]["root"] = _find_root_dir()

        self._config["version"] = self._config.get("version", {})
        self._config["version"]["default"] = -1

        # ---------- end: internal values: -----------------------------

        watchedProperties = {}
        for key in self._property_change_warnings:
            if self.exists(key):
                watchedProperties[key] = self[key]

        # see if user used environment injection to modify what is needed
        # this will get from the os environment or from the currently loaded one
        # first case wins
        self.env = env_name
        if self.env is None:
            self.env = os.environ.get("env", self._config.get("environment", None))
        if self.env is not None:
            self._logger.info(f"Loading environment config: {self.env}")
            self.refresh(self.env)
        else:
            self._logger.info("No environment config specified, using default")
        self.warnSensitiveProperties(watchedProperties)
        self.persistBackup()

    @property
    def userApplicationDataHome(self) -> Path:
        return Path(self["user.application.home"]).expanduser()

    def persistBackup(self) -> None:
        self.userApplicationDataHome.mkdir(parents=True, exist_ok=True)
        backupFile = self.userApplicationDataHome / "application.yml.bak"
        with open(backupFile, "w") as file:
            yaml.dump(self._config, file, default_flow_style=False)
        self._logger.info(f"Backup of application.yml created at {backupFile.absolute()}")

    def loadEnv(self, env_name: str) -> None:
        # load the new environment
        self.reload(env_name)

    @staticmethod
    def _timestamp() -> str:
        return isoFromTimestamp(timestamp())

    def archiveUserYml(self) -> None:
        """Archive the user config for safe keeping."""
        # check if -user.yml exists
        userHome = self._userHome()
        if (userHome / f"{package_name}-user.yml").exists():
            version = self.getUserYmlVersionDisk()

            # generate human readable timestamp
            timestamp = self._timestamp()

            # archive the old -user.yml
            archivePath = userHome / f"{package_name}-user-{version}-{timestamp}.yml.bak"
            shutil.copy(str(userHome / f"{package_name}-user.yml"), str(archivePath))

    @staticmethod
    def _userHome() -> Path:
        return Path.home() / f".{package_name}"

    def packageVersion(self) -> str | None:
        """Get the version of the application."""
        if tavi_version is None or tavi_version == "unknown" or tavi_version == "":
            label = subprocess.check_output(["git", "rev-parse", "HEAD"]).strip()
            if bool(label) and not label == b"":
                return label.decode("utf-8")
            raise ValueError(
                f"The {package_name} Version is not set correctly. "
                f"Please ensure that the {package_name} package is installed correctly."
            )
        return tavi_version

    def getUserYmlVersionDisk(self) -> str | None:
        """Get the associated app version the user config was generated with."""
        # check if -user.yml exists
        userHome = self._userHome()
        if (userHome / f"{package_name}-user.yml").exists():
            with open(str(userHome / f"{package_name}-user.yml"), "r") as f:
                applicationYml = rYAML(typ="safe").load(f)
            version = applicationYml.get("application", {"version": None})["version"]
            return version
        else:
            return None

    def swapToUserYml(self) -> None:
        """Swap to user generated config."""
        # generate root directory for user configurations
        userHome = self._userHome()

        try:
            if not userHome.exists():
                userHome.mkdir(parents=True, exist_ok=True)

            # check if -user.yml exists
            if self.getUserYmlVersionDisk() != self.packageVersion():
                # if the version is not the same, then we need to archive the old one
                if self.getUserYmlVersionDisk() is not None:
                    self._logger.warning(
                        "The user configuration file is out of date. A new configuration file will be generated."
                    )
                # archive the old -user.yml
                self.archiveUserYml()
                # generate a new valid -user.yml
                self._generateUserYml()

        except Exception as e:  # noqa: BLE001
            raise RuntimeError(
                (
                    "Unable to swap to user configuration: "
                    f"{userHome / f'{package_name}-user.yml'}"
                    "\nOriginal configuration maintained."
                )
            ) from e
        else:
            # load the user yml file if the filesystem is ready
            self.loadEnv(str(userHome / f"{package_name}-user.yml"))
            # archive a backup of the current user yml
            self.archiveUserYml()

    def _generateUserYml(self) -> None:
        """Generate a default application config in the user's home dir."""
        userHome = self._userHome()
        # copy commented out application.yml as -user.yml
        # read the application.yml as a dict from resources
        applicationYml = None
        with Resource.open("application.yml", "r") as f:
            applicationYml = rYAML(typ="safe").load(f)

        if applicationYml.get("application") is None:
            applicationYml["application"] = {}
        applicationYml["application"]["version"] = self.packageVersion()
        applicationYml["environment"] = f"{package_name}-user"

        # convert the dict back to a yaml string
        applicationYmlStr = yaml.dump(applicationYml, default_flow_style=False)

        # write the application.yml to the user home as -user.yml
        with open(userHome / f"{package_name}-user.yml", "w") as f:
            f.write(applicationYmlStr)

    def getCurrentEnv(self) -> str:
        if self.env is not None:
            return self.env
        else:
            # this is the default environment
            return "default"

    def refresh(self, env_name: str, clearPrevious: bool = False) -> None:
        """Load and refresh configuration from environment files."""
        # save a copy of previous config if it fails to load
        prevConfig = copy.deepcopy(self._config)

        try:
            if clearPrevious:
                self._config.clear()

            if env_name.endswith(".yml"):
                # name to be put into config
                new_env_name = env_name

                # this is a filename that should be loaded
                filepath = Path(env_name)
                if filepath.exists():
                    self._logger.debug(f"Loading config from {filepath.absolute()}")
                    with open(filepath, "r") as file:
                        envConfig = rYAML(typ="safe").load(file)
                else:
                    # load from the resource
                    with Resource.open(env_name, "r") as file:
                        envConfig = rYAML(typ="safe").load(file)
                    new_env_name = env_name.replace(".yml", "")
                # update the configuration with the  new environment
                self._config = merge_dicts(self._config, envConfig)

                # add the name to the config object if it wasn't specified
                if "environment" not in envConfig:
                    self._config["environment"] = new_env_name

                # do fixups on items that are directories
                self._fix_directory_properties()
            else:
                # recurse this function with a fuller name
                self.refresh(f"{env_name}.yml", clearPrevious)
        except Exception:
            # if it fails, restore the previous config
            self._logger.warning(f"Failed to load {env_name}.yml, restoring previous config")
            self._config = prevConfig
            raise

    def warnSensitiveProperties(self, watchedProperties: dict[str, Any]) -> None:
        for key in watchedProperties:
            msg = self._property_change_warnings[key]
            if watchedProperties[key] != self[key]:
                warningBar = ("/" * 20) + " WARNING " + ("/" * 20)
                self._logger.warning(warningBar)
                self._logger.warning(f"Property '{key}' was changed in {self.env}.yml")
                self._logger.warning(msg)
                self._logger.warning(warningBar)

    # method to regex for string pattern of ${key} and replace with value
    def _replace(self, value: str, remainingKeys: list[str]) -> str:
        # if the value is not a string, then just return it
        if not isinstance(value, str):
            return value

        # Regex all keys of the form ${key.subkey} and store in a list
        regex = r"\$\{([a-zA-Z0-9_\.]+)\}"
        matches = [match for match in re.finditer(regex, value, re.MULTILINE)]
        # replace all keys with their values
        if len(remainingKeys) == 0:
            for match in matches:
                key = match.group()[2:-1]
                if isinstance(self[key], dict):
                    return value
                value = value.replace(f"${{{key}}}", self[key])
        else:
            match = matches[0]
            rootKey = match.group()[2:-1]
            key = rootKey
            val = self[key]

            # while val is a dict, keep appending keys
            while isinstance(val, dict):
                key = key + f".{remainingKeys.pop(0)}"
                val = self[key]

            value = value.replace(f"${{{rootKey}}}", val)
            # if(len(remainingKeys) > 0):
            value = self._replace(value, remainingKeys)

        return value

    def exists(self, key: str) -> bool:
        val = self._find(key)
        return val is not None

    def _find(self, key: str) -> Any:
        keys = key.split(".")
        val = self._config.get(keys[0])
        totalProcessed = 0
        for k in keys[1:]:
            if val is None:
                break
            if isinstance(val, str):
                break
            totalProcessed += 1
            val = val.get(k)

        if val is not None:
            val = self._replace(val, keys[1 + totalProcessed :])
        return val

    # period delimited key lookup
    def __getitem__(self, key: str) -> Any:
        """Lookup config value via . delimited key."""
        val = self._find(key)
        if val is None:
            raise KeyError(f"Key '{key}' not found in configuration")
        if isinstance(val, dict):
            # we need to solve the keys present in dict members
            # we can accomplish this by recursively calling __getitem__
            return {k: self[f"{key}.{k}"] for k in val}
        return val

    def validate(self) -> None:
        """Validate for conflicting property values."""
        pass


Config = _Config()

# `Logging.Level` is not an `Enum`.
_python_logging_level_map = {
    "notset": logging.NOTSET,
    "critical": logging.CRITICAL,
    "error": logging.ERROR,
    "warning": logging.WARNING,
    "info": logging.INFO,
    "debug": logging.DEBUG,
}

_python_logging_level_map_from_mantid = {
    "none": logging.NOTSET,
    "fatal": logging.CRITICAL + 5,  # exists only as a POCO level
    "critical": logging.CRITICAL,
    "error": logging.ERROR,
    "warning": logging.WARNING,
    "notice": logging.INFO + 5,  # exists only as a POCO level
    "information": logging.INFO,
    "debug": logging.DEBUG,
    "trace": logging.DEBUG - 5,  # exists only as a POCO level
}


def fromMantidLoggingLevel(level: str) -> int:
    """Convert mantid log level to python log level."""
    # Python logging level from Mantid logging level

    # Python levels:
    #   logging.NOTSET: 0, logging.DEBUG: 10, logging.INFO: 20,
    #     logging:WARNING: 30, logging.ERROR: 40, logging.CRITICAL: 50

    # Poco levels (implemented as Poco::Message::Priority enum):
    # 'none': , 'fatal', 'critical', 'error', 'warning', 'notice', 'information', 'debug', 'trace'

    pythonLevel = logging.NOTSET
    try:
        pythonLevel = _python_logging_level_map_from_mantid[level.lower()]
    except KeyError as e:
        raise RuntimeError(f"can't convert '{e}' to a Python logging level")
    return pythonLevel


def fromPythonLoggingLevel(level: int | str) -> str:
    """Convert python log level to mantid log level."""
    # Mantid logging level from Python logging level
    try:
        level = level if isinstance(level, int) else _python_logging_level_map[level.lower()]
        mantidLevel = {v: k for k, v in _python_logging_level_map_from_mantid.items()}[level]
    except KeyError as e:
        raise RuntimeError(f"can't convert '{e}' to a Mantid logging level")
    return mantidLevel
