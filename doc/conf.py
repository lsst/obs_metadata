"""Sphinx configuration file for an LSST stack package.

This configuration only affects single-package Sphinx documenation builds.
"""

from documenteer.sphinxconfig.stackconf import build_package_configs
import lsst.obs.metadata
import lsst.obs.metadata.version


_g = globals()
_g.update(build_package_configs(
    project_name="obs_metadata",
    version=lsst.obs.metadata.version.__version__))
