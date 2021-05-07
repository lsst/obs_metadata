# This file is part of astro_metadata_translator.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
# See the LICENSE file at the top-level directory of this distribution
# for details of code ownership.
#
# Use of this source code is governed by a 3-clause BSD-style
# license that can be found in the LICENSE file.

__all__ = ("read_test_file", "MetadataAssertHelper")

import astropy.units as u
import os
import pickle
import yaml
from collections import OrderedDict
from astropy.time import Time
import astropy.utils.exceptions
import warnings
from astropy.io.fits import Header

from astro_metadata_translator import ObservationInfo

# PropertyList is optional
try:
    import lsst.daf.base as daf_base
except ImportError:
    daf_base = None


# For YAML >= 5.1 need a different Loader for the constructor
try:
    Loader = yaml.FullLoader
except AttributeError:
    Loader = yaml.Loader


# Define a YAML loader for lsst.daf.base.PropertySet serializations that
# we can use if daf_base is not available.
def pl_constructor(loader, node):
    """Construct an OrderedDict from a YAML file containing a PropertyList."""
    pl = OrderedDict()
    yield pl
    state = loader.construct_sequence(node, deep=True)
    for key, dtype, value, comment in state:
        if dtype == "Double":
            pl[key] = float(value)
        elif dtype == "Int":
            pl[key] = int(value)
        elif dtype == "Bool":
            pl[key] = True if value == "true" else False
        else:
            pl[key] = value


if daf_base is None:
    yaml.add_constructor("lsst.daf.base.PropertyList", pl_constructor, Loader=Loader)


def read_test_file(filename, dir=None):
    """Read the named test file relative to the location of this helper

    Parameters
    ----------
    filename : `str`
        Name of file in the data directory.
    dir : `str`, optional.
        Directory from which to read file. Current directory used if none
        specified.

    Returns
    -------
    header : `dict`-like
        Header read from file.
    """
    if dir is not None and not os.path.isabs(filename):
        filename = os.path.join(dir, filename)
    with open(filename) as fd:
        header = yaml.load(fd, Loader=Loader)
    # Cannot directly check for Mapping because PropertyList is not one
    if not hasattr(header, "items"):
        raise ValueError(f"Contents of YAML file {filename} are not a mapping, they are {type(header)}")
    return header


class MetadataAssertHelper:
    """Class with helpful asserts that can be used for testing metadata
    translations.
    """

    def assertCoordinatesConsistent(self, obsinfo, max_sep=1.0, amdelta=0.01):  # noqa: N802
        """Check that SkyCoord, AltAz, and airmass are self consistent.

        Parameters
        ----------
        obsinfo : `ObservationInfo`
            Object to check.
        max_sep : `float`, optional
            Maximum separation between AltAz derived from RA/Dec headers
            and that found in the AltAz headers.
        amdelta : `float`, optional
            Max difference between header airmass and derived airmass.

        Raises
        ------
        AssertionError
            Inconsistencies found.
        """
        self.assertIsNotNone(obsinfo.tracking_radec)
        self.assertIsNotNone(obsinfo.altaz_begin)

        # Is airmass from header close to airmass from AltAz headers?
        # In most cases there is uncertainty over whether the elevation
        # and airmass in the header are for the start, end, or middle
        # of the observation.  Sometimes the AltAz is from the start
        # but the airmass is from the middle so accuracy is not certain.
        self.assertAlmostEqual(obsinfo.altaz_begin.secz.to_value(), obsinfo.boresight_airmass, delta=amdelta)

        # Is AltAz from headers close to AltAz from RA/Dec headers?
        # Can trigger warnings from Astropy if date is in future
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=astropy.utils.exceptions.AstropyWarning)
            sep = obsinfo.altaz_begin.separation(obsinfo.tracking_radec.altaz)
        self.assertLess(sep.to_value(unit="arcmin"), max_sep, msg="AltAz inconsistent with RA/Dec")

    def assertObservationInfoIsConsistent(self, file, dir=None, check_wcs=True,  # noqa: N802
                                          wcs_params=None, **kwargs):
        """Check contents of an ObservationInfo instantiated from dict-like and
        Astropy Header objects..

        Parameters
        ----------
        file : `str`
            Path to YAML file representing the header.
        dir : `str`, optional
            Optional directory from which to read ``file``.
        check_wcs : `bool`, optional
            Check the consistency of the RA/Dec and AltAz values.
        wcs_params : `dict`, optional
            Parameters to pass to `assertCoordinatesConsistent`.
        kwargs : `dict`
            Keys matching `ObservationInfo` properties with values
            to be tested.

        Raises
        ------
        AssertionError
            A value in the ObservationInfo derived from the file is
            inconsistent.
        """
        header = read_test_file(file, dir=dir)

        try:
            self.assertObservationInfo(header, filename=file, check_wcs=check_wcs,
                                       wcs_params=wcs_params, **kwargs)
        except AssertionError as error:
            raise AssertionError("ObservationInfo derived from an dict-like "
                                 "type is inconsistent.") from error

        # For Header, `in` operator does not work as expected for a dict-like
        # object. Explicitly verify consistency when instantiating from Header
        if not isinstance(header, dict):
            header = header.toDict()
        # Header can not parse multi-line header keys as lists
        for key, val in header.items():
            if isinstance(val, list):
                header[key] = "".join(val)

        header = Header(header)
        try:
            self.assertObservationInfo(header, filename=file, check_wcs=check_wcs,
                                       wcs_params=wcs_params, **kwargs)
        except AssertionError as error:
            raise AssertionError("ObservationInfo derived from an Astropy "
                                 "Header is inconsistent.") from error

    def assertObservationInfo(self, header, filename=None, check_wcs=True,  # noqa: N802
                              wcs_params=None, **kwargs):
        """Check contents of an ObservationInfo.

        Parameters
        ----------
        header : `dict`-like
            Header to be checked.
        filename : `str`, optional
            Name of the filename associated with this header if known.
        check_wcs : `bool`, optional
            Check the consistency of the RA/Dec and AltAz values.  Checks
            are automatically disabled if the translated header does
            not appear to be "science".
        wcs_params : `dict`, optional
            Parameters to pass to `assertCoordinatesConsistent`.
        kwargs : `dict`
            Keys matching `ObservationInfo` properties with values
            to be tested.

        Raises
        ------
        AssertionError
            A value in the ObservationInfo derived from the file is
            inconsistent.
        """
        # For testing we force pedantic mode since we are in charge
        # of all the translations
        obsinfo = ObservationInfo(header, pedantic=True, filename=filename)
        translator = obsinfo.translator_class_name

        # Check that we can pickle and get back the same properties
        newinfo = pickle.loads(pickle.dumps(obsinfo))
        self.assertEqual(obsinfo, newinfo)

        # Check the properties
        for property, expected in kwargs.items():
            calculated = getattr(obsinfo, property)
            msg = f"Comparing property {property} using translator {translator}"
            if isinstance(expected, u.Quantity) and calculated is not None:
                calculated = calculated.to_value(unit=expected.unit)
                expected = expected.to_value()
                self.assertAlmostEqual(calculated, expected, msg=msg)
            elif isinstance(calculated, u.Quantity):
                # Only happens if the test is not a quantity when it should be
                self.fail(f"Expected {expected!r} but got Quantity '{calculated}': {msg}")
            elif isinstance(expected, float) and calculated is not None:
                self.assertAlmostEqual(calculated, expected, msg=msg)
            else:
                self.assertEqual(calculated, expected, msg=msg)

        # Date comparison error reports will benefit by specifying ISO
        # format.  Generate a new Time object at fixed precision
        # to work around the fact that (as of astropy 3.1) adding 0.0 seconds
        # to a Time results in a new Time object that is a few picoseconds in
        # the past.
        def _format_date_for_testing(date):
            if date is not None:
                date.format = "isot"
                date.precision = 9
                date = Time(str(date), scale=date.scale, format="isot")
            return date

        datetime_begin = _format_date_for_testing(obsinfo.datetime_begin)
        datetime_end = _format_date_for_testing(obsinfo.datetime_end)

        # Check that dates are defined and end is the same or after beginning
        self.assertLessEqual(datetime_begin, datetime_end)

        # Check that exposure time is not outside datetime_end
        self.assertLessEqual(obsinfo.datetime_begin + obsinfo.exposure_time,
                             obsinfo.datetime_end)

        # Check the WCS consistency
        if check_wcs and obsinfo.observation_type == "science":
            if wcs_params is None:
                wcs_params = {}
            self.assertCoordinatesConsistent(obsinfo, **wcs_params)
