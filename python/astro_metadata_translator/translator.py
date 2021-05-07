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

"""Classes and support code for metadata translation"""

__all__ = ("MetadataTranslator", "StubTranslator", "cache_translation")

from abc import abstractmethod
import inspect
import importlib
import logging
import warnings
import math

import astropy.units as u
import astropy.io.fits.card
from astropy.coordinates import Angle

from .properties import PROPERTIES

log = logging.getLogger(__name__)

# Location of the root of the corrections resource files
CORRECTIONS_RESOURCE_ROOT = "corrections"

"""Cache of version strings indexed by class."""
_VERSION_CACHE = dict()


def cache_translation(func, method=None):
    """Decorator to cache the result of a translation method.

    Especially useful when a translation uses many other translation
    methods.  Should be used only on ``to_x()`` methods.

    Parameters
    ----------
    func : `function`
        Translation method to cache.
    method : `str`, optional
        Name of the translation method to cache.  Not needed if the decorator
        is used around a normal method, but necessary when the decorator is
        being used in a metaclass.

    Returns
    -------
    wrapped : `function`
        Method wrapped by the caching function.
    """
    name = func.__name__ if method is None else method

    def func_wrapper(self):
        if name not in self._translation_cache:
            self._translation_cache[name] = func(self)
        return self._translation_cache[name]
    func_wrapper.__doc__ = func.__doc__
    func_wrapper.__name__ = f"{name}_cached"
    return func_wrapper


class MetadataTranslator:
    """Per-instrument metadata translation support

    Parameters
    ----------
    header : `dict`-like
        Representation of an instrument header that can be manipulated
        as if it was a `dict`.
    filename : `str`, optional
        Name of the file whose header is being translated.  For some
        datasets with missing header information this can sometimes
        allow for some fixups in translations.
    """

    # These are all deliberately empty in the base class.
    default_search_path = None
    """Default search path to use to locate header correction files."""

    default_resource_package = __name__.split(".")[0]
    """Module name to use to locate the correction resources."""

    default_resource_root = None
    """Default package resource path root to use to locate header correction
    files within the ``default_resource_package`` package."""

    _trivial_map = {}
    """Dict of one-to-one mappings for header translation from standard
    property to corresponding keyword."""

    _const_map = {}
    """Dict defining a constant for specified standard properties."""

    translators = dict()
    """All registered metadata translation classes."""

    supported_instrument = None
    """Name of instrument understood by this translation class."""

    @classmethod
    def defined_in_this_class(cls, name):
        """Report if the specified class attribute is defined specifically in
        this class.

        Parameters
        ----------
        name : `str`
            Name of the attribute to test.

        Returns
        -------
        in_class : `bool`
            `True` if there is a attribute of that name defined in this
            specific subclass.
            `False` if the method is not defined in this specific subclass
            but is defined in a parent class.
            Returns `None` if the attribute is not defined anywhere
            in the class hierarchy (which can happen if translators have
            typos in their mapping tables).

        Notes
        -----
        Retrieves the attribute associated with the given name.
        Then looks in all the parent classes to determine whether that
        attribute comes from a parent class or from the current class.
        Attributes are compared using `id()`.
        """
        # The attribute to compare.
        if not hasattr(cls, name):
            return None
        attr_id = id(getattr(cls, name))

        # Get all the classes in the hierarchy
        mro = list(inspect.getmro(cls))

        # Remove the first entry from the list since that will be the
        # current class
        mro.pop(0)

        for parent in mro:
            # Some attributes may only exist in subclasses. Skip base classes
            # that are missing the attribute (such as object).
            if hasattr(parent, name):
                if id(getattr(parent, name)) == attr_id:
                    return False
        return True

    @staticmethod
    def _make_const_mapping(property_key, constant):
        """Make a translator method that returns a constant value.

        Parameters
        ----------
        property_key : `str`
            Name of the property to be calculated (for the docstring).
        constant : `str` or `numbers.Number`
            Value to return for this translator.

        Returns
        -------
        f : `function`
            Function returning the constant.
        """
        def constant_translator(self):
            return constant

        if property_key in PROPERTIES:
            property_doc, return_type = PROPERTIES[property_key][:2]
        else:
            return_type = type(constant).__name__
            property_doc = f"Returns constant value for '{property_key}' property"

        constant_translator.__doc__ = f"""{property_doc}

        Returns
        -------
        translation : `{return_type}`
            Translated property.
        """
        return constant_translator

    @staticmethod
    def _make_trivial_mapping(property_key, header_key, default=None, minimum=None, maximum=None,
                              unit=None, checker=None):
        """Make a translator method returning a header value.

        The header value can be converted to a `~astropy.units.Quantity`
        if desired, and can also have its value validated.

        See `MetadataTranslator.validate_value()` for details on the use
        of default parameters.

        Parameters
        ----------
        property_key : `str`
            Name of the translator to be constructed (for the docstring).
        header_key : `str` or `list` of `str`
            Name of the key to look up in the header. If a `list` each
            key will be tested in turn until one matches.  This can deal with
            header styles that evolve over time.
        default : `numbers.Number` or `astropy.units.Quantity`, `str`, optional
            If not `None`, default value to be used if the parameter read from
            the header is not defined or if the header is missing.
        minimum : `numbers.Number` or `astropy.units.Quantity`, optional
            If not `None`, and if ``default`` is not `None`, minimum value
            acceptable for this parameter.
        maximum : `numbers.Number` or `astropy.units.Quantity`, optional
            If not `None`, and if ``default`` is not `None`, maximum value
            acceptable for this parameter.
        unit : `astropy.units.Unit`, optional
            If not `None`, the value read from the header will be converted
            to a `~astropy.units.Quantity`.  Only supported for numeric values.
        checker : `function`, optional
            Callback function to be used by the translator method in case the
            keyword is not present.  Function will be executed as if it is
            a method of the translator class.  Running without raising an
            exception will allow the default to be used. Should usually raise
            `KeyError`.

        Returns
        -------
        t : `function`
            Function implementing a translator with the specified
            parameters.
        """
        if property_key in PROPERTIES:
            property_doc, return_type = PROPERTIES[property_key][:2]
        else:
            return_type = "str` or `numbers.Number"
            property_doc = f"Map '{header_key}' header keyword to '{property_key}' property"

        def trivial_translator(self):
            if unit is not None:
                q = self.quantity_from_card(header_key, unit,
                                            default=default, minimum=minimum, maximum=maximum,
                                            checker=checker)
                # Convert to Angle if this quantity is an angle
                if return_type == "astropy.coordinates.Angle":
                    q = Angle(q)
                return q

            keywords = header_key if isinstance(header_key, list) else [header_key]
            for key in keywords:
                if self.is_key_ok(key):
                    value = self._header[key]
                    if default is not None and not isinstance(value, str):
                        value = self.validate_value(value, default, minimum=minimum, maximum=maximum)
                    self._used_these_cards(key)
                    break
            else:
                # No keywords found, use default, checking first, or raise
                # A None default is only allowed if a checker is provided.
                if checker is not None:
                    try:
                        checker(self)
                        return default
                    except Exception:
                        raise KeyError(f"Could not find {keywords} in header")
                    value = default
                elif default is not None:
                    value = default
                else:
                    raise KeyError(f"Could not find {keywords} in header")

            # If we know this is meant to be a string, force to a string.
            # Sometimes headers represent items as integers which generically
            # we want as strings (eg OBSID).  Sometimes also floats are
            # written as "NaN" strings.
            casts = {"str": str, "float": float, "int": int}
            if return_type in casts and not isinstance(value, casts[return_type]) and value is not None:
                value = casts[return_type](value)

            return value

        # Docstring inheritance means it is confusing to specify here
        # exactly which header value is being used.
        trivial_translator.__doc__ = f"""{property_doc}

        Returns
        -------
        translation : `{return_type}`
            Translated value derived from the header.
        """
        return trivial_translator

    @classmethod
    def __init_subclass__(cls, **kwargs):
        """Register all subclasses with the base class and create dynamic
        translator methods.

        The method provides two facilities.  Firstly, every subclass
        of `MetadataTranslator` that includes a ``name`` class property is
        registered as a translator class that could be selected when automatic
        header translation is attempted.  Only name translator subclasses that
        correspond to a complete instrument.  Translation classes providing
        generic translation support for multiple instrument translators should
        not be named.

        The second feature of this method is to convert simple translations
        to full translator methods.  Sometimes a translation is fixed (for
        example a specific instrument name should be used) and rather than
        provide a full ``to_property()`` translation method the mapping can be
        defined in a class variable named ``_constMap``.  Similarly, for
        one-to-one trivial mappings from a header to a property,
        ``_trivialMap`` can be defined.  Trivial mappings are a dict mapping a
        generic property to either a header keyword, or a tuple consisting of
        the header keyword and a dict containing key value pairs suitable for
        the `MetadataTranslator.quantity_from_card()` method.
        """
        super().__init_subclass__(**kwargs)

        # Only register classes with declared names
        if hasattr(cls, "name") and cls.name is not None:
            if cls.name in MetadataTranslator.translators:
                log.warning("%s: Replacing %s translator with %s",
                            cls.name, MetadataTranslator.translators[cls.name], cls)
            MetadataTranslator.translators[cls.name] = cls

        # Check that we have not inherited constant/trivial mappings from
        # parent class that we have already applied. Empty maps are always
        # assumed okay
        const_map = cls._const_map if cls._const_map and cls.defined_in_this_class("_const_map") else {}
        trivial_map = cls._trivial_map \
            if cls._trivial_map and cls.defined_in_this_class("_trivial_map") else {}

        # Check for shadowing
        trivials = set(trivial_map.keys())
        constants = set(const_map.keys())
        both = trivials & constants
        if both:
            log.warning("%s: defined in both const_map and trivial_map: %s",
                        cls.__name__, ", ".join(both))

        all = trivials | constants
        for name in all:
            if cls.defined_in_this_class(f"to_{name}"):
                # Must be one of trivial or constant. If in both then constant
                # overrides trivial.
                location = "by _trivial_map"
                if name in constants:
                    location = "by _const_map"
                log.warning("%s: %s is defined explicitly but will be replaced %s",
                            cls.__name__, name, location)

        # Go through the trival mappings for this class and create
        # corresponding translator methods
        for property_key, header_key in trivial_map.items():
            kwargs = {}
            if type(header_key) == tuple:
                kwargs = header_key[1]
                header_key = header_key[0]
            translator = cls._make_trivial_mapping(property_key, header_key, **kwargs)
            method = f"to_{property_key}"
            translator.__name__ = f"{method}_trivial_in_{cls.__name__}"
            setattr(cls, method, cache_translation(translator, method=method))
            if property_key not in PROPERTIES:
                log.warning(f"Unexpected trivial translator for '{property_key}' defined in {cls}")

        # Go through the constant mappings for this class and create
        # corresponding translator methods
        for property_key, constant in const_map.items():
            translator = cls._make_const_mapping(property_key, constant)
            method = f"to_{property_key}"
            translator.__name__ = f"{method}_constant_in_{cls.__name__}"
            setattr(cls, method, translator)
            if property_key not in PROPERTIES:
                log.warning(f"Unexpected constant translator for '{property_key}' defined in {cls}")

    def __init__(self, header, filename=None):
        self._header = header
        self.filename = filename
        self._used_cards = set()

        # Prefix to use for warnings about failed translations
        self._log_prefix_cache = None

        # Cache assumes header is read-only once stored in object
        self._translation_cache = {}

    @classmethod
    @abstractmethod
    def can_translate(cls, header, filename=None):
        """Indicate whether this translation class can translate the
        supplied header.

        Parameters
        ----------
        header : `dict`-like
            Header to convert to standardized form.
        filename : `str`, optional
            Name of file being translated.

        Returns
        -------
        can : `bool`
            `True` if the header is recognized by this class. `False`
            otherwise.
        """
        raise NotImplementedError()

    @classmethod
    def can_translate_with_options(cls, header, options, filename=None):
        """Helper method for `can_translate` allowing options.

        Parameters
        ----------
        header : `dict`-like
            Header to convert to standardized form.
        options : `dict`
            Headers to try to determine whether this header can
            be translated by this class.  If a card is found it will
            be compared with the expected value and will return that
            comparison.  Each card will be tried in turn until one is
            found.
        filename : `str`, optional
            Name of file being translated.

        Returns
        -------
        can : `bool`
            `True` if the header is recognized by this class. `False`
            otherwise.

        Notes
        -----
        Intended to be used from within `can_translate` implementations
        for specific translators.  Is not intended to be called directly
        from `determine_translator`.
        """
        for card, value in options.items():
            if card in header:
                return header[card] == value
        return False

    @classmethod
    def determine_translator(cls, header, filename=None):
        """Determine a translation class by examining the header

        Parameters
        ----------
        header : `dict`-like
            Representation of a header.
        filename : `str`, optional
            Name of file being translated.

        Returns
        -------
        translator : `MetadataTranslator`
            Translation class that knows how to extract metadata from
            the supplied header.

        Raises
        ------
        ValueError
            None of the registered translation classes understood the supplied
            header.
        """
        file_msg = ""
        if filename is not None:
            file_msg = f" from {filename}"
        for name, trans in cls.translators.items():
            if trans.can_translate(header, filename=filename):
                log.debug("Using translation class %s%s", name, file_msg)
                return trans
        else:
            raise ValueError(f"None of the registered translation classes {list(cls.translators.keys())}"
                             f" understood this header{file_msg}")

    @classmethod
    def translator_version(cls):
        """Return the version string for this translator class.

        Returns
        -------
        version : `str`
            String identifying the version of this translator.

        Notes
        -----
        Assumes that the version is available from the ``__version__``
        variable in the parent module. If this is not the case a translator
        should subclass this method.
        """
        if cls in _VERSION_CACHE:
            return _VERSION_CACHE[cls]

        version = "unknown"
        module_name = cls.__module__
        components = module_name.split(".")
        while components:
            # This class has already been imported so importing it
            # should work.
            module = importlib.import_module(".".join(components))
            if hasattr(module, v := "__version__"):
                version = getattr(module, v)
                if version == "unknown":
                    # LSST software will have a fingerprint
                    if hasattr(module, v := "__fingerprint__"):
                        version = getattr(module, v)
                break
            else:
                # Remove last component from module name and try again
                components.pop()

        _VERSION_CACHE[cls] = version
        return version

    @classmethod
    def fix_header(cls, header, instrument, obsid, filename=None):
        """Apply global fixes to a supplied header.

        Parameters
        ----------
        header : `dict`
            The header to correct. Correction is in place.
        instrument : `str`
            The name of the instrument.
        obsid : `str`
            Unique observation identifier associated with this header.
            Will always be provided.
        filename : `str`, optional
            Filename associated with this header. May not be set since headers
            can be fixed independently of any filename being known.

        Returns
        -------
        modified : `bool`
            `True` if a correction was applied.

        Notes
        -----
        This method is intended to support major discrepancies in headers
        such as:

        * Periods of time where headers are known to be incorrect in some
          way that can be fixed either by deriving the correct value from
          the existing value or understanding the that correction is static
          for the given time.  This requires that the date header is
          known.
        * The presence of a certain value is always wrong and should be
          corrected with a new static value regardless of date.

        It is assumed that one off problems with headers have been applied
        before this method is called using the per-obsid correction system.

        Usually called from `astro_metadata_translator.fix_header`.

        For log messages, do not assume that the filename will be present.
        Always write log messages to fall back on using the ``obsid`` if
        ``filename`` is `None`.
        """
        return False

    @staticmethod
    def _construct_log_prefix(obsid, filename=None):
        """Construct a log prefix string from the obsid and filename.

        Parameters
        ----------
        obsid : `str`
            The observation identifier.
        filename : `str`, optional
            The filename associated with the header being translated.
            Can be `None`.
        """
        if filename:
            return f"{filename}({obsid})"
        return obsid

    @property
    def _log_prefix(self):
        """Standard prefix that can be used for log messages to report
        useful context.

        Will be either the filename and obsid, or just the obsid depending
        on whether a filename is known.

        Returns
        -------
        prefix : `str`
            The prefix to use.
        """
        if self._log_prefix_cache is None:
            # Protect against the unfortunate event of the obsid failing to
            # be calculated. This should be rare but should not prevent a log
            # message from appearing.
            try:
                obsid = self.to_observation_id()
            except Exception:
                obsid = "unknown_obsid"
            self._log_prefix_cache = self._construct_log_prefix(obsid, self.filename)
        return self._log_prefix_cache

    def _used_these_cards(self, *args):
        """Indicate that the supplied cards have been used for translation.

        Parameters
        ----------
        args : sequence of `str`
            Keywords used to process a translation.
        """
        self._used_cards.update(set(args))

    def cards_used(self):
        """Cards used during metadata extraction.

        Returns
        -------
        used : `frozenset` of `str`
            Cards used when extracting metadata.
        """
        return frozenset(self._used_cards)

    @staticmethod
    def validate_value(value, default, minimum=None, maximum=None):
        """Validate the supplied value, returning a new value if out of range

        Parameters
        ----------
        value : `float`
            Value to be validated.
        default : `float`
            Default value to use if supplied value is invalid or out of range.
            Assumed to be in the same units as the value expected in the
            header.
        minimum : `float`
            Minimum possible valid value, optional.  If the calculated value
            is below this value, the default value will be used.
        maximum : `float`
            Maximum possible valid value, optional.  If the calculated value
            is above this value, the default value will be used.

        Returns
        -------
        value : `float`
            Either the supplied value, or a default value.
        """
        if value is None or math.isnan(value):
            value = default
        else:
            if minimum is not None and value < minimum:
                value = default
            elif maximum is not None and value > maximum:
                value = default
        return value

    @staticmethod
    def is_keyword_defined(header, keyword):
        """Return `True` if the value associated with the named keyword is
        present in the supplied header and defined.

        Parameters
        ----------
        header : `dict`-lik
            Header to use as reference.
        keyword : `str`
            Keyword to check against header.

        Returns
        -------
        is_defined : `bool`
            `True` if the header is present and not-`None`. `False` otherwise.
        """
        if keyword is None or keyword not in header:
            return False

        if header[keyword] is None:
            return False

        # Special case Astropy undefined value
        if isinstance(header[keyword], astropy.io.fits.card.Undefined):
            return False

        return True

    def resource_root(self):
        """Package resource to use to locate correction resources within an
        installed package.

        Returns
        -------
        resource_package : `str`
            Package resource name.  `None` if no package resource are to be
            used.
        resource_root : `str`
            The name of the resource root.  `None` if no package resources
            are to be used.
        """
        return (self.default_resource_package, self.default_resource_root)

    def search_paths(self):
        """Search paths to use when searching for header fix up correction
        files.

        Returns
        -------
        paths : `list`
            Directory paths to search. Can be an empty list if no special
            directories are defined.

        Notes
        -----
        Uses the classes ``default_search_path`` property if defined.
        """
        if self.default_search_path is not None:
            return [self.default_search_path]
        return []

    def is_key_ok(self, keyword):
        """Return `True` if the value associated with the named keyword is
        present in this header and defined.

        Parameters
        ----------
        keyword : `str`
            Keyword to check against header.

        Returns
        -------
        is_ok : `bool`
            `True` if the header is present and not-`None`. `False` otherwise.
        """
        return self.is_keyword_defined(self._header, keyword)

    def are_keys_ok(self, keywords):
        """Are the supplied keys all present and defined?

        Parameters
        ----------
        keywords : iterable of `str`
            Keywords to test.

        Returns
        -------
        all_ok : `bool`
            `True` if all supplied keys are present and defined.
        """
        for k in keywords:
            if not self.is_key_ok(k):
                return False
        return True

    def quantity_from_card(self, keywords, unit, default=None, minimum=None, maximum=None, checker=None):
        """Calculate a Astropy Quantity from a header card and a unit.

        Parameters
        ----------
        keywords : `str` or `list` of `str`
            Keyword to use from header.  If a list each keyword will be tried
            in turn until one matches.
        unit : `astropy.units.UnitBase`
            Unit of the item in the header.
        default : `float`, optional
            Default value to use if the header value is invalid.  Assumed
            to be in the same units as the value expected in the header.  If
            None, no default value is used.
        minimum : `float`, optional
            Minimum possible valid value, optional.  If the calculated value
            is below this value, the default value will be used.
        maximum : `float`, optional
            Maximum possible valid value, optional.  If the calculated value
            is above this value, the default value will be used.
        checker : `function`, optional
            Callback function to be used by the translator method in case the
            keyword is not present.  Function will be executed as if it is
            a method of the translator class.  Running without raising an
            exception will allow the default to be used. Should usually raise
            `KeyError`.

        Returns
        -------
        q : `astropy.units.Quantity`
            Quantity representing the header value.

        Raises
        ------
        KeyError
            The supplied header key is not present.
        """
        keywords = keywords if isinstance(keywords, list) else [keywords]
        for k in keywords:
            if self.is_key_ok(k):
                value = self._header[k]
                keyword = k
                break
        else:
            if checker is not None:
                try:
                    checker(self)
                    value = default
                    if value is not None:
                        value = u.Quantity(value, unit=unit)
                    return value
                except Exception:
                    pass
            raise KeyError(f"Could not find {keywords} in header")
        if isinstance(value, str):
            # Sometimes the header has the wrong type in it but this must
            # be a number if we are creating a quantity.
            value = float(value)
        self._used_these_cards(keyword)
        if default is not None:
            value = self.validate_value(value, default, maximum=maximum, minimum=minimum)
        return u.Quantity(value, unit=unit)

    def _join_keyword_values(self, keywords, delim="+"):
        """Join values of all defined keywords with the specified delimiter.

        Parameters
        ----------
        keywords : iterable of `str`
            Keywords to look for in header.
        delim : `str`, optional
            Character to use to join the values together.

        Returns
        -------
        joined : `str`
            String formed from all the keywords found in the header with
            defined values joined by the delimiter. Empty string if no
            defined keywords found.
        """
        values = []
        for k in keywords:
            if self.is_key_ok(k):
                values.append(self._header[k])
                self._used_these_cards(k)

        if values:
            joined = delim.join(str(v) for v in values)
        else:
            joined = ""

        return joined

    @cache_translation
    def to_detector_unique_name(self):
        """Return a unique name for the detector.

        Base class implementation attempts to combine ``detector_name`` with
        ``detector_group``.  Group is only used if not `None`.

        Can be over-ridden by specialist translator class.

        Returns
        -------
        name : `str`
            ``detector_group``_``detector_name`` if ``detector_group`` is
            defined, else the ``detector_name`` is assumed to be unique.
            If neither return a valid value an exception is raised.

        Raises
        ------
        NotImplementedError
            Raised if neither detector_name nor detector_group is defined.
        """
        name = self.to_detector_name()
        group = self.to_detector_group()

        if group is None and name is None:
            raise NotImplementedError("Can not determine unique name from detector_group and detector_name")

        if group is not None:
            return f"{group}_{name}"

        return name

    @cache_translation
    def to_exposure_group(self):
        """Return the group label associated with this exposure.

        Base class implementation returns the ``exposure_id`` in string
        form.  A subclass may do something different.

        Returns
        -------
        name : `str`
            The ``exposure_id`` converted to a string.
        """
        exposure_id = self.to_exposure_id()
        if exposure_id is None:
            return None
        else:
            return str(exposure_id)

    @cache_translation
    def to_observation_reason(self):
        """Return the reason this observation was taken.

        Base class implementation returns the ``science`` if the
        ``observation_type`` is science, else ``unknown``.
        A subclass may do something different.

        Returns
        -------
        name : `str`
            The reason for this observation.
        """
        obstype = self.to_observation_type()
        if obstype == "science":
            return "science"
        return "unknown"

    @cache_translation
    def to_observing_day(self):
        """Return the YYYYMMDD integer corresponding to the observing day.

        Base class implementation uses the TAI date of the start of the
        observation.

        Returns
        -------
        day : `int`
            The observing day as an integer of form YYYYMMDD. If the header
            is broken and is unable to obtain a date of observation, ``0``
            is returned and the assumption is made that the problem will
            be caught elsewhere.
        """
        datetime_begin = self.to_datetime_begin()
        if datetime_begin is None:
            return 0
        return int(datetime_begin.tai.strftime("%Y%m%d"))

    @cache_translation
    def to_observation_counter(self):
        """Return an integer corresponding to how this observation relates
        to other observations.

        Base class implementation returns ``0`` to indicate that it is not
        known how an observatory will define a counter. Some observatories
        may not use the concept, others may use a counter that increases
        for every observation taken for that instrument, and others may
        define it to be a counter within an observing day.

        Returns
        -------
        sequence : `int`
            The observation counter. Always ``0`` for this implementation.
        """
        return 0

    @classmethod
    def determine_translatable_headers(cls, filename, primary=None):
        """Given a file return all the headers usable for metadata translation.

        This method can optionally be given a header from the file.  This
        header will generally be the primary header or a merge of the first
        two headers.

        In the base class implementation it is assumed that
        this supplied header is the only useful header for metadata translation
        and it will be returned unchanged if given. This can avoid
        unnecesarily re-opening the file and re-reading the header when the
        content is already known.

        If no header is supplied, a header will be read from the supplied
        file using `read_basic_metadata_from_file`, allowing it to merge
        the primary and secondary header of a multi-extension FITS file.
        Subclasses can read the header from the data file using whatever
        technique is best for that instrument.

        Subclasses can return multiple headers and ignore the externally
        supplied header. They can also merge it with another header and return
        a new derived header if that is required by the particular data file.
        There is no requirement for the supplied header to be used.

        Parameters
        ----------
        filename : `str`
            Path to a file in a format understood by this translator.
        primary : `dict`-like, optional
            The primary header obtained by the caller. This is sometimes
            already known, for example if a system is trying to bootstrap
            without already knowing what data is in the file. For many
            instruments where the primary header is the only relevant
            header, the primary header will be returned with no further
            action.

        Yields
        ------
        headers : iterator of `dict`-like
            A header usable for metadata translation. For this base
            implementation it will be either the supplied primary header
            or a header read from the file. This implementation will only
            ever yield a single header.

        Notes
        -----
        Each translator class can have code specifically tailored to its
        own file format. It is important not to call this method with
        an incorrect translator class. The normal paradigm is for the
        caller to have read the first header and then called
        `determine_translator()` on the result to work out which translator
        class to then call to obtain the real headers to be used for
        translation.
        """
        if primary is not None:
            yield primary
        else:
            # Prevent circular import by deferring
            from .file_helpers import read_basic_metadata_from_file

            # Merge primary and secondary header if they exist.
            yield read_basic_metadata_from_file(filename, -1)


def _make_abstract_translator_method(property, doc, return_typedoc, return_type):
    """Create a an abstract translation method for this property.

    Parameters
    ----------
    property : `str`
        Name of the translator for property to be created.
    doc : `str`
        Description of the property.
    return_typedoc : `str`
        Type string of this property (used in the doc string).
    return_type : `class`
        Type of this property.

    Returns
    -------
    m : `function`
        Translator method for this property.
    """
    def to_property(self):
        raise NotImplementedError(f"Translator for '{property}' undefined.")

    to_property.__doc__ = f"""Return value of {property} from headers.

    {doc}

    Returns
    -------
    {property} : `{return_typedoc}`
        The translated property.
    """
    return to_property


# Make abstract methods for all the translators methods.
# Unfortunately registering them as abstractmethods does not work
# as these assignments come after the class has been created.
# Assigning to __abstractmethods__ directly does work but interacts
# poorly with the metaclass automatically generating methods from
# _trivialMap and _constMap.

# Allow for concrete translator methods to exist in the base class
# These translator methods can be defined in terms of other properties
CONCRETE = set()

for name, description in PROPERTIES.items():
    method = f"to_{name}"
    if not MetadataTranslator.defined_in_this_class(method):
        setattr(MetadataTranslator, f"to_{name}",
                abstractmethod(_make_abstract_translator_method(name, *description[:3])))
    else:
        CONCRETE.add(method)


class StubTranslator(MetadataTranslator):
    """Translator where all the translations are stubbed out and issue
    warnings.

    This translator can be used as a base class whilst developing a new
    translator.  It allows testing to proceed without being required to fully
    define all translation methods.  Once complete the class should be
    removed from the inheritance tree.

    """
    pass


def _make_forwarded_stub_translator_method(cls, property, doc, return_typedoc, return_type):
    """Create a stub translation method for this property that calls the
    base method and catches `NotImplementedError`.

    Parameters
    ----------
    cls : `class`
        Class to use when referencing `super()`.  This would usually be
        `StubTranslator`.
    property : `str`
        Name of the translator for property to be created.
    doc : `str`
        Description of the property.
    return_typedoc : `str`
        Type string of this property (used in the doc string).
    return_type : `class`
        Type of this property.

    Returns
    -------
    m : `function`
        Stub translator method for this property.
    """
    method = f"to_{property}"

    def to_stub(self):
        parent = getattr(super(cls, self), method, None)
        try:
            if parent is not None:
                return parent()
        except NotImplementedError:
            pass

        warnings.warn(f"Please implement translator for property '{property}' for translator {self}",
                      stacklevel=3)
        return None

    to_stub.__doc__ = f"""Unimplemented forwarding translator for {property}.

    {doc}

    Calls the base class translation method and if that fails with
    `NotImplementedError` issues a warning reminding the implementer to
    override this method.

    Returns
    -------
    {property} : `None` or `{return_typedoc}`
        Always returns `None`.
    """
    return to_stub


# Create stub translation methods for each property.  These stubs warn
# rather than fail and should be overridden by translators.
for name, description in PROPERTIES.items():
    setattr(StubTranslator, f"to_{name}", _make_forwarded_stub_translator_method(StubTranslator,
                                                                                 name, *description[:3]))
