# This file is part of obs_metadata.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Generically useful translation helpers which translation classes
can use.

They are written as free functions as methods that will be attached
to classes dynamically. They are assumed to be called from
`MetadataTranslator` classes and have access to the standard properties.

"""

__all__ = ("to_location_via_telescope_name", )

from astropy.coordinates import EarthLocation


def to_location_via_telescope_name(self):
    """Calculate the observatory location via the telescope name.

    Returns
    -------
    loc : `astropy.coordinates.EarthLocation`
        Location of the observatory.
    """
    return EarthLocation.of_site(self.to_telescope())