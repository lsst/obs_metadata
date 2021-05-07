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

import os.path
import unittest
import astropy.units as u

from astro_metadata_translator.tests import MetadataAssertHelper

TESTDIR = os.path.abspath(os.path.dirname(__file__))


class SdssTestCase(unittest.TestCase, MetadataAssertHelper):
    datadir = os.path.join(TESTDIR, "data")

    def test_sdss_translator(self):
        test_data = (("fitsheader-sdss-fpC-006377-g4-0399.yaml",
                      dict(telescope="SDSS 2.5m",
                           instrument="Imager on SDSS 2.5m",
                           boresight_rotation_coord="sky",
                           dark_time=0.0*u.s,
                           detector_exposure_id=6377140399,
                           detector_name="54",
                           detector_unique_name="g4",
                           detector_group="4",
                           detector_num=15,
                           detector_serial="UNKNOWN",
                           exposure_id=6377,
                           exposure_group="6377",
                           exposure_time=53.907456*u.s,
                           object="82 S",
                           observation_counter=0,
                           observation_id="6377 4 g 407",
                           observation_type="science",
                           observation_reason="science",
                           observing_day=20060920,
                           physical_filter="g",
                           pressure=None,
                           relative_humidity=None,
                           science_program="82 S",
                           temperature=None,
                           visit_id=6377,
                           wcs_params=dict(max_sep=10.0))),
                     )
        for file, expected in test_data:
            with self.subTest(f"Testing {file}"):
                self.assertObservationInfoIsConsistent(file, dir=self.datadir, **expected)


if __name__ == "__main__":
    unittest.main()
