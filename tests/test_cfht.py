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

from astro_metadata_translator.tests import read_test_file, MetadataAssertHelper
from astro_metadata_translator import ObservationInfo

TESTDIR = os.path.abspath(os.path.dirname(__file__))


class MegaPrimeTestCase(unittest.TestCase, MetadataAssertHelper):
    datadir = os.path.join(TESTDIR, "data")

    def test_megaprime_translator(self):
        test_data = (("fitsheader-megaprime.yaml",
                      dict(telescope="CFHT 3.6m",
                           instrument="MegaPrime",
                           boresight_rotation_coord="sky",
                           boresight_rotation_angle=0*u.degree,
                           dark_time=615.0*u.s,
                           detector_exposure_id=37398350,
                           detector_name="ccd02",
                           detector_unique_name="ccd02",
                           detector_num=2,
                           detector_serial="8352-15-3",
                           exposure_id=1038843,
                           exposure_group="1038843",
                           exposure_time=615.037*u.s,
                           object="w2.+2+2",
                           observation_counter=1038843,
                           observation_id="1038843",
                           observation_type="science",
                           observation_reason="science",
                           observing_day=20081101,
                           physical_filter="i.MP9702",
                           pressure=617.65*u.hPa,
                           relative_humidity=39.77,
                           science_program="08BL05",
                           temperature=0.9*u.deg_C,
                           visit_id=1038843,
                           )),
                     ("fitsheader-megaprime-calexp-849375-14.yaml",
                      dict(telescope="CFHT 3.6m",
                           instrument="MegaPrime",
                           boresight_rotation_coord="sky",
                           boresight_rotation_angle=0*u.degree,
                           dark_time=300.0*u.s,
                           detector_exposure_id=30577599,
                           detector_name="ccd99",
                           detector_unique_name="ccd99",
                           detector_num=99,
                           detector_serial="8434-13-5",
                           exposure_id=849375,
                           exposure_group="849375",
                           exposure_time=300.202*u.s,
                           object="D3",
                           observation_counter=849375,
                           observation_id="849375",
                           observation_type="science",
                           observation_reason="science",
                           observing_day=20060520,
                           physical_filter="r",
                           pressure=615.79*u.hPa,
                           relative_humidity=15.76,
                           science_program="06AL01",
                           temperature=1.75*u.deg_C,
                           visit_id=849375,
                           )),
                     )
        for file, expected in test_data:
            with self.subTest(f"Testing {file}"):
                self.assertObservationInfoIsConsistent(file, dir=self.datadir, **expected)

    def test_megaprime_stripping(self):
        header = read_test_file("fitsheader-megaprime.yaml", dir=self.datadir)
        v1 = ObservationInfo(header)

        # Check that headers have been removed
        new_hdr = v1.stripped_header()
        self.assertNotIn("INSTRUME", new_hdr)
        self.assertNotIn("TELESCOP", new_hdr)
        self.assertIn("CCD", new_hdr)


if __name__ == "__main__":
    unittest.main()
