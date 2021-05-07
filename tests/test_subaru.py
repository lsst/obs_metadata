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

from astro_metadata_translator import merge_headers
from astro_metadata_translator.tests import MetadataAssertHelper, read_test_file

TESTDIR = os.path.abspath(os.path.dirname(__file__))


class HscTestCase(unittest.TestCase, MetadataAssertHelper):
    datadir = os.path.join(TESTDIR, "data")

    def test_hsc_translator(self):
        test_data = (("fitsheader-hsc.yaml",
                      dict(telescope="Subaru",
                           instrument="HSC",
                           boresight_rotation_coord="sky",
                           dark_time=30.0*u.s,
                           detector_exposure_id=180804850,
                           detector_name="12",
                           detector_unique_name="1_12",
                           detector_num=50,
                           detector_serial="120",
                           exposure_id=904024,
                           exposure_group="904024",
                           exposure_time=30.0*u.s,
                           object="STRIPE82L",
                           observation_counter=904024,
                           observation_id="HSCA90402400",
                           observation_type="science",
                           observation_reason="science",
                           observing_day=20131102,
                           physical_filter="HSC-I",
                           pressure=621.7*u.hPa,
                           relative_humidity=33.1,
                           science_program="o13015",
                           temperature=272.35*u.K,
                           visit_id=904024,
                           )),
                     ("fitsheader-hsc-HSCA04090107.yaml",
                      dict(telescope="Subaru",
                           instrument="HSC",
                           boresight_rotation_coord="sky",
                           dark_time=150.0*u.s,
                           detector_exposure_id=8180037,
                           detector_name="07",
                           detector_unique_name="1_07",
                           detector_num=37,
                           detector_serial="061",
                           exposure_id=40900,
                           exposure_group="40900",
                           exposure_time=150.0*u.s,
                           object="SSP-Wide",
                           observation_counter=40900,
                           observation_id="HSCA04090000",
                           observation_type="science",
                           observation_reason="science",
                           observing_day=20151010,
                           physical_filter="HSC-R",
                           pressure=625.4*u.hPa,
                           relative_humidity=8.6,
                           science_program="o15426",
                           temperature=278.35*u.K,
                           visit_id=40900,
                           )),
                     )
        for file, expected in test_data:
            with self.subTest(f"Testing {file}"):
                self.assertObservationInfoIsConsistent(file, dir=self.datadir, **expected)

    def test_suprimecam_translator(self):
        # In this case the airmass is average during observation
        # but it looks like ALTITUDE is from a different time so loosen amdelta
        test_data = (("fitsheader-suprimecam-CORR40535770.yaml",
                      dict(telescope="Subaru",
                           instrument="SuprimeCam",
                           boresight_rotation_coord="unknown",
                           dark_time=200.0*u.s,
                           detector_exposure_id=535770,
                           detector_name="nausicaa",
                           detector_unique_name="nausicaa",
                           detector_num=0,
                           detector_serial="w67c1",
                           exposure_id=53577,
                           exposure_group="53577",
                           exposure_time=200.0*u.s,
                           object="Ecliptic Deep Field",
                           observation_counter=53577,
                           observation_id="SUPE00535770",
                           observation_type="science",
                           observation_reason="science",
                           observing_day=20070423,
                           physical_filter="W-S-R+",
                           pressure=621.5*u.hPa,
                           relative_humidity=4.9,
                           science_program="o07222",
                           temperature=273.15*u.K,
                           visit_id=53577,
                           wcs_params=dict(amdelta=0.015))),
                     )
        for file, expected in test_data:
            with self.subTest(f"Testing {file}"):
                self.assertObservationInfoIsConsistent(file, dir=self.datadir, **expected)

    def test_merging_hsc(self):
        files = ("fitsheader-hsc-HSCA04090107.yaml", "fitsheader-hsc.yaml")
        headers = [read_test_file(f, dir=self.datadir) for f in files]
        merged = merge_headers(headers, mode="first", sort=False)

        # The MJD-STR should come from the first file
        self.assertAlmostEqual(merged["MJD-STR"], 57305.34729859)

        # If we sort then MJD-STR should come from the oldest file
        merged = merge_headers(headers, mode="first", sort=True)
        self.assertAlmostEqual(merged["MJD-STR"], 56598.26106374757)

        # Drop headers that differ, MJD-STR should not appear
        merged = merge_headers(headers, mode="drop", sort=True)
        self.assertNotIn("MJD-STR", merged)

        # Drop but retain first MJD-STR without sorting
        merged = merge_headers(headers, mode="drop", sort=False, first=["MJD-STR", "UT-STR"])
        self.assertAlmostEqual(merged["MJD-STR"], 57305.34729859)
        self.assertEqual(merged["UT-STR"], "08:20:06.598")

        # Drop but retain first MJD-STR
        merged = merge_headers(headers, mode="drop", sort=True, first=["MJD-STR", "UT-STR"])
        self.assertAlmostEqual(merged["MJD-STR"], 56598.26106374757)
        self.assertEqual(merged["UT-STR"], "06:15:55.908")


if __name__ == "__main__":
    unittest.main()
