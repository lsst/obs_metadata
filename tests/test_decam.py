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


class DecamTestCase(unittest.TestCase, MetadataAssertHelper):
    datadir = os.path.join(TESTDIR, "data")

    def test_decam_translator(self):
        test_data = (("fitsheader-decam.yaml",
                      dict(telescope="CTIO 4.0-m telescope",
                           instrument="DECam",
                           boresight_rotation_coord="sky",
                           dark_time=201.15662*u.s,
                           detector_exposure_id=22938825,
                           detector_name="1",
                           detector_unique_name="S1",
                           detector_group="S",
                           detector_num=25,
                           detector_serial="S3-111_107419-8-3",
                           exposure_id=229388,
                           exposure_group="229388",
                           exposure_time=200.0*u.s,
                           object="DES supernova hex SN-S1 tiling 22",
                           observation_counter=229388,
                           observation_id="ct4m20130901t060255",
                           observation_type="science",
                           observation_reason="science",
                           observing_day=20130901,
                           physical_filter="z DECam SDSS c0004 9260.0 1520.0",
                           pressure=779.0*u.hPa,
                           relative_humidity=23.0,
                           science_program="2012B-0001",
                           temperature=11.9*u.deg_C,
                           visit_id=229388,
                           wcs_params=dict(max_sep=1.5))),
                     ("fitsheader-decam-0160496.yaml",
                      dict(telescope="CTIO 4.0-m telescope",
                           instrument="DECam",
                           boresight_rotation_coord="sky",
                           boresight_rotation_angle=90*u.degree,
                           dark_time=0.0407898*u.s,
                           detector_exposure_id=16049625,
                           detector_name="1",
                           detector_unique_name="S1",
                           detector_group="S",
                           detector_num=25,
                           detector_serial="S3-111_107419-8-3",
                           exposure_id=160496,
                           exposure_group="160496",
                           exposure_time=0.0*u.s,
                           object="postflats-BIAS",
                           observation_counter=160496,
                           observation_id="ct4m20121211t220632",
                           observation_type="zero",
                           observation_reason="unknown",
                           observing_day=20121211,
                           physical_filter="solid plate 0.0 0.0",  # corrected value
                           pressure=777.0*u.hPa,
                           relative_humidity=38.0,
                           science_program="2012B-0416",
                           temperature=17.0*u.deg_C,
                           visit_id=160496,
                           wcs_params=dict(max_sep=1.5))),
                     ("fitsheader-decam-calexp-0412037_10.yaml",
                      dict(telescope="CTIO 4.0-m telescope",
                           instrument="DECam",
                           boresight_rotation_coord="sky",
                           boresight_rotation_angle=90*u.degree,
                           dark_time=87.1054702*u.s,
                           detector_exposure_id=41203701,
                           detector_name="29",
                           detector_unique_name="S29",
                           detector_group="S",
                           detector_num=1,
                           detector_serial="S3-06_123195-15-3",
                           exposure_id=412037,
                           exposure_group="412037",
                           exposure_time=86.0*u.s,
                           object="Blind15A_03",
                           observation_counter=412037,
                           observation_id="ct4m20150220t004721",
                           observation_type="science",
                           observation_reason="science",
                           observing_day=20150220,
                           physical_filter="g",
                           pressure=777.0*u.hPa,
                           relative_humidity=76.0,
                           science_program="2015A-0608",
                           temperature=9.0*u.deg_C,
                           visit_id=412037,
                           wcs_params=dict(max_sep=5.0))),
                     ("fitsheader-decam-instcal-c4d_190402_050618_ooi_VR_v1.yaml",
                      dict(telescope="CTIO 4.0-m telescope",
                           instrument="DECam",
                           boresight_rotation_coord="sky",
                           boresight_rotation_angle=90*u.degree,
                           dark_time=120.7646399*u.s,
                           detector_exposure_id=84529101,
                           detector_name="29",
                           detector_unique_name="S29",
                           detector_group="S",
                           detector_num=1,
                           detector_serial="S3-06_123195-15-3",
                           exposure_id=845291,
                           exposure_group="845291",
                           exposure_time=120.0*u.s,
                           object="",
                           observation_counter=845291,
                           observation_id="ct4m20190402t050618",
                           observation_type="science",
                           observation_reason="science",
                           observing_day=20190402,
                           physical_filter="VR DECam c0007 6300.0 2600.0",
                           pressure=779.0*u.hPa,
                           relative_humidity=38.0,
                           science_program="2019A-0337",
                           temperature=15.1*u.deg_C,
                           visit_id=845291,
                           wcs_params=dict(max_sep=5.0))),
                     )
        for file, expected in test_data:
            with self.subTest(f"Testing {file} from dict-like"):
                self.assertObservationInfoIsConsistent(file, dir=self.datadir, **expected)


if __name__ == "__main__":
    unittest.main()
