import unittest

import numpy as np

import CrossSectOuterNotchedRotor


class RecordingDrawer:
    def __init__(self):
        self.sketches = []

    def getSketch(self, name, color):
        self.sketches.append((name, color))

    @staticmethod
    def drawLine(start, end):
        return [("line", tuple(start), tuple(end))]

    @staticmethod
    def drawArc(center, start, end):
        return [("arc", tuple(center), tuple(start), tuple(end))]


class CrossSectOuterNotchedRotorTests(unittest.TestCase):
    def setUp(self):
        self.rotor = CrossSectOuterNotchedRotor.CrossSectOuterNotchedRotor(
            mm_r_ro=74.0,
            mm_d_ri=8.0,
            mm_d_pm=4.0,
            mm_d_rp=4.0,
            mm_d_rs=0.0,
            deg_alpha_rm=30.0,
            deg_alpha_rs=30.0,
            p=5,
            s=1,
        )
        self.magnet = CrossSectOuterNotchedRotor.CrossSectOuterNotchedMagnet(
            notched_rotor=self.rotor
        )

    def test_radial_order_matches_outer_rotor_topology(self):
        self.assertGreater(self.rotor.mm_r_ro, self.rotor.mm_r_ry_inner)
        self.assertGreater(self.rotor.mm_r_ry_inner, self.rotor.mm_r_magnet_inner)

    def test_full_pole_arc_does_not_retrace_the_first_radial_edge(self):
        rotor = CrossSectOuterNotchedRotor.CrossSectOuterNotchedRotor(
            mm_d_pm=4.0,
            deg_alpha_rm=36.0,
            deg_alpha_rs=36.0,
            mm_d_ri=8.0,
            mm_r_ro=74.0,
            mm_d_rp=4.0,
            mm_d_rs=0.0,
            p=5,
            s=1,
        )

        token = rotor.draw(RecordingDrawer())

        radial_lines = [
            segment
            for segment in token["list_regions"][0]
            if segment[0] == "line"
        ]
        self.assertEqual(2, len(radial_lines))
        self.assertAlmostEqual(rotor.mm_r_ro, np.hypot(*radial_lines[0][1]))
        self.assertAlmostEqual(
            rotor.mm_r_ry_inner,
            np.hypot(*radial_lines[0][2]),
        )

    def test_core_draws_sector_and_whole_model(self):
        sector = self.rotor.draw(RecordingDrawer())
        whole = self.rotor.draw(RecordingDrawer(), bool_draw_whole_model=True)
        self.assertEqual(len(sector["list_regions"]), 1)
        self.assertGreater(len(sector["list_regions"][0]), 0)
        self.assertGreater(
            len(whole["list_regions"][0]),
            len(sector["list_regions"][0]),
        )
        radial_lines_to_outer_radius = []
        for segment in whole["list_regions"][0]:
            if segment[0] != "line":
                continue
            radii = [np.hypot(*segment[1]), np.hypot(*segment[2])]
            if max(radii) >= self.rotor.mm_r_ro - 1e-9:
                radial_lines_to_outer_radius.append(segment)
        self.assertEqual(radial_lines_to_outer_radius, [])

    def test_magnet_area_and_region_count(self):
        expected_area = (
            self.rotor.deg_alpha_rm / (180.0 / self.rotor.p)
            * np.pi
            * (
                self.rotor.mm_r_ry_inner**2
                - self.rotor.mm_r_magnet_inner**2
            )
        )
        area = self.magnet.draw(RecordingDrawer(), bool_re_evaluate=True)
        self.assertAlmostEqual(area, expected_area)

        sector = self.magnet.draw(RecordingDrawer())
        whole = self.magnet.draw(RecordingDrawer(), bool_draw_whole_model=True)
        self.assertEqual(len(sector["list_regions"]), self.rotor.s)
        self.assertEqual(len(whole["list_regions"]), 2 * self.rotor.p * self.rotor.s)

    def test_segmented_magnets_are_supported(self):
        rotor = CrossSectOuterNotchedRotor.CrossSectOuterNotchedRotor(
            mm_r_ro=74.0,
            mm_d_ri=8.0,
            mm_d_pm=4.0,
            mm_d_rp=3.0,
            mm_d_rs=1.0,
            deg_alpha_rm=30.0,
            deg_alpha_rs=9.0,
            p=5,
            s=3,
        )
        magnet = CrossSectOuterNotchedRotor.CrossSectOuterNotchedMagnet(
            notched_rotor=rotor
        )
        self.assertEqual(len(magnet.draw(RecordingDrawer())["list_regions"]), 3)
        self.assertEqual(
            len(
                magnet.draw(
                    RecordingDrawer(), bool_draw_whole_model=True
                )["list_regions"]
            ),
            30,
        )

    def test_inner_sleeve_is_inside_magnets(self):
        sleeve = CrossSectOuterNotchedRotor.CrossSectOuterRotorSleeve(
            notched_magnet=self.magnet,
            d_sleeve=1.0,
        )
        self.assertLess(sleeve.mm_r_si, sleeve.mm_r_so)
        self.assertEqual(sleeve.mm_r_so, self.rotor.mm_r_magnet_inner)
        self.assertGreater(len(sleeve.draw(RecordingDrawer())["list_regions"][0]), 0)
        whole = sleeve.draw(RecordingDrawer(), bool_draw_whole_model=True)
        self.assertEqual(len(whole["list_regions"][0]), 4)
        self.assertTrue(all(segment[0] == "arc" for segment in whole["list_regions"][0]))

    def test_invalid_notch_depth_is_rejected(self):
        with self.assertRaises(ValueError):
            CrossSectOuterNotchedRotor.CrossSectOuterNotchedRotor(
                mm_d_pm=4.0,
                mm_d_rp=5.0,
            )


if __name__ == "__main__":
    unittest.main()
