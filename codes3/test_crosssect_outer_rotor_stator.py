import unittest

import CrossSectOuterRotorStator
import CrossSectStator


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


class CrossSectOuterRotorStatorTests(unittest.TestCase):
    def setUp(self):
        self.core = CrossSectOuterRotorStator.CrossSectOuterRotorStator(
            deg_alpha_st=12.0,
            deg_alpha_sto=6.0,
            mm_r_so=60.0,
            mm_d_sto=3.0,
            mm_d_stt=5.0,
            mm_d_st=20.0,
            mm_d_sy=10.0,
            mm_w_st=8.0,
            Q=12,
        )

    def test_sector_core_has_one_closed_region(self):
        token = self.core.draw(RecordingDrawer())
        self.assertEqual(len(token["list_regions"]), 1)
        self.assertEqual(len(token["list_regions"][0]), 8)
        self.assertGreater(token["innerCoord"][0], self.core.mm_r_sy_inner)
        self.assertLess(token["innerCoord"][0], self.core.mm_r_so)

    def test_whole_core_contains_all_teeth(self):
        token = self.core.draw(RecordingDrawer(), bool_draw_whole_model=True)
        self.assertEqual(len(token["list_regions"]), 1)
        self.assertEqual(len(token["list_regions"][0]), self.core.Q * 8 + 2)

    def test_winding_area_and_region_count(self):
        winding = CrossSectOuterRotorStator.CrossSectOuterRotorStatorWinding(
            stator_core=self.core
        )
        area = winding.draw(RecordingDrawer(), bool_re_evaluate=True)
        self.assertGreater(area, 0.0)

        sector_token = winding.draw(RecordingDrawer())
        self.assertEqual(len(sector_token["list_regions"]), 2)

        whole_token = winding.draw(RecordingDrawer(), bool_draw_whole_model=True)
        self.assertEqual(len(whole_token["list_regions"]), 2 * self.core.Q)

    def test_invalid_radial_stack_is_rejected(self):
        with self.assertRaises(ValueError):
            CrossSectOuterRotorStator.CrossSectOuterRotorStator(
                mm_r_so=30.0,
                mm_d_stt=5.0,
                mm_d_st=15.0,
                mm_d_sy=10.0,
            )

    def test_tooth_body_must_leave_the_slot_bottom_open(self):
        core = CrossSectOuterRotorStator.CrossSectOuterRotorStator(
            deg_alpha_st=28.0,
            deg_alpha_sto=14.0,
            mm_r_so=60.0,
            mm_d_sto=1.4137,
            mm_d_stt=2.1206,
            mm_d_st=23.7423,
            mm_d_sy=14.1372,
            mm_w_st=18.8496,
            Q=12,
        )
        with self.assertRaisesRegex(ValueError, "intersect at the slot bottom"):
            core._build_points()

    def test_valid_stator_has_a_finite_slot_bottom_span(self):
        points = self.core._build_points()
        self.assertGreater(points["slot_bottom_span"], 0.0)

    def test_existing_stator_module_reexports_outer_rotor_classes(self):
        self.assertIs(
            CrossSectStator.CrossSectOuterRotorStator,
            CrossSectOuterRotorStator.CrossSectOuterRotorStator,
        )
        self.assertIs(
            CrossSectStator.CrossSectInnerStator,
            CrossSectOuterRotorStator.CrossSectOuterRotorStator,
        )


if __name__ == "__main__":
    unittest.main()
