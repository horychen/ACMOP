import json
import os
import unittest

import bearingless_outer_rotor_spmsm_design


HERE = os.path.dirname(os.path.abspath(__file__))


class BearinglessOuterRotorSPMSMDesignTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open(os.path.join(HERE, "machine_specifications.json"), "r") as stream:
            specifications = json.load(stream)
        with open(os.path.join(HERE, "machine_simulation.json"), "r") as stream:
            simulations = json.load(stream)

        cls.spec = specifications[
            "OuterRotor SPMSM Q12p5ps4y1 Prototype"
        ]["Inputs"]
        cls.config = simulations["#0301 JMAG Non-Bearingless"]

    def build_template_and_variant(self):
        template = bearingless_outer_rotor_spmsm_design.bearingless_outer_rotor_spmsm_template(
            self.config,
            self.spec,
        )
        variant = bearingless_outer_rotor_spmsm_design.bearingless_outer_rotor_spmsm_design_variant(
            template=template,
            x_denorm=template.build_x_denorm(),
            counter=0,
            counter_loop=1,
        )
        return template, variant

    def test_template_uses_outer_rotor_radius_relations(self):
        template, _ = self.build_template_and_variant()
        GP = template.d["GP"]

        expected_stator_outer_radius = (
            GP["mm_r_ro"].value
            - GP["mm_d_ri"].value
            - GP["mm_d_pm"].value
            - GP["mm_d_sleeve"].value
            - GP["mm_d_mech_air_gap"].value
        )
        self.assertAlmostEqual(GP["mm_r_so"].value, expected_stator_outer_radius)
        self.assertAlmostEqual(
            GP["split_ratio"].value,
            GP["mm_r_so"].value / GP["mm_r_ro"].value,
        )
        self.assertGreater(GP["mm_r_so"].value, GP["mm_r_si"].value)

    def test_variant_assembles_outer_rotor_and_inner_stator(self):
        template, variant = self.build_template_and_variant()
        GP = template.d["GP"]

        self.assertIsNone(variant.shaft)
        self.assertAlmostEqual(variant.rotorCore.mm_r_ro, GP["mm_r_ro"].value)
        self.assertAlmostEqual(variant.stator_core.mm_r_so, GP["mm_r_so"].value)
        self.assertAlmostEqual(
            variant.sleeve.mm_r_si - variant.stator_core.mm_r_so,
            GP["mm_d_mech_air_gap"].value,
        )
        self.assertGreater(variant.rotorCore.mm_r_ro, variant.rotorCore.mm_r_ry_inner)
        self.assertGreater(
            variant.rotorCore.mm_r_ry_inner,
            variant.rotorCore.mm_r_magnet_inner,
        )

    def test_rotor_volume_is_an_annulus(self):
        template, _ = self.build_template_and_variant()
        GP = template.d["GP"]
        stack_length = 10.0
        expected = (
            3.141592653589793
            * (
                (GP["mm_r_ro"].value * 1e-3) ** 2
                - (GP["mm_r_ri"].value * 1e-3) ** 2
            )
            * stack_length
            * 1e-3
        )
        self.assertAlmostEqual(
            template.get_rotor_volume(stack_length=stack_length),
            expected,
        )

    def test_analytical_properties_use_si_units(self):
        template, _ = self.build_template_and_variant()
        EX = template.d["EX"]

        self.assertGreater(EX["stator_slot_area"], 0.0)
        self.assertLess(EX["stator_slot_area"], 1e-3)
        self.assertGreater(EX["end_winding_length_Lew"], 0.0)
        self.assertLess(EX["end_winding_length_Lew"], 1.0)

    def test_segment_span_is_clamped_to_one_segment_pitch(self):
        template, _ = self.build_template_and_variant()
        GP = template.d["GP"]
        segmented_spec = dict(self.spec)
        segmented_spec["no_segmented_magnets"] = 3
        GP["deg_alpha_rm"].value = 30.0
        GP["deg_alpha_rs"].value = 20.0
        GP["mm_d_rs"].value = 0.1

        bearingless_outer_rotor_spmsm_design.bearingless_outer_rotor_spmsm_design_variant.check_invalid_design(
            GP,
            segmented_spec,
        )
        self.assertEqual(GP["deg_alpha_rs"].value, 10.0)


if __name__ == "__main__":
    unittest.main()
