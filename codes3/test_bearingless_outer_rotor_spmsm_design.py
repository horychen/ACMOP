import json
import math
import os
import unittest
from copy import deepcopy
from itertools import product

import bearingless_outer_rotor_spmsm_design
import pyrhonen_procedure_as_function


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

    def build_template_and_variant(self, spec=None, x_denorm=None):
        spec = deepcopy(self.spec if spec is None else spec)
        template = bearingless_outer_rotor_spmsm_design.bearingless_outer_rotor_spmsm_template(
            deepcopy(self.config),
            spec,
        )
        if x_denorm is None:
            x_denorm = template.build_x_denorm()
        variant = bearingless_outer_rotor_spmsm_design.bearingless_outer_rotor_spmsm_design_variant(
            template=template,
            x_denorm=x_denorm,
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

    def test_fixed_envelope_has_one_radial_closure_parameter(self):
        template, _ = self.build_template_and_variant()
        GP = template.d["GP"]

        self.assertEqual(GP["mm_r_ro"].type, "fixed")
        self.assertEqual(GP["mm_r_si"].type, "fixed")
        self.assertEqual(GP["mm_d_sy"].type, "derived")
        self.assertEqual(GP["mm_r_so"].type, "derived")
        self.assertEqual(GP["mm_r_ri"].type, "derived")

    def test_initial_yoke_uses_magnet_pole_arc_ratio(self):
        template, _ = self.build_template_and_variant()
        GP = template.d["GP"]
        SI = template.SI
        expected_yoke_depth = (
            SI["guess_air_gap_flux_density_Bg"]
            * 3.141592653589793
            * (2.0 * GP["mm_r_so"].value * 1e-3)
            * SI["rotor_magnet_pole_arc_ratio"]
            / (4.0 * SI["guess_stator_yoke_flux_density_Bsy"] * SI["p"])
            * 1e3
        )

        self.assertAlmostEqual(expected_yoke_depth, GP["mm_d_sy"].value)

    def test_initial_tooth_width_uses_slot_flux_distribution(self):
        template, _ = self.build_template_and_variant()
        GP = template.d["GP"]
        SI = template.SI
        slot_flux_angle = SI["p"] * 3.141592653589793 / SI["Qs"]
        slot_flux_factor = abs(math.sin(slot_flux_angle) / slot_flux_angle)
        expected_tooth_width = (
            SI["guess_air_gap_flux_density_Bg"]
            * 3.141592653589793
            * 2.0
            * GP["mm_r_so"].value
            * slot_flux_factor
            / (
                SI["guess_stator_tooth_flux_density_Bst"]
                * SI["Qs"]
                * SI["lamination_stacking_factor_kFe"]
            )
        )

        self.assertAlmostEqual(expected_tooth_width, GP["mm_w_st"].value)
        self.assertAlmostEqual(
            slot_flux_factor,
            template.d["EX"]["stator_slot_flux_factor"],
        )

    def test_configured_circumferential_slot_opening_is_two_mm(self):
        template, _ = self.build_template_and_variant()

        self.assertAlmostEqual(
            2.0,
            template.d["EX"]["mm_stator_slot_opening_width"],
        )

    def test_small_stator_bore_has_nonempty_geometry_bounds(self):
        spec = deepcopy(self.spec)
        spec["mm_stator_inner_radius"] = 15.0
        template, variant = self.build_template_and_variant(spec=spec)
        GP = template.d["GP"]

        self.assertLessEqual(
            template.original_template_neighbor_bounds["mm_w_st"][0],
            GP["mm_w_st"].value,
        )
        self.assertGreaterEqual(
            template.original_template_neighbor_bounds["mm_w_st"][1],
            GP["mm_w_st"].value,
        )
        self.assertLessEqual(
            template.original_template_neighbor_bounds["mm_d_st"][0],
            GP["mm_d_st"].value,
        )
        self.assertGreaterEqual(
            template.original_template_neighbor_bounds["mm_d_st"][1],
            GP["mm_d_st"].value,
        )
        variant.stator_core._build_points()

    def test_configured_stack_length_is_used_for_zq(self):
        template, _ = self.build_template_and_variant()
        EX = template.d["EX"]

        self.assertEqual(self.spec["mm_stack_length"], EX["mm_template_stack_length"])
        self.assertIsInstance(EX["DriveW_zQ"], int)
        self.assertGreater(EX["DriveW_zQ"], 0)

    def test_zq_is_not_limited_to_990(self):
        template, _ = self.build_template_and_variant()
        GP = template.d["GP"]
        EX = template.d["EX"]
        magnet_inner_diameter = 2e-3 * (
            GP["mm_r_ro"].value
            - GP["mm_d_ri"].value
            - GP["mm_d_pm"].value
        )

        zQ = pyrhonen_procedure_as_function.get_zQ(
            template.SI,
            EX["wily"],
            magnet_inner_diameter,
            GP["mm_r_so"].value * 2e-3,
            specified_mm_stack_length=0.1,
        )

        self.assertIsInstance(zQ, int)
        self.assertGreater(zQ, 990)

    def test_only_stator_slot_parameters_are_free(self):
        spec = deepcopy(self.spec)
        spec["optimization_parameterization"] = "stator_slot_only"
        template, _ = self.build_template_and_variant(spec=spec)
        template.build_x_denorm()
        GP = template.d["GP"]

        self.assertEqual(
            list(template.x_denorm_dict),
            ["deg_alpha_st", "mm_w_st", "mm_d_sto", "mm_d_st"],
        )
        for key in ("mm_d_pm", "mm_d_ri", "mm_d_mech_air_gap", "mm_d_sleeve"):
            self.assertEqual(GP[key].type, "fixed")
        for key in ("deg_alpha_rm", "mm_d_rp", "deg_alpha_rs"):
            self.assertEqual(GP[key].type, "derived")

    def test_notched_rotor_parameters_are_free(self):
        template, variant = self.build_template_and_variant()
        template.build_x_denorm()
        GP = template.d["GP"]

        self.assertEqual(
            list(template.x_denorm_dict),
            [
                "deg_alpha_st",
                "mm_w_st",
                "mm_d_sto",
                "mm_d_st",
                "deg_alpha_rm",
                "mm_d_rp",
            ],
        )
        self.assertLess(GP["deg_alpha_rm"].value, 180.0 / template.SI["p"])
        self.assertLess(GP["mm_d_rp"].value, GP["mm_d_pm"].value)
        self.assertEqual(GP["deg_alpha_rm"].type, "free")
        self.assertEqual(GP["mm_d_rp"].type, "free")
        self.assertAlmostEqual(
            variant.rotorCore.deg_alpha_rm,
            GP["deg_alpha_rm"].value,
        )
        self.assertAlmostEqual(variant.rotorCore.mm_d_rp, GP["mm_d_rp"].value)

    def test_all_search_space_corners_are_geometrically_valid(self):
        reference = bearingless_outer_rotor_spmsm_design.bearingless_outer_rotor_spmsm_template(
            deepcopy(self.config),
            deepcopy(self.spec),
        )
        for choices in product((0, 1), repeat=len(reference.bounds_denorm)):
            with self.subTest(choices=choices):
                template = bearingless_outer_rotor_spmsm_design.bearingless_outer_rotor_spmsm_template(
                    deepcopy(self.config),
                    deepcopy(self.spec),
                )
                x_denorm = [
                    bounds[choice]
                    for bounds, choice in zip(template.bounds_denorm, choices)
                ]
                bearingless_outer_rotor_spmsm_design.bearingless_outer_rotor_spmsm_design_variant(
                    template=template,
                    x_denorm=x_denorm,
                    counter=0,
                    counter_loop=1,
                )

    def test_fixed_envelope_keeps_bore_and_outer_radius_fixed(self):
        template = bearingless_outer_rotor_spmsm_design.bearingless_outer_rotor_spmsm_template(
            deepcopy(self.config),
            deepcopy(self.spec),
        )
        initial_r_ro = template.d["GP"]["mm_r_ro"].value
        initial_r_si = template.d["GP"]["mm_r_si"].value
        initial_d_sy = template.d["GP"]["mm_d_sy"].value
        x_denorm = template.build_x_denorm()
        tooth_depth_index = list(template.x_denorm_dict).index("mm_d_st")
        x_denorm[tooth_depth_index] += 0.5

        bearingless_outer_rotor_spmsm_design.bearingless_outer_rotor_spmsm_design_variant(
            template=template,
            x_denorm=x_denorm,
            counter=0,
            counter_loop=1,
        )
        GP = template.d["GP"]
        self.assertAlmostEqual(GP["mm_r_ro"].value, initial_r_ro)
        self.assertAlmostEqual(GP["mm_r_si"].value, initial_r_si)
        self.assertAlmostEqual(GP["mm_d_sy"].value, initial_d_sy - 0.5)

    def test_sleeve_can_be_disabled(self):
        spec = deepcopy(self.spec)
        spec["use_sleeve"] = False
        template, variant = self.build_template_and_variant(spec=spec)
        GP = template.d["GP"]

        self.assertEqual(GP["mm_d_sleeve"].value, 0.0)
        self.assertEqual(GP["mm_d_sleeve"].type, "fixed")
        self.assertIsNone(variant.sleeve)
        self.assertAlmostEqual(
            variant.rotorCore.mm_r_magnet_inner - variant.stator_core.mm_r_so,
            GP["mm_d_mech_air_gap"].value,
        )

    def test_fixed_stator_bore_derives_rotor_outer_radius(self):
        spec = deepcopy(self.spec)
        spec["radial_constraint_mode"] = "fixed_stator_bore"
        template = bearingless_outer_rotor_spmsm_design.bearingless_outer_rotor_spmsm_template(
            deepcopy(self.config),
            spec,
        )
        GP = template.d["GP"]
        initial_r_ro = GP["mm_r_ro"].value
        initial_r_so = GP["mm_r_so"].value
        x_denorm = template.build_x_denorm()
        tooth_depth_index = list(template.x_denorm_dict).index("mm_d_st")
        x_denorm[tooth_depth_index] += 0.5

        bearingless_outer_rotor_spmsm_design.bearingless_outer_rotor_spmsm_design_variant(
            template=template,
            x_denorm=x_denorm,
            counter=0,
            counter_loop=1,
        )
        self.assertEqual(GP["mm_r_ro"].type, "derived")
        self.assertEqual(GP["mm_r_si"].type, "fixed")
        self.assertAlmostEqual(GP["mm_r_so"].value, initial_r_so + 0.5)
        self.assertAlmostEqual(GP["mm_r_ro"].value, initial_r_ro + 0.5)

    def test_fixed_outer_radius_derives_stator_bore(self):
        spec = deepcopy(self.spec)
        spec["radial_constraint_mode"] = "fixed_outer_radius"
        template = bearingless_outer_rotor_spmsm_design.bearingless_outer_rotor_spmsm_template(
            deepcopy(self.config),
            spec,
        )
        GP = template.d["GP"]
        initial_r_ro = GP["mm_r_ro"].value
        initial_r_si = GP["mm_r_si"].value
        initial_d_sy = GP["mm_d_sy"].value
        x_denorm = template.build_x_denorm()
        tooth_depth_index = list(template.x_denorm_dict).index("mm_d_st")
        x_denorm[tooth_depth_index] += 0.5

        bearingless_outer_rotor_spmsm_design.bearingless_outer_rotor_spmsm_design_variant(
            template=template,
            x_denorm=x_denorm,
            counter=0,
            counter_loop=1,
        )
        self.assertEqual(GP["mm_r_ro"].type, "fixed")
        self.assertEqual(GP["mm_r_si"].type, "derived")
        self.assertEqual(GP["mm_d_sy"].type, "fixed")
        self.assertAlmostEqual(GP["mm_r_ro"].value, initial_r_ro)
        self.assertAlmostEqual(GP["mm_d_sy"].value, initial_d_sy)
        self.assertAlmostEqual(GP["mm_r_si"].value, initial_r_si - 0.5)

    def test_variant_assembles_outer_rotor_and_inner_stator(self):
        template, variant = self.build_template_and_variant()
        GP = template.d["GP"]

        self.assertIsNone(variant.shaft)
        self.assertAlmostEqual(variant.rotorCore.mm_r_ro, GP["mm_r_ro"].value)
        self.assertAlmostEqual(variant.stator_core.mm_r_so, GP["mm_r_so"].value)
        rotating_inner_radius = (
            variant.sleeve.mm_r_si
            if variant.sleeve is not None
            else variant.rotorCore.mm_r_magnet_inner
        )
        self.assertAlmostEqual(
            rotating_inner_radius - variant.stator_core.mm_r_so,
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
