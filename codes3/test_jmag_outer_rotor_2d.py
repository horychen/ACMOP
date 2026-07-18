import json
import os
import unittest
from copy import deepcopy
from types import SimpleNamespace

import JMAG
import bearingless_outer_rotor_spmsm_design
import utility


HERE = os.path.dirname(os.path.abspath(__file__))


class RecordingJMAGDrawer:
    def __init__(self):
        self.sections = []
        self.sketches = []
        self.current_calculated = False
        self.saved = False
        self.bMirror = False
        self.iRotateCopy = 0

    def getSketch(self, name, color=None):
        self.sketches.append((name, color))

    @staticmethod
    def drawLine(start, end):
        return [("line", tuple(start), tuple(end))]

    @staticmethod
    def drawArc(center, start, end):
        return [("arc", tuple(center), tuple(start), tuple(end))]

    def prepareSection(self, token, **kwargs):
        self.sections.append(
            {
                "regions": len(token["list_regions"]),
                "mirror": self.bMirror,
                "copies": self.iRotateCopy,
                "kwargs": kwargs,
            }
        )
        return token["list_regions"]

    def calculate_excitation_current(self, acm_variant):
        self.current_calculated = True

    def save(self, name, description):
        self.saved = True

    @staticmethod
    def show(acm_variant, toString=False):
        return ""


class FakeSelection:
    def __init__(self):
        self.part_ids = []
        self.positions = []

    def SelectPart(self, part_id):
        self.part_ids.append(part_id)

    def SelectPartByPosition(self, x, y, z):
        self.positions.append((x, y, z))


class FakePartSet:
    def __init__(self):
        self.selection = FakeSelection()
        self.part_ids = []
        self.positions = []

    def SetMatcherType(self, matcher_type):
        self.matcher_type = matcher_type

    def ClearParts(self):
        self.part_ids = []
        self.positions = []

    def GetSelection(self):
        return self.selection

    def AddSelected(self, selection):
        self.part_ids = list(selection.part_ids)
        self.positions = list(selection.positions)


class FakeSetList:
    def __init__(self):
        self.sets = {}

    def CreatePartSet(self, name):
        self.sets[name] = FakePartSet()

    def GetSet(self, name):
        return self.sets[name]


class FakeGroupList:
    def __init__(self):
        self.groups = {}

    def CreateGroup(self, name):
        self.groups[name] = []

    def AddPartToGroup(self, name, part_id):
        self.groups[name].append(part_id)


class FakeModel:
    def __init__(self, part_ids):
        self.part_ids = part_ids
        self.set_list = FakeSetList()
        self.group_list = FakeGroupList()

    def GetPartIDs(self):
        return self.part_ids

    def GetSetList(self):
        return self.set_list

    def GetGroupList(self):
        return self.group_list


class FakeMaterial:
    def __init__(self):
        self.values = {}

    def SetValue(self, key, value):
        self.values[key] = value


class FakeStudy:
    def __init__(self):
        self.assignments = {}
        self.materials = {}

    def SetMaterialByName(self, part_name, material_name):
        self.assignments[part_name] = material_name
        self.materials.setdefault(part_name, FakeMaterial())

    def GetMaterial(self, part_name):
        return self.materials[part_name]


class JMAGOuterRotor2DTests(unittest.TestCase):
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

    def build_variant(self, use_sleeve=True):
        spec = deepcopy(self.spec)
        spec["use_sleeve"] = use_sleeve
        template = bearingless_outer_rotor_spmsm_design.bearingless_outer_rotor_spmsm_template(
            deepcopy(self.config),
            spec,
        )
        return bearingless_outer_rotor_spmsm_design.bearingless_outer_rotor_spmsm_design_variant(
            template,
            template.build_x_denorm(),
            counter=0,
            counter_loop=1,
        )

    def test_part_id_parser_supports_optional_sleeve(self):
        with_sleeve = list(range(1, 14))
        layout = JMAG.JMAG.parse_outer_rotor_part_ids(
            with_sleeve,
            p=2,
            s=1,
            Q=3,
            has_sleeve=True,
        )
        self.assertEqual(layout["rotor_core"], 1)
        self.assertEqual(layout["magnets"], [2, 3, 4, 5])
        self.assertEqual(layout["sleeve"], 6)
        self.assertEqual(layout["stator_core"], 7)
        self.assertEqual(layout["coils"], [8, 9, 10, 11, 12, 13])

        without_sleeve = list(range(1, 13))
        layout = JMAG.JMAG.parse_outer_rotor_part_ids(
            without_sleeve,
            p=2,
            s=1,
            Q=3,
            has_sleeve=False,
        )
        self.assertIsNone(layout["sleeve"])
        self.assertEqual(layout["stator_core"], 6)

    def test_part_id_parser_rejects_unexpected_count(self):
        with self.assertRaises(utility.ExceptionBadNumberOfParts):
            JMAG.JMAG.parse_outer_rotor_part_ids(
                [1, 2, 3],
                p=2,
                s=1,
                Q=3,
                has_sleeve=True,
            )

    def test_outer_rotor_draw_skips_shaft_and_supports_sleeve_switch(self):
        for use_sleeve, expected_sections in ((True, 5), (False, 4)):
            with self.subTest(use_sleeve=use_sleeve):
                variant = self.build_variant(use_sleeve=use_sleeve)
                drawer = RecordingJMAGDrawer()
                result = JMAG.JMAG.draw_outer_rotor_spmsm(drawer, variant)

                self.assertTrue(result)
                self.assertEqual(len(drawer.sections), expected_sections)
                self.assertTrue(drawer.current_calculated)
                self.assertTrue(drawer.saved)
                self.assertNotIn("Shaft", [name for name, _ in drawer.sketches])

    def test_preprocess_motion_region_contains_only_rotating_outer_parts(self):
        variant = self.build_variant(use_sleeve=True)
        variant.coils.draw(None, bool_re_evaluate=True)
        p = variant.template.SI["p"]
        s = variant.template.SI["no_segmented_magnets"]
        Q = variant.template.SI["Qs"]
        expected_count = 1 + 2 * p * s + 1 + 1 + 2 * Q
        model = FakeModel(list(range(1, expected_count + 1)))
        tool = JMAG.JMAG(deepcopy(self.config), deepcopy(self.spec))

        self.assertTrue(tool.pre_process_outer_rotor_spmsm(None, model, variant))
        layout = tool.outer_rotor_part_layout
        expected_motion = (
            [layout["rotor_core"]]
            + layout["magnets"]
            + [layout["sleeve"]]
        )
        self.assertEqual(
            model.set_list.GetSet("Motion_Region").part_ids,
            expected_motion,
        )
        self.assertNotIn("ShaftSet", model.set_list.sets)
        self.assertEqual(
            model.set_list.GetSet("MagnetSet").part_ids,
            layout["magnets"],
        )
        self.assertEqual(len(model.group_list.groups["Coils"]), 2 * Q)

    def test_material_assignment_uses_actual_inner_stator_name(self):
        template = SimpleNamespace(
            name="__OuterRotor",
            spec_input_dict={"Steel": "M19Gauge29", "Temperature": 35},
        )
        variant = SimpleNamespace(
            template=template,
            rotorCore=SimpleNamespace(name="OuterNotchedRotor"),
            stator_core=SimpleNamespace(name="InnerStatorCore"),
        )
        study = FakeStudy()
        tool = JMAG.JMAG(deepcopy(self.config), deepcopy(self.spec))

        tool.add_material(study, variant)

        self.assertEqual(
            study.assignments["InnerStatorCore"],
            "M-19 Steel Gauge-29",
        )
        self.assertNotIn("StatorCore", study.assignments)


if __name__ == "__main__":
    unittest.main()
