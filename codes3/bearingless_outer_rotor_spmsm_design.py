import logging
from collections import OrderedDict

import numpy as np

import CrossSectOuterNotchedRotor
import CrossSectOuterRotorStator
import Location2D
import pyrhonen_procedure_as_function
import inner_rotor_motor
import winding_layout
from utility import acmop_parameter


RADIAL_CONSTRAINT_MODES = {
    "fixed_outer_radius",
    "fixed_stator_bore",
    "fixed_envelope",
}
OPTIMIZATION_PARAMETERIZATIONS = {
    "legacy",
    "stator_slot_only",
    "stator_slot_and_notched_rotor",
}


def get_radial_constraint_mode(SI):
    mode = SI.get("radial_constraint_mode", "fixed_outer_radius")
    if mode not in RADIAL_CONSTRAINT_MODES:
        raise ValueError(
            "radial_constraint_mode must be one of %s, got %r."
            % (sorted(RADIAL_CONSTRAINT_MODES), mode)
        )
    return mode


def get_optimization_parameterization(SI):
    parameterization = SI.get("optimization_parameterization", "legacy")
    if parameterization not in OPTIMIZATION_PARAMETERIZATIONS:
        raise ValueError(
            "optimization_parameterization must be one of %s, got %r."
            % (sorted(OPTIMIZATION_PARAMETERIZATIONS), parameterization)
        )
    return parameterization


def uses_sleeve(SI):
    value = SI.get("use_sleeve", SI.get("mm_sleeve_length", 0.0) > 0.0)
    if not isinstance(value, bool):
        raise TypeError("use_sleeve must be a JSON boolean.")
    return value


def get_effective_sleeve_depth(GP, SI):
    return GP["mm_d_sleeve"].value if uses_sleeve(SI) else 0.0


def get_magnet_pole_arc_ratio(SI):
    ratio = float(SI.get("rotor_magnet_pole_arc_ratio", 1.0))
    if not 0.0 < ratio <= 1.0:
        raise ValueError("rotor_magnet_pole_arc_ratio must be in (0, 1].")
    return ratio


def get_inter_polar_iron_depth_ratio(SI):
    ratio = float(SI.get("inter_polar_iron_depth_ratio", 1.0))
    if not 0.0 <= ratio <= 1.0:
        raise ValueError("inter_polar_iron_depth_ratio must be in [0, 1].")
    return ratio


def _slot_opening_width_from_tooth_span(radius, tooth_span_deg, Q):
    slot_opening_angle = 2.0 * np.pi / Q - np.deg2rad(tooth_span_deg)
    return 2.0 * radius * np.sin(0.5 * slot_opening_angle)


def _tooth_span_from_slot_opening_width(radius, slot_opening_width, Q):
    maximum_width = 2.0 * radius * np.sin(np.pi / Q)
    if not 0.0 < slot_opening_width < maximum_width:
        raise ValueError(
            "The circumferential stator slot opening must be in (0, %.6g) mm."
            % maximum_width
        )
    slot_opening_angle = 2.0 * np.arcsin(
        slot_opening_width / (2.0 * radius)
    )
    return np.rad2deg(2.0 * np.pi / Q - slot_opening_angle)


def _maximum_tooth_depth(radius, tooth_width, Q, clearance=0.0):
    half_tooth_width = 0.5 * tooth_width
    if half_tooth_width >= radius:
        return 0.0
    return (
        np.sqrt(radius**2 - half_tooth_width**2)
        - half_tooth_width / np.tan(np.pi / Q)
        - clearance
    )


def _maximum_tooth_width(radius, tooth_depth, Q, clearance=0.0):
    effective_depth = tooth_depth + clearance
    half_slot_pitch = np.pi / Q
    radicand = radius**2 - (
        effective_depth * np.sin(half_slot_pitch)
    ) ** 2
    if radicand <= 0.0:
        return 0.0
    maximum_half_width = np.sin(half_slot_pitch) * (
        np.sqrt(radicand)
        - effective_depth * np.cos(half_slot_pitch)
    )
    return max(0.0, 2.0 * maximum_half_width)


def _ensure_mm_d_stt(GP, SI):
    if GP["mm_d_stt"].value is None:
        GP["mm_d_stt"].value = GP["mm_d_stt"].calc(GP, SI)
    return GP["mm_d_stt"].value


def derive_mm_r_ro(GP, SI):
    """Derive rotor outer radius when the stator bore is the radial anchor."""
    if GP["mm_r_so"].value is None:
        derive_mm_r_so(GP, SI)
    GP["mm_r_ro"].value = (
        GP["mm_r_so"].value
        + GP["mm_d_mech_air_gap"].value
        + get_effective_sleeve_depth(GP, SI)
        + GP["mm_d_pm"].value
        + GP["mm_d_ri"].value
    )
    return GP["mm_r_ro"].value


def derive_mm_r_so(GP, SI):
    """Derive the air-gap-facing radius of the inner stator."""
    mode = get_radial_constraint_mode(SI)
    if mode == "fixed_stator_bore":
        GP["mm_r_so"].value = (
            GP["mm_r_si"].value
            + GP["mm_d_sy"].value
            + GP["mm_d_st"].value
            + _ensure_mm_d_stt(GP, SI)
        )
    else:
        GP["mm_r_so"].value = (
            GP["mm_r_ro"].value
            - GP["mm_d_ri"].value
            - GP["mm_d_pm"].value
            - get_effective_sleeve_depth(GP, SI)
            - GP["mm_d_mech_air_gap"].value
        )
    return GP["mm_r_so"].value


def derive_mm_r_si(GP, SI):
    """Derive the inner bore radius of the inner stator."""
    if GP["mm_r_so"].value is None:
        derive_mm_r_so(GP, SI)
    GP["mm_r_si"].value = (
        GP["mm_r_so"].value
        - GP["mm_d_stt"].value
        - GP["mm_d_st"].value
        - GP["mm_d_sy"].value
    )
    return GP["mm_r_si"].value


def derive_mm_d_sy(GP, SI):
    """Close a fixed inner/outer envelope through the stator yoke depth."""
    GP["mm_d_sy"].value = (
        GP["mm_r_ro"].value
        - GP["mm_r_si"].value
        - _ensure_mm_d_stt(GP, SI)
        - GP["mm_d_st"].value
        - GP["mm_d_mech_air_gap"].value
        - get_effective_sleeve_depth(GP, SI)
        - GP["mm_d_pm"].value
        - GP["mm_d_ri"].value
    )
    return GP["mm_d_sy"].value


def derive_deg_alpha_rm(GP, SI):
    GP["deg_alpha_rm"].value = (
        get_magnet_pole_arc_ratio(SI) * 360.0 / (2.0 * SI["p"])
    )
    return GP["deg_alpha_rm"].value


def derive_mm_d_rp(GP, SI):
    GP["mm_d_rp"].value = (
        get_inter_polar_iron_depth_ratio(SI) * GP["mm_d_pm"].value
    )
    return GP["mm_d_rp"].value


def derive_deg_alpha_rs(GP, SI):
    if GP["deg_alpha_rm"].value is None:
        derive_deg_alpha_rm(GP, SI)
    GP["deg_alpha_rs"].value = (
        GP["deg_alpha_rm"].value / SI["no_segmented_magnets"]
    )
    return GP["deg_alpha_rs"].value


def derive_mm_d_rs(GP, SI):
    GP["mm_d_rs"].value = 0.2 * GP["mm_d_pm"].value
    return GP["mm_d_rs"].value


def derive_split_ratio(GP, SI):
    if GP["mm_r_ro"].value is None:
        derive_mm_r_ro(GP, SI)
    if GP["mm_r_so"].value is None:
        derive_mm_r_so(GP, SI)
    GP["split_ratio"].value = GP["mm_r_so"].value / GP["mm_r_ro"].value
    return GP["split_ratio"].value


def derive_mm_r_ri(GP, SI):
    """Derive the inner radius of the complete rotating assembly."""
    if GP["mm_r_so"].value is None:
        derive_mm_r_so(GP, SI)
    GP["mm_r_ri"].value = (
        GP["mm_r_so"].value + GP["mm_d_mech_air_gap"].value
    )
    return GP["mm_r_ri"].value


def validate_radial_geometry(GP, SI):
    expected_r_so_from_rotor = (
        GP["mm_r_ro"].value
        - GP["mm_d_ri"].value
        - GP["mm_d_pm"].value
        - get_effective_sleeve_depth(GP, SI)
        - GP["mm_d_mech_air_gap"].value
    )
    expected_r_so_from_stator = (
        GP["mm_r_si"].value
        + GP["mm_d_sy"].value
        + GP["mm_d_st"].value
        + GP["mm_d_stt"].value
    )
    if not np.isclose(GP["mm_r_so"].value, expected_r_so_from_rotor):
        raise ValueError("The rotor-side radial dimensions do not close.")
    if not np.isclose(GP["mm_r_so"].value, expected_r_so_from_stator):
        raise ValueError("The stator-side radial dimensions do not close.")
    if not (
        0.0
        < GP["mm_r_si"].value
        < GP["mm_r_so"].value
        < GP["mm_r_ri"].value
        < GP["mm_r_ro"].value
    ):
        raise ValueError(
            "Expected mm_r_si < mm_r_so < mm_r_ri < mm_r_ro for an outer rotor."
        )
    if GP["mm_d_sy"].value <= 0.0:
        raise ValueError("The radial closure produces a non-positive stator yoke.")
    if uses_sleeve(SI):
        if GP["mm_d_sleeve"].value <= 0.0:
            raise ValueError("A sleeve-enabled design requires positive sleeve depth.")
    elif not np.isclose(GP["mm_d_sleeve"].value, 0.0):
        raise ValueError("A sleeve-disabled design must use zero sleeve depth.")


class bearingless_outer_rotor_spmsm_template(
    inner_rotor_motor.template_machine_as_numbers
):
    """Numerical template for a surface-PM outer-rotor machine."""

    def __init__(self, fea_config_dict, spec_input_dict):
        super(bearingless_outer_rotor_spmsm_template, self).__init__(
            fea_config_dict,
            spec_input_dict,
        )

        self.machine_type = "OuterRotorSPMSM"
        self.name = "__OuterRotorSPMSM"

        GP = self.d["GP"]
        EX = self.d["EX"]
        SI = self.SI

        # Replace the topology-dependent parameters inherited from the
        # inner-rotor base class. Assignment preserves OrderedDict order.
        GP["mm_r_ro"] = acmop_parameter(
            "fixed",
            "outer_rotor_outer_radius",
            None,
            [None, None],
            derive_mm_r_ro,
        )
        GP["split_ratio"] = acmop_parameter(
            "derived",
            "split_ratio_r_os_slash_r_or",
            None,
            [None, None],
            derive_split_ratio,
        )
        GP["mm_r_si"] = acmop_parameter(
            "derived",
            "inner_stator_bore_radius",
            None,
            [None, None],
            derive_mm_r_si,
        )
        GP["mm_r_so"] = acmop_parameter(
            "derived",
            "inner_stator_outer_radius",
            None,
            [None, None],
            derive_mm_r_so,
        )
        GP["mm_d_sy"] = acmop_parameter(
            "fixed",
            "stator_yoke_depth",
            None,
            [None, None],
            derive_mm_d_sy,
        )

        child_gp = OrderedDict(
            {
                "mm_d_pm": acmop_parameter(
                    "free", "magnet_depth", None, [None, None], lambda GP, SI: None
                ),
                "mm_d_ri": acmop_parameter(
                    "free",
                    "outer_rotor_yoke_depth",
                    None,
                    [None, None],
                    lambda GP, SI: None,
                ),
                "deg_alpha_rm": acmop_parameter(
                    "fixed",
                    "magnet_pole_span_angle",
                    None,
                    [None, None],
                    derive_deg_alpha_rm,
                ),
                "mm_d_rp": acmop_parameter(
                    "free",
                    "inter_polar_iron_thickness",
                    None,
                    [None, None],
                    derive_mm_d_rp,
                ),
                "deg_alpha_rs": acmop_parameter(
                    "free" if SI["no_segmented_magnets"] != 1 else "fixed",
                    "magnet_segment_span_angle",
                    None,
                    [None, None],
                    derive_deg_alpha_rs,
                ),
                "mm_d_rs": acmop_parameter(
                    "free" if SI["no_segmented_magnets"] != 1 else "fixed",
                    "inter_segment_iron_thickness",
                    None,
                    [None, None],
                    derive_mm_d_rs,
                ),
                "mm_r_ri": acmop_parameter(
                    "derived",
                    "outer_rotor_inner_radius",
                    None,
                    [None, None],
                    derive_mm_r_ri,
                ),
            }
        )
        GP.update(child_gp)

        self.PracticalInitialDesign(fea_config_dict, SI, GP, EX)
        self.set_gp_values_and_types_based_on_fea_config(fea_config_dict)
        self.configure_radial_parameter_types(GP, SI)
        self.configure_optimization_parameter_types(GP, SI)
        validate_radial_geometry(GP, SI)

        self.original_template_neighbor_bounds = self.get_template_neighbor_bounds()
        self.bounds_denorm = self.define_search_space(
            GP,
            self.original_template_neighbor_bounds,
        )

        self.get_other_properties_after_geometric_parameters_are_initialized(GP, SI)
        EX["BeariW_zQ"] = EX["DriveW_zQ"]
        EX["BeariW_CurrentAmp"] = fea_config_dict[
            "circuit.SUSPENSION_CURRENT_RATIO"
        ] * (
            EX["DriveW_CurrentAmp"]
            / fea_config_dict["circuit.TORQUE_CURRENT_RATIO"]
        )
        EX["BeariW_Freq"] = EX["DriveW_Freq"]
        EX["BeariW_Rs"] = EX["DriveW_Rs"]
        EX["BeariW_poles"] = SI["ps"] * 2
        EX["slot_current_utilizing_ratio"] = (
            fea_config_dict["circuit.SUSPENSION_CURRENT_RATIO"]
            + fea_config_dict["circuit.TORQUE_CURRENT_RATIO"]
        )

    @staticmethod
    def configure_radial_parameter_types(GP, SI):
        mode = get_radial_constraint_mode(SI)
        GP["mm_r_ro"].type = (
            "derived" if mode == "fixed_stator_bore" else "fixed"
        )
        GP["mm_r_si"].type = (
            "derived" if mode == "fixed_outer_radius" else "fixed"
        )
        GP["mm_d_sy"].type = (
            "derived" if mode == "fixed_envelope" else "fixed"
        )
        GP["mm_r_so"].type = "derived"
        GP["mm_r_ri"].type = "derived"
        GP["split_ratio"].type = "derived"

        if not uses_sleeve(SI):
            GP["mm_d_sleeve"].type = "fixed"
            GP["mm_d_sleeve"].value = 0.0

    @staticmethod
    def configure_optimization_parameter_types(GP, SI):
        parameterization = get_optimization_parameterization(SI)
        if parameterization == "legacy":
            return

        for key in ("deg_alpha_st", "mm_w_st", "mm_d_sto", "mm_d_st"):
            GP[key].type = "free"

        for key in (
            "mm_d_mech_air_gap",
            "mm_d_sleeve",
            "mm_d_pm",
            "mm_d_ri",
        ):
            GP[key].type = "fixed"

        rotor_parameter_type = (
            "free"
            if parameterization == "stator_slot_and_notched_rotor"
            else "derived"
        )
        GP["deg_alpha_rm"].type = rotor_parameter_type
        GP["mm_d_rp"].type = rotor_parameter_type
        GP["deg_alpha_rs"].type = "derived"
        GP["mm_d_rs"].type = (
            "fixed" if SI["no_segmented_magnets"] == 1 else "derived"
        )
        if SI["no_segmented_magnets"] == 1:
            GP["mm_d_rs"].value = 0.0

    @staticmethod
    def _machine_outer_diameter(SI):
        for key in (
            "mm_rotor_outer_diameter",
            "mm_machine_outer_diameter",
            "mm_stator_outer_diameter",
        ):
            if key in SI:
                return SI[key]
        raise KeyError(
            "Outer-rotor design requires mm_rotor_outer_diameter, "
            "mm_machine_outer_diameter, or mm_stator_outer_diameter."
        )

    def PracticalInitialDesign(self, fea_config_dict, SI, GP, EX):
        Bg = SI["guess_air_gap_flux_density_Bg"]
        Bst = SI["guess_stator_tooth_flux_density_Bst"]
        Bsy = SI["guess_stator_yoke_flux_density_Bsy"]
        p = SI["p"]
        Q = SI["Qs"]

        rotor_outer_radius = 0.5 * self._machine_outer_diameter(SI)
        rotor_yoke_depth = SI["mm_d_ri"]
        magnet_depth = SI["mm_d_pm"]
        sleeve_depth = SI["mm_sleeve_length"] if uses_sleeve(SI) else 0.0
        mechanical_air_gap = SI["minimum_mechanical_air_gap_length_mm"]

        stator_outer_radius = (
            rotor_outer_radius
            - rotor_yoke_depth
            - magnet_depth
            - sleeve_depth
            - mechanical_air_gap
        )
        if stator_outer_radius <= 0:
            raise ValueError("The rotor radial stack leaves no room for the stator.")

        alpha_rm_over_alpha_rp = get_magnet_pole_arc_ratio(SI)
        stator_air_gap_diameter_m = 2.0 * stator_outer_radius * 1e-3
        stator_yoke_depth = (
            Bg
            * np.pi
            * stator_air_gap_diameter_m
            * alpha_rm_over_alpha_rp
            / (2.0 * Bsy * 2.0 * p)
            * 1e3
        )
        stator_bore_radius = SI.get(
            "mm_stator_inner_radius",
            SI.get("mm_radius_shaft", 0.2 * stator_outer_radius),
        )
        stator_tooth_open_depth = float(
            SI.get("mm_stator_slot_opening_depth", 2.0)
        )
        if stator_tooth_open_depth <= 0:
            raise ValueError("mm_stator_slot_opening_depth must be positive.")
        stator_tooth_tip_depth = 1.5 * stator_tooth_open_depth
        stator_tooth_depth = (
            stator_outer_radius
            - stator_bore_radius
            - stator_yoke_depth
            - stator_tooth_tip_depth
        )
        if stator_tooth_depth <= 0:
            raise ValueError(
                "The inner-stator radial stack is invalid: "
                "r_so=%.6g mm, r_si=%.6g mm, d_sy=%.6g mm, "
                "d_stt=%.6g mm, leaving d_st=%.6g mm. Reduce the stator "
                "bore, rotor yoke, magnet, sleeve, air-gap, or assumed flux "
                "loading."
                % (
                    stator_outer_radius,
                    stator_bore_radius,
                    stator_yoke_depth,
                    stator_tooth_tip_depth,
                    stator_tooth_depth,
                )
            )

        slot_flux_angle = p * np.pi / Q
        slot_flux_factor = abs(np.sin(slot_flux_angle) / slot_flux_angle)
        stator_iron_fill_factor = float(
            SI.get("lamination_stacking_factor_kFe", 1.0)
        )
        if not 0.0 < stator_iron_fill_factor <= 1.0:
            raise ValueError("lamination_stacking_factor_kFe must be in (0, 1].")
        stator_tooth_width = (
            Bg
            * np.pi
            * stator_air_gap_diameter_m
            * slot_flux_factor
            / (Bst * Q * stator_iron_fill_factor)
            * 1e3
        )

        configured_slot_opening_width = SI.get("mm_stator_slot_opening_width")
        if configured_slot_opening_width is None:
            stator_tooth_span = 360.0 / Q - 2.0
            slot_opening_width = _slot_opening_width_from_tooth_span(
                stator_outer_radius,
                stator_tooth_span,
                Q,
            )
        else:
            slot_opening_width = float(configured_slot_opening_width)
            stator_tooth_span = _tooth_span_from_slot_opening_width(
                stator_outer_radius,
                slot_opening_width,
                Q,
            )

        GP["mm_r_ro"].value = rotor_outer_radius
        GP["mm_d_mech_air_gap"].value = mechanical_air_gap
        GP["mm_d_sleeve"].value = sleeve_depth
        GP["deg_alpha_st"].value = stator_tooth_span
        GP["mm_w_st"].value = stator_tooth_width
        GP["mm_d_sto"].value = stator_tooth_open_depth
        GP["deg_alpha_sto"].value = 0.5 * GP["deg_alpha_st"].value
        GP["mm_d_stt"].value = stator_tooth_tip_depth
        GP["mm_d_sy"].value = stator_yoke_depth
        GP["mm_d_st"].value = stator_tooth_depth
        GP["mm_d_pm"].value = magnet_depth
        GP["mm_d_ri"].value = rotor_yoke_depth
        GP["mm_r_so"].value = stator_outer_radius
        GP["mm_r_si"].value = stator_bore_radius
        GP["split_ratio"].value = stator_outer_radius / rotor_outer_radius
        GP["mm_r_ri"].value = stator_outer_radius + mechanical_air_gap

        EX["stator_slot_flux_factor"] = slot_flux_factor
        EX["mm_stator_slot_opening_width"] = slot_opening_width

        GP["deg_alpha_rm"].value = (
            get_magnet_pole_arc_ratio(SI) * 360.0 / (2.0 * p)
        )
        GP["mm_d_rp"].value = (
            get_inter_polar_iron_depth_ratio(SI) * magnet_depth
        )
        GP["deg_alpha_rs"].value = (
            GP["deg_alpha_rm"].value / SI["no_segmented_magnets"]
        )
        GP["mm_d_rs"].value = (
            0.0
            if SI["no_segmented_magnets"] == 1
            else 0.2 * magnet_depth
        )

        radius_slot_outer = stator_outer_radius - stator_tooth_tip_depth
        radius_yoke_outer = stator_bore_radius + stator_yoke_depth
        EX["stator_slot_area"] = 1e-6 * (
            np.pi
            / Q
            * (radius_slot_outer**2 - radius_yoke_outer**2)
            - stator_tooth_width * stator_tooth_depth
        )
        mean_slot_radius = 0.5 * (radius_slot_outer + radius_yoke_outer)
        slot_pitch = 2.0 * np.pi * mean_slot_radius / Q
        EX["end_winding_length_Lew"] = 1e-3 * (
            0.5 * np.pi * (slot_pitch + stator_tooth_width)
            + slot_pitch * 1.8 * (SI["coil_pitch_y"] - 1)
        )

    def get_other_properties_after_geometric_parameters_are_initialized(
        self,
        GP,
        SI,
        specified_mm_stack_length=None,
    ):
        EX = self.d["EX"]
        if specified_mm_stack_length is None:
            specified_mm_stack_length = SI.get("mm_stack_length")
        if "Wrap_Around" in SI:
            wily = winding_layout.winding_layout_v2(
                SI["DPNV_or_SEPA"],
                SI["Qs"],
                SI["p"],
                SI["ps"],
                SI["coil_pitch_y"],
                m=SI["m"],
                Wrap_Around=SI["Wrap_Around"],
            )
        else:
            wily = winding_layout.winding_layout_v2(
                SI["DPNV_or_SEPA"],
                SI["Qs"],
                SI["p"],
                SI["ps"],
                SI["coil_pitch_y"],
                m=SI["m"],
            )
        EX["wily"] = wily

        if specified_mm_stack_length is None:
            EX["mm_template_stack_length"] = (
                pyrhonen_procedure_as_function.get_mm_template_stack_length(
                    SI,
                    GP["mm_r_so"].value * 1e-3,
                )
            )
        else:
            if specified_mm_stack_length <= 0:
                raise ValueError("mm_stack_length must be positive.")
            EX["mm_template_stack_length"] = specified_mm_stack_length

        EX["mm_mechanical_air_gap_length"] = SI[
            "minimum_mechanical_air_gap_length_mm"
        ]
        EX["Js"] = SI["Js"]
        EX["WindingFill"] = SI["WindingFill"]

        magnet_inner_radius = (
            GP["mm_r_ro"].value
            - GP["mm_d_ri"].value
            - GP["mm_d_pm"].value
        )
        EX["DriveW_zQ"] = pyrhonen_procedure_as_function.get_zQ(
            SI,
            wily,
            magnet_inner_radius * 2e-3,
            GP["mm_r_so"].value * 2e-3,
            specified_mm_stack_length=EX["mm_template_stack_length"],
        )
        EX["DriveW_CurrentAmp"] = (
            np.sqrt(2.0)
            * pyrhonen_procedure_as_function.get_stator_phase_current_rms(SI)
        )
        EX["DriveW_Freq"] = SI["ExcitationFreqSimulated"]
        EX["DriveW_Rs"] = 1.0
        EX["DriveW_poles"] = SI["p"] * 2
        return EX

    def get_template_neighbor_bounds(self):
        GP = self.d["GP"]
        Q = self.SI["Qs"]
        p = self.SI["p"]
        s = self.SI["no_segmented_magnets"]

        def around(value, lower_factor=0.8, upper_factor=1.2, floor=1e-3):
            return [max(floor, lower_factor * value), upper_factor * value]

        slot_opening_width = _slot_opening_width_from_tooth_span(
            GP["mm_r_so"].value,
            GP["deg_alpha_st"].value,
            Q,
        )
        slot_opening_width_bounds = around(slot_opening_width, 0.8, 1.2)
        tooth_span_bounds = sorted(
            _tooth_span_from_slot_opening_width(
                GP["mm_r_so"].value,
                width,
                Q,
            )
            for width in slot_opening_width_bounds
        )

        tooth_width_bounds = around(GP["mm_w_st"].value, 0.7, 1.2)
        tooth_open_depth_bounds = around(
            GP["mm_d_sto"].value,
            0.8,
            1.2,
            0.05,
        )
        maximum_tooth_tip_depth = 1.5 * tooth_open_depth_bounds[1]
        minimum_slot_opening_radius = (
            GP["mm_r_so"].value - maximum_tooth_tip_depth
        )
        geometry_clearance = 1e-3
        maximum_tooth_depth_from_radial_closure = (
            GP["mm_r_so"].value
            - GP["mm_r_si"].value
            - maximum_tooth_tip_depth
            - geometry_clearance
        )
        tooth_depth_bounds = around(GP["mm_d_st"].value, 0.8, 1.15)
        maximum_tooth_depth_at_initial_width = _maximum_tooth_depth(
            minimum_slot_opening_radius,
            GP["mm_w_st"].value,
            Q,
            geometry_clearance,
        )
        tooth_depth_bounds[1] = min(
            tooth_depth_bounds[1],
            maximum_tooth_depth_at_initial_width,
            maximum_tooth_depth_from_radial_closure,
        )
        if tooth_depth_bounds[1] + 1e-9 < GP["mm_d_st"].value:
            raise ValueError(
                "The initial stator tooth cannot remain valid over the requested "
                "slot-opening-depth range. Reduce mm_stator_slot_opening_depth, "
                "mm_w_st, or increase mm_stator_inner_radius."
            )

        maximum_width_at_depth_upper = _maximum_tooth_width(
            minimum_slot_opening_radius,
            tooth_depth_bounds[1],
            Q,
            geometry_clearance,
        )
        tooth_width_bounds[1] = min(
            tooth_width_bounds[1],
            maximum_width_at_depth_upper,
        )
        if tooth_width_bounds[1] + 1e-9 < GP["mm_w_st"].value:
            raise ValueError(
                "The initial stator tooth width is incompatible with the tooth-depth "
                "optimization range. Reduce the upper tooth-depth factor."
            )

        return {
            "mm_r_ro": around(GP["mm_r_ro"].value, 0.9, 1.1),
            "mm_d_mech_air_gap": around(
                GP["mm_d_mech_air_gap"].value,
                0.7,
                1.5,
            ),
            "mm_d_sleeve": (
                around(GP["mm_d_sleeve"].value, 0.5, 2.0)
                if uses_sleeve(self.SI)
                else [0.0, 0.0]
            ),
            "split_ratio": [0.4, 0.95],
            "deg_alpha_st": tooth_span_bounds,
            "mm_w_st": tooth_width_bounds,
            "mm_d_sto": tooth_open_depth_bounds,
            "deg_alpha_sto": [0.0, 0.5 * 360.0 / Q],
            "mm_d_stt": around(GP["mm_d_stt"].value, 0.7, 1.3),
            "mm_r_si": around(GP["mm_r_si"].value, 0.7, 1.2),
            "mm_r_so": around(GP["mm_r_so"].value, 0.8, 1.05),
            "mm_d_sy": around(GP["mm_d_sy"].value, 0.8, 1.2),
            "mm_d_st": tooth_depth_bounds,
            "mm_d_pm": around(GP["mm_d_pm"].value, 0.6, 1.4),
            "mm_d_ri": around(GP["mm_d_ri"].value, 0.7, 1.3),
            "deg_alpha_rm": [0.6 * 360.0 / (2 * p), 360.0 / (2 * p)],
            "mm_d_rp": [max(0.05, 0.4 * GP["mm_d_pm"].value), GP["mm_d_pm"].value],
            "deg_alpha_rs": [
                0.5 * GP["deg_alpha_rm"].value / s,
                GP["deg_alpha_rm"].value / s,
            ],
            "mm_d_rs": [
                0.0 if s == 1 else max(0.01, 0.05 * GP["mm_d_pm"].value),
                GP["mm_d_pm"].value,
            ],
            "mm_r_ri": around(GP["mm_r_ri"].value, 0.8, 1.05),
        }

    def get_rotor_volume(self, stack_length=None):
        GP = self.d["GP"]
        if stack_length is None:
            stack_length = self.d["EX"]["mm_template_stack_length"]
        outer_radius_m = GP["mm_r_ro"].value * 1e-3
        inner_radius_m = GP["mm_r_ri"].value * 1e-3
        return (
            np.pi
            * (outer_radius_m**2 - inner_radius_m**2)
            * stack_length
            * 1e-3
        )


class bearingless_outer_rotor_spmsm_design_variant(
    inner_rotor_motor.variant_machine_as_objects
):
    """Object-based outer-rotor design assembled from cross sections."""

    def __init__(self, template=None, x_denorm=None, counter=None, counter_loop=None):
        if x_denorm is None:
            raise ValueError("x_denorm is required.")

        super(bearingless_outer_rotor_spmsm_design_variant, self).__init__(
            template,
            x_denorm,
            counter,
            counter_loop,
        )
        self.x_denorm = x_denorm
        self.name = "ind%s" % counter
        if counter_loop is not None and counter_loop > 1:
            self.name += "-redo%s" % counter_loop

        GP = self.template.d["GP"]
        SI = self.template.SI
        self.check_invalid_design(GP, SI)
        validate_radial_geometry(GP, SI)

        origin = Location2D.Location2D(anchor_xy=[0, 0], deg_theta=0)
        self.rotorCore = CrossSectOuterNotchedRotor.CrossSectOuterNotchedRotor(
            name="OuterNotchedRotor",
            mm_d_pm=GP["mm_d_pm"].value,
            deg_alpha_rm=GP["deg_alpha_rm"].value,
            deg_alpha_rs=GP["deg_alpha_rs"].value,
            mm_d_ri=GP["mm_d_ri"].value,
            mm_r_ro=GP["mm_r_ro"].value,
            mm_d_rp=GP["mm_d_rp"].value,
            mm_d_rs=GP["mm_d_rs"].value,
            p=SI["p"],
            s=SI["no_segmented_magnets"],
            location=origin,
        )
        self.rotorMagnet = (
            CrossSectOuterNotchedRotor.CrossSectOuterNotchedMagnet(
                name="OuterRotorMagnet",
                notched_rotor=self.rotorCore,
            )
        )
        self.sleeve = None
        if uses_sleeve(SI):
            self.sleeve = CrossSectOuterNotchedRotor.CrossSectOuterRotorSleeve(
                name="OuterRotorSleeve",
                notched_magnet=self.rotorMagnet,
                d_sleeve=GP["mm_d_sleeve"].value,
            )
        self.stator_core = CrossSectOuterRotorStator.CrossSectOuterRotorStator(
            name="InnerStatorCore",
            deg_alpha_st=GP["deg_alpha_st"].value,
            deg_alpha_sto=GP["deg_alpha_sto"].value,
            mm_r_so=GP["mm_r_so"].value,
            mm_d_sto=GP["mm_d_sto"].value,
            mm_d_stt=GP["mm_d_stt"].value,
            mm_d_st=GP["mm_d_st"].value,
            mm_d_sy=GP["mm_d_sy"].value,
            mm_w_st=GP["mm_w_st"].value,
            mm_r_st=0.0,
            mm_r_sf=0.0,
            mm_r_sb=0.0,
            Q=SI["Qs"],
            location=origin,
        )
        self.coils = CrossSectOuterRotorStator.CrossSectOuterRotorStatorWinding(
            name="InnerStatorCoils",
            stator_core=self.stator_core,
        )
        self.stator_core._build_points()

        # An outer-rotor model has no rotating center shaft. JMAG needs a
        # topology-specific motion-region implementation before FEA dispatch.
        self.shaft = None

        self.update_mechanical_parameters()
        deg_pole_span = 180.0 / SI["p"]
        winding_axis = self.template.d["EX"][
            "wily"
        ].deg_winding_U_phase_phase_axis_angle
        self.InitialRotationAngle = (
            (deg_pole_span - GP["deg_alpha_rm"].value) * 0.5
            - deg_pole_span * 0.5
            + winding_axis
            + deg_pole_span
        )
        self.boolCustomizedCircuit = False

    @staticmethod
    def check_invalid_design(GP, SI):
        if GP["mm_d_rp"].value > GP["mm_d_pm"].value:
            logging.getLogger(__name__).warning(
                "mm_d_rp exceeds mm_d_pm; clamping it to the magnet depth."
            )
            GP["mm_d_rp"].value = GP["mm_d_pm"].value

        if GP["mm_d_rs"].value > GP["mm_d_pm"].value:
            logging.getLogger(__name__).warning(
                "mm_d_rs exceeds mm_d_pm; clamping it to the magnet depth."
            )
            GP["mm_d_rs"].value = GP["mm_d_pm"].value

        segment_count = SI["no_segmented_magnets"]
        if segment_count == 1:
            GP["deg_alpha_rs"].value = GP["deg_alpha_rm"].value
            GP["mm_d_rs"].value = 0.0
        else:
            maximum_segment_span = GP["deg_alpha_rm"].value / segment_count
            if GP["deg_alpha_rs"].value > maximum_segment_span:
                logging.getLogger(__name__).warning(
                    "deg_alpha_rs exceeds deg_alpha_rm/s; clamping it."
                )
                GP["deg_alpha_rs"].value = maximum_segment_span
            if GP["mm_d_rs"].value <= 0:
                GP["mm_d_rs"].value = min(
                    GP["mm_d_pm"].value,
                    max(0.01, 0.05 * GP["mm_d_pm"].value),
                )


# Familiar names when this module is selected directly by the framework.
bearingless_spmsm_template = bearingless_outer_rotor_spmsm_template
bearingless_spmsm_design_variant = bearingless_outer_rotor_spmsm_design_variant
