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


def derive_mm_r_so(GP, SI):
    """Derive the air-gap-facing radius of the inner stator."""
    GP["mm_r_so"].value = (
        GP["mm_r_ro"].value
        - GP["mm_d_ri"].value
        - GP["mm_d_pm"].value
        - GP["mm_d_sleeve"].value
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


def derive_split_ratio(GP, SI):
    if GP["mm_r_so"].value is None:
        derive_mm_r_so(GP, SI)
    GP["split_ratio"].value = GP["mm_r_so"].value / GP["mm_r_ro"].value
    return GP["split_ratio"].value


def derive_mm_r_ri(GP, SI):
    """Derive the inner radius of the complete rotating assembly."""
    GP["mm_r_ri"].value = (
        GP["mm_r_ro"].value
        - GP["mm_d_ri"].value
        - GP["mm_d_pm"].value
        - GP["mm_d_sleeve"].value
    )
    return GP["mm_r_ri"].value


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
                    lambda GP, SI: None,
                ),
                "mm_d_rp": acmop_parameter(
                    "free",
                    "inter_polar_iron_thickness",
                    None,
                    [None, None],
                    lambda GP, SI: None,
                ),
                "deg_alpha_rs": acmop_parameter(
                    "free" if SI["no_segmented_magnets"] != 1 else "fixed",
                    "magnet_segment_span_angle",
                    None,
                    [None, None],
                    lambda GP, SI: None,
                ),
                "mm_d_rs": acmop_parameter(
                    "free" if SI["no_segmented_magnets"] != 1 else "fixed",
                    "inter_segment_iron_thickness",
                    None,
                    [None, None],
                    lambda GP, SI: None,
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
        sleeve_depth = SI["mm_sleeve_length"]
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

        alpha_rm_over_alpha_rp = 1.0 if p >= 2 else 0.75
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
        stator_tooth_open_depth = max(0.05, 0.1 * stator_yoke_depth)
        stator_tooth_tip_depth = 1.5 * stator_tooth_open_depth
        stator_tooth_depth = (
            stator_outer_radius
            - stator_bore_radius
            - stator_yoke_depth
            - stator_tooth_tip_depth
        )
        if stator_tooth_depth <= 0:
            raise ValueError(
                "The inner-stator radial stack is invalid. Reduce the stator bore, "
                "rotor yoke, magnet, sleeve, or air-gap depth."
            )

        stator_tooth_width = (
            Bg
            * np.pi
            * stator_air_gap_diameter_m
            / (Bst * Q)
            * 1e3
        )

        GP["mm_r_ro"].value = rotor_outer_radius
        GP["mm_d_mech_air_gap"].value = mechanical_air_gap
        GP["mm_d_sleeve"].value = sleeve_depth
        GP["deg_alpha_st"].value = 360.0 / Q - 2.0
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

        GP["deg_alpha_rm"].value = 360.0 / (2.0 * p)
        GP["mm_d_rp"].value = magnet_depth
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

        return {
            "mm_r_ro": around(GP["mm_r_ro"].value, 0.9, 1.1),
            "mm_d_mech_air_gap": around(
                GP["mm_d_mech_air_gap"].value,
                0.7,
                1.5,
            ),
            "mm_d_sleeve": around(GP["mm_d_sleeve"].value, 0.5, 2.0),
            "split_ratio": [0.4, 0.95],
            "deg_alpha_st": [0.35 * 360.0 / Q, 0.95 * 360.0 / Q],
            "mm_w_st": around(GP["mm_w_st"].value, 0.7, 1.2),
            "mm_d_sto": around(GP["mm_d_sto"].value, 0.5, 2.0, 0.05),
            "deg_alpha_sto": [0.0, 0.5 * 360.0 / Q],
            "mm_d_stt": around(GP["mm_d_stt"].value, 0.7, 1.3),
            "mm_r_si": around(GP["mm_r_si"].value, 0.7, 1.2),
            "mm_r_so": around(GP["mm_r_so"].value, 0.8, 1.05),
            "mm_d_sy": around(GP["mm_d_sy"].value, 0.8, 1.2),
            "mm_d_st": around(GP["mm_d_st"].value, 0.8, 1.15),
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
