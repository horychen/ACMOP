import numpy as np


def _polar(radius, angle):
    return [radius * np.cos(angle), radius * np.sin(angle)]


def _rotate(point, angle):
    cos_angle = np.cos(angle)
    sin_angle = np.sin(angle)
    return [
        point[0] * cos_angle - point[1] * sin_angle,
        point[0] * sin_angle + point[1] * cos_angle,
    ]


def _mirror_x(point):
    return [point[0], -point[1]]


def _polygon_area(points):
    x = np.asarray([point[0] for point in points], dtype=float)
    y = np.asarray([point[1] for point in points], dtype=float)
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


class CrossSectOuterRotorStator:
    """Inner stator cross-section for an outer-rotor radial-flux machine.

    The teeth face radially outward. ``mm_r_so`` is the air-gap-facing
    stator radius, while the yoke is located at the smaller radii.
    """

    def __init__(
        self,
        name="StatorCore",
        color="#BAFD01",
        deg_alpha_st=12.0,
        deg_alpha_sto=6.0,
        mm_r_so=60.0,
        mm_d_sto=3.0,
        mm_d_stt=5.0,
        mm_d_st=20.0,
        mm_d_sy=10.0,
        mm_w_st=8.0,
        mm_r_st=0.0,
        mm_r_sf=0.0,
        mm_r_sb=0.0,
        Q=12,
        location=None,
    ):
        self.name = name
        self.color = color
        self.deg_alpha_st = deg_alpha_st
        self.deg_alpha_sto = deg_alpha_sto
        self.mm_r_so = mm_r_so
        self.mm_d_sto = mm_d_sto
        self.mm_d_stt = mm_d_stt
        self.mm_d_st = mm_d_st
        self.mm_d_sy = mm_d_sy
        self.mm_w_st = mm_w_st
        self.mm_r_st = mm_r_st
        self.mm_r_sf = mm_r_sf
        self.mm_r_sb = mm_r_sb
        self.Q = Q
        self.location = location

        self._validate_parameters()

    @property
    def mm_r_sy_inner(self):
        return self.mm_r_so - self.mm_d_stt - self.mm_d_st - self.mm_d_sy

    @property
    def mm_r_slot_opening(self):
        return self.mm_r_so - self.mm_d_stt

    def _validate_parameters(self):
        if not isinstance(self.Q, int) or self.Q < 3:
            raise ValueError("Q must be an integer greater than or equal to 3.")

        positive_parameters = {
            "mm_r_so": self.mm_r_so,
            "mm_d_sto": self.mm_d_sto,
            "mm_d_stt": self.mm_d_stt,
            "mm_d_st": self.mm_d_st,
            "mm_d_sy": self.mm_d_sy,
            "mm_w_st": self.mm_w_st,
        }
        for name, value in positive_parameters.items():
            if value <= 0:
                raise ValueError(f"{name} must be positive.")

        for name, value in {
            "mm_r_st": self.mm_r_st,
            "mm_r_sf": self.mm_r_sf,
            "mm_r_sb": self.mm_r_sb,
        }.items():
            if value < 0:
                raise ValueError(f"{name} cannot be negative.")

        slot_span_deg = 360.0 / self.Q
        if not 0 < self.deg_alpha_st < slot_span_deg:
            raise ValueError("deg_alpha_st must be smaller than one stator slot pitch.")
        if not 0 <= self.deg_alpha_sto < 90:
            raise ValueError("deg_alpha_sto must be in the range [0, 90).")
        if self.mm_r_sy_inner <= 0:
            raise ValueError(
                "The radial stack is invalid: mm_r_so must exceed "
                "mm_d_stt + mm_d_st + mm_d_sy."
            )

        half_tooth_angle = np.arctan2(0.5 * self.mm_w_st, self.mm_r_slot_opening)
        if half_tooth_angle >= np.pi / self.Q:
            raise ValueError("mm_w_st is too large for the selected stator slot pitch.")

    def _build_points(self):
        alpha_st = np.deg2rad(self.deg_alpha_st)
        alpha_sto = np.deg2rad(self.deg_alpha_sto)
        alpha_slot_span = 2.0 * np.pi / self.Q

        p1 = [self.mm_r_so, 0.0]
        p2 = _polar(self.mm_r_so, -0.5 * alpha_st)

        edge_angle = np.pi - 0.5 * alpha_st + alpha_sto
        p3 = [
            p2[0] + self.mm_d_sto * np.cos(edge_angle),
            p2[1] + self.mm_d_sto * np.sin(edge_angle),
        ]

        half_tooth_angle = np.arctan2(0.5 * self.mm_w_st, self.mm_r_slot_opening)
        p4 = _polar(self.mm_r_slot_opening, -half_tooth_angle)
        p5 = [p4[0] - self.mm_d_st, p4[1]]

        radius_yoke_outer = float(np.hypot(p5[0], p5[1]))
        tooth_body_half_angle = abs(np.arctan2(p5[1], p5[0]))
        half_slot_span = 0.5 * alpha_slot_span
        if p5[0] <= 0 or tooth_body_half_angle >= half_slot_span:
            raise ValueError(
                "The stator teeth intersect at the slot bottom. Reduce mm_w_st "
                "or mm_d_st, or increase the stator inner radius."
            )

        p6 = _polar(radius_yoke_outer, -0.5 * alpha_slot_span)
        p7 = _polar(self.mm_r_sy_inner, -0.5 * alpha_slot_span)
        p8 = [self.mm_r_sy_inner, 0.0]

        if np.hypot(p3[0], p3[1]) <= self.mm_r_slot_opening:
            raise ValueError(
                "The tooth-edge geometry intersects the tooth body. "
                "Reduce mm_d_sto or deg_alpha_sto."
            )
        if radius_yoke_outer <= self.mm_r_sy_inner:
            raise ValueError("The tooth body leaves no room for the stator yoke.")

        return {
            "P1": p1,
            "P2": p2,
            "P3": p3,
            "P4": p4,
            "P5": p5,
            "P6": p6,
            "P7": p7,
            "P8": p8,
            "alpha_slot_span": alpha_slot_span,
            "radius_yoke_outer": radius_yoke_outer,
            "slot_bottom_span": 2.0 * (half_slot_span - tooth_body_half_angle),
        }

    def draw(self, drawer, bool_draw_whole_model=False):
        drawer.getSketch(self.name, self.color)
        points = self._build_points()
        p1 = points["P1"]
        p2 = points["P2"]
        p3 = points["P3"]
        p4 = points["P4"]
        p5 = points["P5"]
        p6 = points["P6"]
        p7 = points["P7"]
        p8 = points["P8"]
        alpha_slot_span = points["alpha_slot_span"]

        list_segments = []
        if bool_draw_whole_model:
            p2_mirror = _mirror_x(p2)
            p3_mirror = _mirror_x(p3)
            p4_mirror = _mirror_x(p4)
            p5_mirror = _mirror_x(p5)

            for index in range(self.Q):
                angle = index * alpha_slot_span
                p5_next = _rotate(p5, angle + alpha_slot_span)
                list_segments += drawer.drawArc(
                    [0, 0], _rotate(p2, angle), _rotate(p2_mirror, angle)
                )
                list_segments += drawer.drawLine(
                    _rotate(p2, angle), _rotate(p3, angle)
                )
                list_segments += drawer.drawLine(
                    _rotate(p2_mirror, angle), _rotate(p3_mirror, angle)
                )
                list_segments += drawer.drawLine(
                    _rotate(p3, angle), _rotate(p4, angle)
                )
                list_segments += drawer.drawLine(
                    _rotate(p3_mirror, angle), _rotate(p4_mirror, angle)
                )
                list_segments += drawer.drawLine(
                    _rotate(p4, angle), _rotate(p5, angle)
                )
                list_segments += drawer.drawLine(
                    _rotate(p4_mirror, angle), _rotate(p5_mirror, angle)
                )
                list_segments += drawer.drawArc(
                    [0, 0], _rotate(p5_mirror, angle), p5_next
                )

            p8_negative = [-p8[0], p8[1]]
            list_segments += drawer.drawArc([0, 0], p8, p8_negative)
            list_segments += drawer.drawArc([0, 0], p8_negative, p8)
        else:
            list_segments += drawer.drawArc([0, 0], p2, p1)
            list_segments += drawer.drawLine(p2, p3)
            list_segments += drawer.drawLine(p3, p4)
            list_segments += drawer.drawLine(p4, p5)
            list_segments += drawer.drawArc([0, 0], p6, p5)
            list_segments += drawer.drawLine(p6, p7)
            list_segments += drawer.drawArc([0, 0], p7, p8)
            list_segments += drawer.drawLine(p8, p1)

        self.innerCoord = (
            0.5 * (self.mm_r_so + points["radius_yoke_outer"]),
            0.0,
        )
        return {
            "innerCoord": self.innerCoord,
            "list_regions": [list_segments],
            "mirrorAxis": [
                (self.mm_r_so + 5.0, 0.0),
                (self.mm_r_so + 15.0, 0.0),
            ],
        }


class CrossSectOuterRotorStatorWinding:
    """Winding regions associated with :class:`CrossSectOuterRotorStator`."""

    def __init__(self, name="Coils", color="#3D9970", stator_core=None):
        if stator_core is None:
            raise ValueError("stator_core is required.")
        self.type = "Cynlinder"
        self.name = name
        self.color = color
        self.stator_core = stator_core

    @staticmethod
    def _shrink(center, point, factor=0.9):
        return [
            center[0] + factor * (point[0] - center[0]),
            center[1] + factor * (point[1] - center[1]),
        ]

    def _build_slot_points(self):
        core_points = self.stator_core._build_points()
        alpha_slot_span = core_points["alpha_slot_span"]
        p_open = _polar(
            self.stator_core.mm_r_slot_opening,
            -0.5 * alpha_slot_span,
        )
        p4 = core_points["P4"]
        p5 = core_points["P5"]
        p6 = core_points["P6"]

        midpoint_45 = [0.5 * (p4[0] + p5[0]), 0.5 * (p4[1] + p5[1])]
        midpoint_6_open = [
            0.5 * (p6[0] + p_open[0]),
            0.5 * (p6[1] + p_open[1]),
        ]
        coil_center = [
            0.5 * (midpoint_45[0] + midpoint_6_open[0]),
            0.5 * (midpoint_45[1] + midpoint_6_open[1]),
        ]

        return {
            "POpen": p_open,
            "P4": p4,
            "P5": p5,
            "P6": p6,
            "PCoil": coil_center,
            "alpha_slot_span": alpha_slot_span,
        }

    def draw(self, drawer, bool_re_evaluate=False, bool_draw_whole_model=False):
        slot_points = self._build_slot_points()
        p_open = slot_points["POpen"]
        p4 = slot_points["P4"]
        p5 = slot_points["P5"]
        p6 = slot_points["P6"]
        self.PCoil = slot_points["PCoil"]
        alpha_slot_span = slot_points["alpha_slot_span"]

        self.mm2_slot_area = 2.0 * _polygon_area([p4, p5, p6, p_open])
        if bool_re_evaluate:
            return self.mm2_slot_area

        drawer.getSketch(self.name, self.color)
        p6_shrink = self._shrink(self.PCoil, p6)
        p5_shrink = self._shrink(self.PCoil, p5)
        p4_shrink = self._shrink(self.PCoil, p4)
        p_open_shrink = self._shrink(self.PCoil, p_open)

        list_regions = []
        if bool_draw_whole_model:
            lower_points = [p6_shrink, p5_shrink, p4_shrink, p_open_shrink]
            upper_points = [_mirror_x(point) for point in lower_points]
            for index in range(self.stator_core.Q):
                angle = index * alpha_slot_span

                lower = [_rotate(point, angle) for point in lower_points]
                lower_segments = []
                lower_segments += drawer.drawArc([0, 0], lower[0], lower[1])
                lower_segments += drawer.drawLine(lower[1], lower[2])
                lower_segments += drawer.drawLine(lower[2], lower[3])
                lower_segments += drawer.drawLine(lower[3], lower[0])
                list_regions.append(lower_segments)

                upper = [_rotate(point, angle) for point in upper_points]
                upper_segments = []
                upper_segments += drawer.drawArc([0, 0], upper[1], upper[0])
                upper_segments += drawer.drawLine(upper[1], upper[2])
                upper_segments += drawer.drawLine(upper[2], upper[3])
                upper_segments += drawer.drawLine(upper[3], upper[0])
                list_regions.append(upper_segments)
        else:
            lower_segments = []
            lower_segments += drawer.drawArc([0, 0], p6_shrink, p5_shrink)
            lower_segments += drawer.drawLine(p5_shrink, p4_shrink)
            lower_segments += drawer.drawLine(p4_shrink, p_open_shrink)
            lower_segments += drawer.drawLine(p_open_shrink, p6_shrink)
            list_regions.append(lower_segments)

            p6_upper = _mirror_x(p6_shrink)
            p5_upper = _mirror_x(p5_shrink)
            p4_upper = _mirror_x(p4_shrink)
            p_open_upper = _mirror_x(p_open_shrink)
            upper_segments = []
            upper_segments += drawer.drawArc([0, 0], p5_upper, p6_upper)
            upper_segments += drawer.drawLine(p5_upper, p4_upper)
            upper_segments += drawer.drawLine(p4_upper, p_open_upper)
            upper_segments += drawer.drawLine(p_open_upper, p6_upper)
            list_regions.append(upper_segments)

        self.innerCoord = (
            0.5 * (p_open[0] + p6[0]),
            0.5 * (p_open[1] + p6[1]),
        )
        return {
            "innerCoord": self.innerCoord,
            "list_regions": list_regions,
            "mirrorAxis": None,
        }


# Topology-oriented and geometry-oriented names are both supported.
CrossSectInnerStator = CrossSectOuterRotorStator
CrossSectInnerStatorWinding = CrossSectOuterRotorStatorWinding
