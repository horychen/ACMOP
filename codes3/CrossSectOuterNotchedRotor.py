import logging

import numpy as np


EPS = 1e-3


class ExceptionBadDesign(Exception):
    """Raised when an outer-rotor geometry cannot form valid regions."""


def _polar(radius, angle):
    return [radius * np.cos(angle), radius * np.sin(angle)]


def _rotate(point, angle):
    cos_angle = np.cos(angle)
    sin_angle = np.sin(angle)
    return [
        point[0] * cos_angle - point[1] * sin_angle,
        point[0] * sin_angle + point[1] * cos_angle,
    ]


def _distance(point_a, point_b):
    return float(np.hypot(point_a[0] - point_b[0], point_a[1] - point_b[1]))


def _append_line(drawer, segments, start, end):
    if _distance(start, end) > EPS:
        segments += drawer.drawLine(start, end)


def _append_arc(drawer, segments, start, end):
    start_angle = np.arctan2(start[1], start[0])
    end_angle = np.arctan2(end[1], end[0])
    angle_difference = abs((end_angle - start_angle + np.pi) % (2 * np.pi) - np.pi)
    if angle_difference > EPS:
        segments += drawer.drawArc([0, 0], start, end)


class CrossSectOuterNotchedRotor:
    """Notched outer rotor whose magnets face radially inward.

    ``mm_r_ro`` is the rotor outer radius. The steel yoke extends inward by
    ``mm_d_ri``. Magnets are placed on the yoke inner surface and extend
    inward by ``mm_d_pm`` toward the stator.
    """

    def __init__(
        self,
        name="OuterNotchedRotor",
        color="#FE840E",
        mm_d_pm=4.0,
        deg_alpha_rm=30.0,
        deg_alpha_rs=30.0,
        mm_d_ri=8.0,
        mm_r_ro=74.0,
        mm_d_rp=4.0,
        mm_d_rs=0.0,
        p=5,
        s=1,
        location=None,
    ):
        self.name = name
        self.color = color
        self.mm_d_pm = mm_d_pm
        self.deg_alpha_rm = deg_alpha_rm
        self.deg_alpha_rs = deg_alpha_rs
        self.mm_d_ri = mm_d_ri
        self.mm_r_ro = mm_r_ro
        self.mm_d_rp = mm_d_rp
        self.mm_d_rs = mm_d_rs
        self.p = p
        self.s = s
        self.location = location

        self._validate_parameters()

    @property
    def mm_r_ry_inner(self):
        """Inner radius of the continuous outer-rotor steel yoke."""
        return self.mm_r_ro - self.mm_d_ri

    @property
    def mm_r_magnet_inner(self):
        """Air-gap-facing radius of the magnets."""
        return self.mm_r_ry_inner - self.mm_d_pm

    @property
    def alpha_rp(self):
        return np.pi / self.p

    def _validate_parameters(self):
        if not isinstance(self.p, int) or self.p < 1:
            raise ValueError("p must be a positive integer.")
        if not isinstance(self.s, int) or self.s < 1:
            raise ValueError("s must be a positive integer.")

        for name, value in {
            "mm_r_ro": self.mm_r_ro,
            "mm_d_ri": self.mm_d_ri,
            "mm_d_pm": self.mm_d_pm,
        }.items():
            if value <= 0:
                raise ValueError(f"{name} must be positive.")
        for name, value in {
            "mm_d_rp": self.mm_d_rp,
            "mm_d_rs": self.mm_d_rs,
        }.items():
            if value < 0:
                raise ValueError(f"{name} cannot be negative.")

        if self.mm_r_magnet_inner <= 0:
            raise ValueError("The rotor yoke and magnet depths exceed mm_r_ro.")
        if self.mm_d_rp > self.mm_d_pm:
            raise ValueError("mm_d_rp cannot exceed mm_d_pm.")
        if self.mm_d_rs > self.mm_d_pm:
            raise ValueError("mm_d_rs cannot exceed mm_d_pm.")

        pole_pitch_deg = 180.0 / self.p
        if not 0 < self.deg_alpha_rm <= pole_pitch_deg:
            raise ValueError("deg_alpha_rm must be in the range (0, 180/p].")
        if self.s == 1:
            if not np.isclose(self.deg_alpha_rs, self.deg_alpha_rm):
                raise ValueError("For s=1, deg_alpha_rs must equal deg_alpha_rm.")
            if not np.isclose(self.mm_d_rs, 0.0):
                raise ValueError("For s=1, mm_d_rs must be zero.")
        else:
            if self.mm_d_rs <= 0:
                raise ValueError("For s>1, mm_d_rs must be positive.")
            if self.deg_alpha_rs * self.s > self.deg_alpha_rm + 1e-9:
                raise ValueError("The magnet segments exceed deg_alpha_rm.")

    def _effective_angles(self):
        alpha_rm = np.deg2rad(self.deg_alpha_rm)
        alpha_rs = np.deg2rad(self.deg_alpha_rs)
        if abs(self.alpha_rp - alpha_rm) <= np.deg2rad(2.0):
            alpha_rm = self.alpha_rp
            if self.s == 1:
                alpha_rs = alpha_rm
        if self.s == 1:
            alpha_notch = 0.0
        else:
            alpha_notch = max(0.0, (alpha_rm - self.s * alpha_rs) / (self.s - 1))
        return alpha_rm, alpha_rs, alpha_notch

    def _draw_one_pole(self, drawer, angle_offset=0.0, close_to_outer_radius=True):
        alpha_rm, alpha_rs, alpha_notch = self._effective_angles()
        alpha_start = self.alpha_rp - alpha_rm
        radius_yoke_inner = self.mm_r_ry_inner
        radius_interpole = radius_yoke_inner - self.mm_d_rp
        radius_intersegment = radius_yoke_inner - self.mm_d_rs

        def point(radius, local_angle):
            return _polar(radius, angle_offset - local_angle)

        segments = []
        p2 = point(radius_interpole, 0.0)
        p3 = point(radius_interpole, alpha_start)
        p4 = point(radius_yoke_inner, alpha_start)

        if close_to_outer_radius:
            p1 = point(self.mm_r_ro, 0.0)
            if alpha_start <= EPS:
                # At full pole arc p2 and p3 coincide. Drawing p1->p2->p4
                # would retrace the radial edge and JMAG drops the region.
                _append_line(drawer, segments, p1, p4)
            else:
                _append_line(drawer, segments, p1, p2)
        if alpha_start > EPS:
            _append_arc(drawer, segments, p3, p2)
            _append_line(drawer, segments, p3, p4)

        segment_start = alpha_start
        current_yoke_point = p4
        for segment_index in range(self.s):
            segment_end = segment_start + alpha_rs
            segment_end_yoke = point(radius_yoke_inner, segment_end)
            _append_arc(drawer, segments, segment_end_yoke, current_yoke_point)

            if segment_index < self.s - 1:
                notch_start = point(radius_intersegment, segment_end)
                notch_end_angle = segment_end + alpha_notch
                notch_end = point(radius_intersegment, notch_end_angle)
                next_yoke_point = point(radius_yoke_inner, notch_end_angle)
                _append_line(drawer, segments, segment_end_yoke, notch_start)
                _append_arc(drawer, segments, notch_end, notch_start)
                _append_line(drawer, segments, notch_end, next_yoke_point)
                current_yoke_point = next_yoke_point
                segment_start = notch_end_angle
            else:
                current_yoke_point = segment_end_yoke

        if close_to_outer_radius:
            p_outer_end = point(self.mm_r_ro, self.alpha_rp)
            _append_line(drawer, segments, current_yoke_point, p_outer_end)
            _append_arc(drawer, segments, p_outer_end, p1)
        else:
            next_interpole_point = point(radius_interpole, self.alpha_rp)
            _append_line(
                drawer,
                segments,
                current_yoke_point,
                next_interpole_point,
            )
        return segments

    def draw(self, drawer, bool_draw_whole_model=False):
        drawer.getSketch(self.name, self.color)

        if bool_draw_whole_model:
            list_segments = []
            for pole_index in range(2 * self.p):
                list_segments.extend(
                    self._draw_one_pole(
                        drawer,
                        pole_index * self.alpha_rp,
                        close_to_outer_radius=False,
                    )
                )
            outer_positive = [self.mm_r_ro, 0.0]
            outer_negative = [-self.mm_r_ro, 0.0]
            list_segments += drawer.drawArc(
                [0, 0], outer_positive, outer_negative
            )
            list_segments += drawer.drawArc(
                [0, 0], outer_negative, outer_positive
            )
        else:
            list_segments = self._draw_one_pole(drawer)

        self.innerCoord = (
            0.5 * (self.mm_r_ro + self.mm_r_ry_inner),
            0.0,
        )
        return {
            "innerCoord": self.innerCoord,
            "list_regions": [list_segments],
            "mirrorAxis": None,
        }


class CrossSectOuterNotchedMagnet:
    """Surface magnets attached to the inner face of an outer rotor."""

    def __init__(
        self,
        name="OuterRotorMagnet",
        color="#0BA0E2",
        notched_rotor=None,
    ):
        if notched_rotor is None:
            raise ValueError("notched_rotor is required.")
        self.name = name
        self.color = color
        self.notched_rotor = notched_rotor

    def _segment_regions(self, drawer, pole_offset=0.0):
        rotor = self.notched_rotor
        alpha_rm, alpha_rs, alpha_notch = rotor._effective_angles()
        segment_start = rotor.alpha_rp - alpha_rm
        radius_outer = rotor.mm_r_ry_inner
        radius_inner = rotor.mm_r_magnet_inner
        regions = []

        def point(radius, local_angle):
            return _polar(radius, pole_offset - local_angle)

        for segment_index in range(rotor.s):
            segment_end = segment_start + alpha_rs
            outer_start = point(radius_outer, segment_start)
            inner_start = point(radius_inner, segment_start)
            inner_end = point(radius_inner, segment_end)
            outer_end = point(radius_outer, segment_end)

            segments = []
            _append_line(drawer, segments, inner_start, outer_start)
            _append_arc(drawer, segments, outer_end, outer_start)
            _append_line(drawer, segments, outer_end, inner_end)
            _append_arc(drawer, segments, inner_end, inner_start)
            regions.append(segments)

            segment_start = segment_end + alpha_notch
        return regions

    def draw(self, drawer, bool_re_evaluate=False, bool_draw_whole_model=False):
        rotor = self.notched_rotor
        alpha_rm, _, _ = rotor._effective_angles()
        radius_outer = rotor.mm_r_ry_inner
        radius_inner = rotor.mm_r_magnet_inner
        self.mm2_magnet_area = (
            alpha_rm / rotor.alpha_rp
            * np.pi
            * (radius_outer**2 - radius_inner**2)
        )
        logging.getLogger(__name__).info(
            "Outer-rotor magnet area in total is %g mm^2",
            self.mm2_magnet_area,
        )
        if bool_re_evaluate:
            return self.mm2_magnet_area

        drawer.getSketch(self.name, self.color)
        list_regions = []
        if bool_draw_whole_model:
            for pole_index in range(2 * rotor.p):
                list_regions.extend(
                    self._segment_regions(drawer, pole_index * rotor.alpha_rp)
                )
        else:
            list_regions = self._segment_regions(drawer)

        first_region_angle = rotor.alpha_rp - alpha_rm + 0.5 * np.deg2rad(
            rotor.deg_alpha_rs
        )
        radius_middle = 0.5 * (radius_outer + radius_inner)
        self.innerCoord = tuple(_polar(radius_middle, -first_region_angle))
        return {
            "innerCoord": self.innerCoord,
            "list_regions": list_regions,
            "mirrorAxis": None,
        }


class CrossSectOuterRotorSleeve:
    """Retaining sleeve located on the air-gap side of the outer-rotor magnets."""

    def __init__(
        self,
        name="OuterRotorSleeve",
        color="#11E322",
        notched_magnet=None,
        d_sleeve=1.0,
    ):
        if notched_magnet is None:
            raise ValueError("notched_magnet is required.")
        if d_sleeve <= 0:
            raise ValueError("d_sleeve must be positive.")
        self.name = name
        self.color = color
        self.notched_magnet = notched_magnet
        self.d_sleeve = d_sleeve

        if self.mm_r_si <= 0:
            raise ValueError("The sleeve depth leaves no positive inner radius.")

    @property
    def mm_r_so(self):
        return self.notched_magnet.notched_rotor.mm_r_magnet_inner

    @property
    def mm_r_si(self):
        return self.mm_r_so - self.d_sleeve

    def _draw_sector(self, drawer, angle_offset=0.0):
        pole_span = self.notched_magnet.notched_rotor.alpha_rp
        p1 = _polar(self.mm_r_si, angle_offset)
        p2 = _polar(self.mm_r_so, angle_offset)
        p3 = _polar(self.mm_r_so, angle_offset - pole_span)
        p4 = _polar(self.mm_r_si, angle_offset - pole_span)

        segments = []
        _append_line(drawer, segments, p1, p2)
        _append_arc(drawer, segments, p3, p2)
        _append_line(drawer, segments, p3, p4)
        _append_arc(drawer, segments, p4, p1)
        return segments

    def draw(self, drawer, bool_draw_whole_model=False):
        drawer.getSketch(self.name, self.color)
        rotor = self.notched_magnet.notched_rotor

        if bool_draw_whole_model:
            inner_positive = [self.mm_r_si, 0.0]
            inner_negative = [-self.mm_r_si, 0.0]
            outer_positive = [self.mm_r_so, 0.0]
            outer_negative = [-self.mm_r_so, 0.0]
            list_segments = []
            list_segments += drawer.drawArc(
                [0, 0], inner_positive, inner_negative
            )
            list_segments += drawer.drawArc(
                [0, 0], inner_negative, inner_positive
            )
            list_segments += drawer.drawArc(
                [0, 0], outer_positive, outer_negative
            )
            list_segments += drawer.drawArc(
                [0, 0], outer_negative, outer_positive
            )
        else:
            list_segments = self._draw_sector(drawer)

        self.innerCoord = (0.5 * (self.mm_r_so + self.mm_r_si), 0.0)
        return {
            "innerCoord": self.innerCoord,
            "list_regions": [list_segments],
            "mirrorAxis": None,
        }


# Shorter aliases for template code.
CrossSectOuterRotor = CrossSectOuterNotchedRotor
CrossSectOuterMagnet = CrossSectOuterNotchedMagnet
CrossSectSleeve = CrossSectOuterRotorSleeve
