"use client";

import React, { useEffect, useState, useMemo, useRef, useCallback } from "react";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/alert";
import { Loader2, AlertTriangle, RefreshCw } from "lucide-react";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from "@/components/ui/select";

/** Sections and param keys aligned with backend/codes4/user_minitureMachine.py (GeometrySpecs, WindingSpecs, MaterialSpecs, PerformanceTargets) */
const PARAM_SECTIONS: {
  section: "geometry" | "winding" | "materials" | "targets";
  title: string;
  keys: { key: string; label: string; type?: "number" | "text" | "select"; optionsKey?: string }[];
}[] = [
    {
      section: "geometry",
      title: "Geometry",
      keys: [
        { key: "r_stator_outer", label: "Stator OR (mm)" },
        { key: "r_rotor_outer", label: "Rotor OR (mm)" },
        { key: "r_shaft", label: "Shaft R (mm)" },
        { key: "d_air_gap", label: "Air gap (mm)" },
        { key: "d_magnet", label: "Magnet thick. (mm)" },
        { key: "w_tooth", label: "Tooth width (mm)" },
        { key: "d_tooth", label: "Tooth depth (mm)" },
        { key: "d_stator_yoke", label: "Yoke depth (mm)" },
        { key: "d_tooth_shoe", label: "Tooth shoe (mm)" },
        { key: "tooth_shape", label: "Tooth Shape", type: "select", optionsKey: "tooth_shape_options" },
        { key: "l_stack", label: "Stack (mm)", type: "select", optionsKey: "stack_length_options" },
        { key: "split_ratio", label: "Split Ratio" },
      ],
    },
    {
      section: "winding",
      title: "Winding",
      keys: [
        { key: "num_slots", label: "Slots" },
        { key: "num_poles", label: "Poles" },
        { key: "coil_pitch_y", label: "Coil pitch y" },
        { key: "conductors_per_slot", label: "Conductors/slot" },
        { key: "wire_diameter_with_insulation", label: "Wire diam. (mm)" },
        { key: "rated_current_density", label: "J_s (A/mm²)" },
        { key: "rated_speed", label: "Rated Speed (rpm)" },
        { key: "dc_bus_voltage", label: "DC Bus (V)" },
      ],
    },
    {
      section: "materials",
      title: "Materials",
      keys: [
        { key: "magnet_grade", label: "Magnet grade", type: "text" },
        { key: "magnet_br", label: "B_r (T)" },
        { key: "magnet_h_cj", label: "H_cj (kA/m)" },
        { key: "magnet_temp_max", label: "Magnet T_max (°C)" },
        { key: "magnet_temperature", label: "Operation Temp (°C)" },
        { key: "stator_steel", label: "Stator steel (legacy)", type: "text" },
        { key: "stator_core_material", label: "Stator core", type: "text" },
        { key: "rotor_core_material", label: "Rotor core", type: "text" },
        { key: "steel_thickness", label: "Steel thick. (mm)" },
        { key: "steel_stack_factor", label: "Stack factor" },
        { key: "lamination_factor", label: "Lamination (%)" },
        { key: "steel_max_flux_density", label: "B_sat (T)" },
      ],
    },
    {
      section: "targets",
      title: "Performance",
      keys: [
        { key: "torque_target", label: "Torque Goal (mNm)" },
        { key: "efficiency_target", label: "Eff. Goal" },
        { key: "torque_ripple", label: "Ripple (%)", type: "text" },
      ],
    },
  ];

const PARAM_KEYS = PARAM_SECTIONS.flatMap(({ section, keys }) =>
  keys.map((k) => ({ section, key: k.key, label: k.label, type: k.type ?? "number" }))
);

function paramKey(section: string, key: string) {
  return `${section}.${key}`;
}

interface ParameterObj {
  name: string;
  type: "free" | "fixed" | "derived";
  value: number | string;
  bounds?: number[];
  unit?: string;
}

function specsToInputValues(specs: MachineSpecs): Record<string, string> {
  const out: Record<string, string> = {};
  for (const { section, key } of PARAM_KEYS) {
    const sec = specs[section as keyof MachineSpecs];
    if (!sec || typeof sec !== "object") continue;
    const v = (sec as Record<string, unknown>)[key];
    if (v !== undefined && v !== null) {
      if (typeof v === "object" && v !== null && "value" in v) {
        out[paramKey(section, key)] = String((v as any).value);
      } else {
        out[paramKey(section, key)] = String(v);
      }
    }
  }
  return out;
}

/** Matches backend MotorSpecs (GeometrySpecs, WindingSpecs, MaterialSpecs, PerformanceTargets) */
interface MachineSpecs {
  geometry: {
    r_stator_outer: ParameterObj;
    r_rotor_outer: ParameterObj;
    r_shaft: ParameterObj;
    d_air_gap: ParameterObj;
    d_magnet: ParameterObj;
    w_tooth: ParameterObj;
    d_tooth: ParameterObj;
    d_stator_yoke: ParameterObj;
    d_tooth_shoe: ParameterObj;
    tooth_shape: string;
    tooth_shape_options?: string[];
    l_stack: ParameterObj;
    stack_length_options?: number[];
    split_ratio: ParameterObj;
  };
  winding: {
    num_slots: number;
    num_poles: number;
    coil_pitch_y: number;
    conductors_per_slot: number;
    wire_diameter_with_insulation: number;
    rated_speed: number;
    dc_bus_voltage: number;
    rated_current_density: number;
    winding_factor: number;
  };
  materials: {
    magnet_grade: string;
    magnet_br: number;
    magnet_h_cj: number;
    magnet_temp_max: number;
    magnet_temperature: number;
    stator_steel: string;
    stator_core_material: string;
    rotor_core_material: string;
    steel_thickness: number;
    steel_stack_factor: number;
    lamination_factor: number;
    steel_max_flux_density: number;
  };
  targets: {
    torque_target: number;
    efficiency_target: number;
    torque_ripple: number;
    thd_voltage: number;
  };
  validations?: Record<string, {
    status: "pass" | "warn" | "fail";
    errors: string[];
    warnings: string[];
    metrics?: Record<string, string>;
    fill_factor?: number | string;
  }>;
}

const MachineCrossSection = React.memo(({
  geometry,
  winding,
  size = 400,
  showLabels = true
}: {
  geometry: MachineSpecs["geometry"],
  winding: MachineSpecs["winding"],
  size?: number,
  showLabels?: boolean
}) => {
  const R_so = (typeof geometry.r_stator_outer.value === "number" ? geometry.r_stator_outer.value : 13.0 / 2);
  const R_ro = (typeof geometry.r_rotor_outer.value === "number" ? geometry.r_rotor_outer.value : 8.0 / 2);
  const gap = (typeof geometry.d_air_gap.value === "number" ? geometry.d_air_gap.value : 0.15);
  const R_si = R_so * (typeof geometry.split_ratio.value === "number" ? geometry.split_ratio.value : 0.615);
  const toothDepth = (typeof geometry.d_tooth.value === "number" ? geometry.d_tooth.value : 2.0);
  const R_sb = R_si + toothDepth; // Slot bottom
  const toothWidth = (typeof geometry.w_tooth.value === "number" ? geometry.w_tooth.value : 1.2);
  const shoeDepth = (typeof geometry.d_tooth_shoe.value === "number" ? geometry.d_tooth_shoe.value : 0.5);
  const toothShape = geometry.tooth_shape;

  const numSlots = winding.num_slots;
  const numPoles = winding.num_poles;
  const conductorsPerSlot = winding.conductors_per_slot;
  const wireDiam = winding.wire_diameter_with_insulation;
  const wireRad = wireDiam / 2;
  const magnetThickness = (typeof geometry.d_magnet.value === "number" ? geometry.d_magnet.value : 3.0);
  const R_shaft = (typeof geometry.r_shaft.value === "number" ? geometry.r_shaft.value : 2.0) / 2;

  const getPt = (r: number, theta: number) => {
    const rad = (theta * Math.PI) / 180;
    return { x: r * Math.cos(rad), y: r * Math.sin(rad) };
  };

  const R_mag_in = R_ro - magnetThickness;
  const hasBackIron = R_mag_in > R_shaft + 0.01;

  const magnetPaths: { d: string, color: string }[] = [];
  if (numPoles > 0) {
    const poleAngle = 360 / numPoles;
    for (let i = 0; i < numPoles; i++) {
      const startAngle = i * poleAngle;
      const endAngle = (i + 1) * poleAngle;
      const p1 = getPt(R_ro, startAngle);
      const p2 = getPt(R_ro, endAngle);
      const p3 = getPt(R_mag_in, endAngle);
      const p4 = getPt(R_mag_in, startAngle);
      let d = `M ${p1.x} ${p1.y} A ${R_ro} ${R_ro} 0 0 1 ${p2.x} ${p2.y} L ${p3.x} ${p3.y} A ${R_mag_in} ${R_mag_in} 0 0 0 ${p4.x} ${p4.y} Z`;
      magnetPaths.push({ d, color: i % 2 === 0 ? "#fca5a5" : "#93c5fd" });
    }
  }

  // Radial Partition Lines (Slot Center Axis)
  const partitionLines: { x1: number, y1: number, x2: number, y2: number }[] = [];
  for (let i = 0; i < numSlots; i++) {
    const angle = (i + 0.5) * 360 / numSlots;
    const p1 = getPt(R_si, angle);
    const p2 = getPt(R_sb, angle);
    partitionLines.push({ x1: p1.x, y1: p1.y, x2: p2.x, y2: p2.y });
  }

  const getToothPoints = (angle: number, r: number, w: number) => {
    const rad = (angle * Math.PI) / 180;
    const u = { x: Math.cos(rad), y: Math.sin(rad) };
    const n = { x: -Math.sin(rad), y: Math.cos(rad) };
    const h = w / 2;
    if (r < h) return null;
    const d = Math.sqrt(r * r - h * h);
    return { cw: { x: d * u.x - h * n.x, y: d * u.y - h * n.y }, ccw: { x: d * u.x + h * n.x, y: d * u.y + h * n.y } };
  };

  let d_path = "";
  for (let i = 0; i < numSlots; i++) {
    const angle = i * 360 / numSlots;
    const root = getToothPoints(angle, R_sb, toothWidth);
    const stalk = getToothPoints(angle, R_si + shoeDepth, toothWidth);
    let curW = toothWidth;
    if (toothShape === "semi-closed") curW = (2 * Math.PI * R_si / numSlots) - (toothWidth * 0.4);
    else if (toothShape === "closed") curW = (2 * Math.PI * R_si / numSlots) - 0.05;
    const tip = getToothPoints(angle, R_si, Math.min(curW, 2 * Math.PI * R_si / numSlots - 0.1));
    if (!tip || !root || !stalk) continue;

    const nextAngle = (i + 1) * 360 / numSlots;
    const nextRoot = getToothPoints(nextAngle, R_sb, toothWidth);
    const nextStalk = getToothPoints(nextAngle, R_si + shoeDepth, toothWidth);
    let ncurW = toothWidth;
    if (toothShape === "semi-closed") ncurW = (2 * Math.PI * R_si / numSlots) - (toothWidth * 0.4);
    else if (toothShape === "closed") ncurW = (2 * Math.PI * R_si / numSlots) - 0.05;
    const nextTip = getToothPoints(nextAngle, R_si, Math.min(ncurW, 2 * Math.PI * R_si / numSlots - 0.1));
    if (!nextRoot || !nextStalk || !nextTip) continue;

    if (i === 0) d_path += `M ${tip.cw.x} ${tip.cw.y} `;
    d_path += `A ${R_si} ${R_si} 0 0 1 ${tip.ccw.x} ${tip.ccw.y} L ${stalk.ccw.x} ${stalk.ccw.y} L ${root.ccw.x} ${root.ccw.y} A ${R_sb} ${R_sb} 0 0 1 ${nextRoot.cw.x} ${nextRoot.cw.y} L ${nextStalk.cw.x} ${nextStalk.cw.y} L ${nextTip.cw.x} ${nextTip.cw.y} `;
  }
  d_path += "Z";
  const fullStatorPath = `M ${R_so} 0 A ${R_so} ${R_so} 0 1 1 ${-R_so} 0 A ${R_so} ${R_so} 0 1 1 ${R_so} 0 Z ` + d_path;

  // Concentrated Winding Arrangement (Packing against tooth walls)
  const conductorsPerSide = Math.floor(conductorsPerSlot / 2);
  const effR = wireRad + 0.005;
  const dy = effR * 2;
  const dx = effR * Math.sqrt(3);

  const points: { x: number, y: number }[] = [];
  const rMin = R_si + shoeDepth + wireRad + 0.05;
  const rMax = R_sb - wireRad - 0.05;
  const slotCenterAngle = (180 / numSlots); // Half slot angle in degrees

  // For one tooth side (the CCW side of the tooth at angle 0)
  // This side is in the slot at positive angles.
  for (let layer = 0; points.length < conductorsPerSide; layer++) {
    const xDist = (toothWidth / 2) + wireRad + layer * dx;
    // Check if this layer exceeds the slot center line at any radius
    // xDist / r < sin(slotCenterAngle) => r > xDist / sin(slotCenterAngle)
    const sinCenter = Math.sin(slotCenterAngle * Math.PI / 180);
    const rLimit = xDist / sinCenter;

    if (rLimit > rMax) break; // Entire layer is beyond slot center

    const endR = Math.max(rMin, rLimit);
    const startR = rMax - (layer % 2) * wireRad;
    for (let r = startR; r >= endR && points.length < conductorsPerSide; r -= dy) {
      const theta = Math.asin(xDist / r) * 180 / Math.PI;
      points.push(getPt(r, theta));
    }
    if (layer > 20) break; // Safety
  }

  const pad = 1.1;
  const viewBoxSize = R_so * 2 * pad;
  return (
    <svg viewBox={`${-viewBoxSize / 2} ${-viewBoxSize / 2} ${viewBoxSize} ${viewBoxSize}`} className="w-full h-full">
      {hasBackIron && <circle cx="0" cy="0" r={R_ro - magnetThickness} fill="#e5e7eb" stroke="#9ca3af" strokeWidth="0.01" />}
      <circle cx="0" cy="0" r={R_shaft} fill="#fff" stroke="#9ca3af" strokeWidth="0.01" />
      {magnetPaths.map((m, i) => <path key={`m-${i}`} d={m.d} fill={m.color} stroke="#fff" strokeWidth="0.01" />)}
      <path d={fullStatorPath} fill="#4b5563" stroke="#1f2937" strokeWidth="0.01" fillRule="evenodd" />

      {/* Slot Partition Lines (Radial center axis) */}
      {partitionLines.map((l, i) => (
        <line key={`pl-${i}`} x1={l.x1} y1={l.y1} x2={l.x2} y2={l.y2} stroke="#9ca3af" strokeWidth="0.01" strokeDasharray="0.05,0.05" />
      ))}

      {numSlots > 0 && Array.from({ length: numSlots }).map((_, s) => (
        <g key={`s-${s}`} transform={`rotate(${s * 360 / numSlots})`}>
          {/* Conductors belonging to tooth "s" CCW side (in slot "s") */}
          {points.map((p, j) => (
            <circle key={`ccw-${j}`} cx={p.x} cy={p.y} r={wireRad} fill="#fbbf24" stroke="#d97706" strokeWidth="0.005" />
          ))}
          {/* Conductors belonging to tooth "s+1" CW side (in slot "s") */}
          {/* Slot center is (180/numSlots). Tooth center is at 0 and (360/numSlots). */}
          {/* Tooth "s+1" CW side mirrored across slot center line at (180/numSlots) */}
          <g transform={`rotate(${360 / numSlots}) scale(1,-1)`}>
            {points.map((p, j) => (
              <circle key={`cw-${j}`} cx={p.x} cy={p.y} r={wireRad} fill="#f59e0b" stroke="#b45309" strokeWidth="0.005" />
            ))}
          </g>
        </g>
      ))}
      {showLabels && (
        <g><line x1={R_so} y1="0" x2={R_so + 1} y2="0" stroke="#9ca3af" strokeWidth="0.01" /><text x={R_so + 1.2} y="0.5" style={{ fontSize: '1.2px' }} className="fill-muted-foreground">OD: {(R_so * 2).toFixed(1)}mm</text></g>
      )}
    </svg>
  );
});
MachineCrossSection.displayName = "MachineCrossSection";

const COUNTDOWN_SEC = 3;

interface MachineVisualizerProps {
  currentStepId?: string;
  onValidationChange?: (isValid: boolean) => void;
}

export default function MachineVisualizer({ currentStepId = "geometry", onValidationChange }: MachineVisualizerProps) {
  const [initialSpecs, setInitialSpecs] = useState<MachineSpecs | null>(null);
  const [specs, setSpecs] = useState<MachineSpecs | null>(null);
  const [inputValues, setInputValues] = useState<Record<string, string>>({});
  const [pendingUpdate, setPendingUpdate] = useState<{ section: "geometry" | "winding" | "materials" | "targets"; key: string; value: string } | null>(null);
  const [countdownSec, setCountdownSec] = useState(0);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [conductorError, setConductorError] = useState<string | null>(null);
  const [apiDebug, setApiDebug] = useState<{ tooth_width_from_py?: number; path_used?: string } | null>(null);
  const countdownRef = useRef<ReturnType<typeof setInterval> | null>(null);

  const apiBase = typeof window !== "undefined" ? "" : process.env.BACKEND_URL || "http://localhost:8000";
  const specsUrl = `${apiBase}/api/machine-specs?t=${Date.now()}`;

  const fetchSpecs = () => {
    setLoading(true);
    fetch(specsUrl)
      .then((res) => {
        if (!res.ok) throw new Error("Failed to fetch machine specs");
        return res.json();
      })
      .then((data) => {
        const normalized = {
          geometry: { ...data?.geometry },
          winding: { ...data?.winding },
          materials: { ...data?.materials },
          targets: { ...data?.targets },
          validations: data?.validations,
        } as MachineSpecs;
        const inputV = specsToInputValues(normalized);
        setApiDebug(data?._debug ?? null);
        setInitialSpecs(normalized);
        setSpecs(normalized);
        setInputValues(inputV);
        setPendingUpdate(null);
        setCountdownSec(0);
        if (countdownRef.current) {
          clearInterval(countdownRef.current);
          countdownRef.current = null;
        }
        setLoading(false);
      })
      .catch((err) => {
        setError(err.message);
        setLoading(false);
      });
  };

  useEffect(() => {
    fetchSpecs();
  }, []);

  const applyPending = useCallback(() => {
    if (!pendingUpdate || !specs) return;
    const { section, key, value } = pendingUpdate;
    const sec = specs[section as keyof MachineSpecs] as Record<string, unknown>;
    if (!sec) return;
    const current = sec[key];
    const isText = section === "materials" && (key === "magnet_grade" || key === "stator_steel");
    const numValue = parseFloat(value);
    const newVal = isText ? value : (Number.isNaN(numValue) ? current : numValue);
    setSpecs({
      ...specs,
      [section]: { ...sec, [key]: newVal }
    });
    setInputValues((prev) => ({ ...prev, [paramKey(section, key)]: value.trim() === "" ? String(sec[key] ?? "") : value }));
    setPendingUpdate(null);
    setCountdownSec(0);
    if (countdownRef.current) {
      clearInterval(countdownRef.current);
      countdownRef.current = null;
    }
  }, [pendingUpdate, specs]);

  useEffect(() => {
    if (pendingUpdate != null && countdownSec === 0) {
      applyPending();
      return;
    }
  }, [pendingUpdate, countdownSec, applyPending]);

  useEffect(() => {
    if (pendingUpdate == null || countdownSec <= 0) return;
    const id = setInterval(() => {
      setCountdownSec((prev) => (prev <= 1 ? 0 : prev - 1));
    }, 1000);
    return () => clearInterval(id);
  }, [pendingUpdate, countdownSec]);

  const handleInputChange = useCallback((section: "geometry" | "winding" | "materials" | "targets", key: string, value: string) => {
    setInputValues((prev) => ({ ...prev, [paramKey(section, key)]: value }));
    setPendingUpdate({ section, key, value });
    setCountdownSec(COUNTDOWN_SEC);
  }, []);

  useEffect(() => {
    if (specs && onValidationChange) {
      const currentValidation = specs.validations?.[currentStepId ?? "geometry"];
      onValidationChange(currentValidation?.status !== "fail");
    }
  }, [specs, currentStepId, onValidationChange]);

  useEffect(() => {
    if (specs && Object.keys(inputValues).length === 0) {
      setInputValues(specsToInputValues(specs));
    }
  }, [specs]);

  const apiStatus =
    error
      ? `请求失败: ${error}`
      : apiDebug
        ? `API Mode`
        : specs
          ? `OK`
          : loading
            ? "..."
            : "No data";

  if (loading) return <div className="flex flex-col items-center gap-4 p-8"><Loader2 className="animate-spin" /></div>;
  if (error) return (
    <div className="space-y-2 p-4">
      <Alert variant="destructive"><AlertTitle>Error</AlertTitle><AlertDescription>{error}</AlertDescription></Alert>
    </div>
  );
  if (!specs) return <div className="p-4">No data.</div>;

  const conductorsPerSlot = specs.winding.conductors_per_slot;

  const renderStepContent = () => {
    const section = currentStepId as "geometry" | "winding" | "materials" | "targets";
    const sectionConfig = PARAM_SECTIONS.find(s => s.section === section);

    if (section === "geometry") {
      return (
        <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
          <Card className="md:col-span-1">
            <CardHeader>
              <CardTitle className="flex items-center justify-between text-base">
                Geometry Parameters
                <Button variant="outline" size="icon" className="h-8 w-8" onClick={fetchSpecs} title="Reload">
                  <RefreshCw className="h-4 w-4" />
                </Button>
              </CardTitle>
            </CardHeader>
            <CardContent className="space-y-3">
              <div className="space-y-3 max-h-[600px] overflow-y-auto pr-1">
                {sectionConfig && (
                  <div className="space-y-1.5">
                    <div className="grid gap-x-2 gap-y-1.5 grid-cols-[1fr,auto] items-center">
                      {sectionConfig.keys.map(({ key, label, type, optionsKey }) => {
                        const k = paramKey(section, key);
                        const v = (specs.geometry as any)?.[key];
                        const param = typeof v === "object" ? v as ParameterObj : null;
                        const displayValue = inputValues[k] ?? (param ? String(param.value) : String(v ?? ""));
                        const options = (optionsKey && (specs.geometry as any)?.[optionsKey] as (string | number)[]) || [];

                        return (
                          <React.Fragment key={k}>
                            <div className="flex flex-col">
                              <Label className="text-xs">{label}</Label>
                              {param && (
                                <span className="text-[8px] uppercase font-bold text-muted-foreground">{param.type}</span>
                              )}
                            </div>
                            <div className="flex items-center gap-1">
                              {type === "select" ? (
                                <Select value={displayValue} onValueChange={(val) => handleInputChange(section, key, val)}>
                                  <SelectTrigger className="h-7 text-xs w-24"><SelectValue /></SelectTrigger>
                                  <SelectContent>
                                    {options.map((opt) => <SelectItem key={String(opt)} value={String(opt)}>{String(opt)}</SelectItem>)}
                                  </SelectContent>
                                </Select>
                              ) : (
                                <Input value={displayValue} onChange={(e) => handleInputChange(section, key, e.target.value)} className="h-7 text-xs w-24" />
                              )}
                            </div>
                          </React.Fragment>
                        );
                      })}
                    </div>
                  </div>
                )}
              </div>
            </CardContent>
          </Card>

          <Card className="md:col-span-2">
            <CardHeader>
              <CardTitle className="text-base">Machine Cross-Section</CardTitle>
            </CardHeader>
            <CardContent className="flex flex-col items-center">
              <div className="relative border rounded p-4 bg-white shadow-inner">
                <div className="w-[450px] h-[450px]">
                  <MachineCrossSection geometry={specs.geometry} winding={specs.winding} />
                </div>
              </div>
              <div className="w-full mt-4">
                {specs.validations?.[section] && (
                  <Alert variant={specs.validations[section].status === "fail" ? "destructive" : "default"} className={specs.validations[section].status === "warn" ? "border-amber-500 bg-amber-50" : ""}>
                    <AlertTitle>{section.toUpperCase()} Status: {specs.validations[section].status}</AlertTitle>
                    <AlertDescription>
                      {specs.validations[section].errors.map((e, i) => <div key={i}>• {e}</div>)}
                      {specs.validations[section].warnings.map((w, i) => <div key={i}>• {w}</div>)}
                    </AlertDescription>
                  </Alert>
                )}
              </div>
            </CardContent>
          </Card>
          <GeometryMatrix specs={specs} />
        </div>
      );
    }

    // Default view for other steps
    return (
      <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
        <Card>
          <CardHeader>
            <CardTitle className="text-base">{section.toUpperCase()} Specifications</CardTitle>
          </CardHeader>
          <CardContent>
            <div className="grid gap-y-4">
              {sectionConfig?.keys.map(({ key, label, type, optionsKey }) => {
                const k = paramKey(section, key);
                const sectionData = specs[section as keyof MachineSpecs];
                const v = (sectionData as any)?.[key];
                const displayValue = inputValues[k] ?? String(v ?? "");
                const options = (optionsKey && (sectionData as any)?.[optionsKey] as (string | number)[]) || [];

                return (
                  <div key={k} className="flex items-center justify-between">
                    <Label className="text-sm text-muted-foreground">{label}</Label>
                    <div className="w-48">
                      {type === "select" ? (
                        <Select value={displayValue} onValueChange={(val) => handleInputChange(section, key, val)}>
                          <SelectTrigger className="h-9"><SelectValue /></SelectTrigger>
                          <SelectContent>
                            {options.map((opt) => <SelectItem key={String(opt)} value={String(opt)}>{String(opt)}</SelectItem>)}
                          </SelectContent>
                        </Select>
                      ) : (
                        <Input value={displayValue} onChange={(e) => handleInputChange(section, key, e.target.value)} />
                      )}
                    </div>
                  </div>
                );
              })}
            </div>
          </CardContent>
        </Card>

        <Card>
          <CardHeader>
            <CardTitle className="text-base">Technical Validation & Metrics</CardTitle>
          </CardHeader>
          <CardContent className="space-y-4">
            {specs.validations?.[section] ? (
              <div className="space-y-3">
                <div className={`p-4 rounded-lg border ${specs.validations[section].status === "fail" ? "bg-red-50 border-red-200" : specs.validations[section].status === "warn" ? "bg-amber-50 border-amber-200" : "bg-green-50 border-green-200"
                  }`}>
                  <p className="font-bold mb-1">Status: {specs.validations[section].status.toUpperCase()}</p>
                  {specs.validations[section].errors.length > 0 && (
                    <div className="text-sm text-red-700">
                      {specs.validations[section].errors.map((e, i) => <p key={i}>• {e}</p>)}
                    </div>
                  )}
                </div>
                {specs.validations[section].metrics && (
                  <div className="grid grid-cols-1 gap-2">
                    {Object.entries(specs.validations[section].metrics).map(([mK, mV]) => (
                      <div key={mK} className="flex justify-between items-center p-2 rounded bg-slate-50 border">
                        <span className="text-xs font-medium text-slate-500">{mK}</span>
                        <span className="text-xs font-mono font-bold">{String(mV)}</span>
                      </div>
                    ))}
                  </div>
                )}
              </div>
            ) : (
              <p className="text-sm text-muted-foreground">No validation data available for this step.</p>
            )}
          </CardContent>
        </Card>
      </div>
    );
  };

  return (
    <div className="p-6">
      {renderStepContent()}
    </div>
  );
}
const GeometryMatrix = React.memo(({ specs }: { specs: MachineSpecs }) => {
  const geoKeys = PARAM_SECTIONS.find(s => s.section === "geometry")?.keys || [];
  const explorationKeys = geoKeys.filter(k => k.type !== "select" && k.key !== "split_ratio");

  return (
    <Card className="col-span-full mt-6 shadow-sm border-slate-200">
      <CardHeader className="py-3 bg-slate-50 border-b">
        <CardTitle className="text-sm font-bold flex items-center gap-2 text-slate-700">
          <RefreshCw className="h-4 w-4 text-primary" />
          几何尺寸变分矩阵 (±20% Exploration)
        </CardTitle>
      </CardHeader>
      <CardContent className="p-0 overflow-x-auto">
        <table className="w-full border-collapse text-[10px]">
          <thead className="bg-slate-50/50">
            <tr>
              <th className="p-2 border text-left font-semibold text-slate-500 w-32">参数名</th>
              <th className="p-2 border font-semibold text-slate-500 w-20">变分 -20%</th>
              <th className="p-2 border font-semibold text-blue-600 w-24">当前设计</th>
              <th className="p-2 border font-semibold text-slate-500 w-20">变分 +20%</th>
            </tr>
          </thead>
          <tbody>
            {explorationKeys.map(({ key, label }) => {
              const param = (specs.geometry as any)[key] as ParameterObj;
              if (!param || typeof param.value !== "number") return null;

              const baseVal = param.value;
              const variations = [baseVal * 0.8, baseVal, baseVal * 1.2];

              return (
                <tr key={key} className="hover:bg-slate-50/50 transition-colors">
                  <td className="p-2 border font-medium text-slate-700 bg-slate-50/30">
                    <div className="flex flex-col">
                      <span>{label}</span>
                      <span className="text-[9px] text-muted-foreground uppercase">{param.type}</span>
                    </div>
                  </td>
                  {variations.map((v, idx) => {
                    // Create a modified geometry object for preview
                    const modifiedGeo = {
                      ...specs.geometry,
                      [key]: { ...param, value: v }
                    };

                    return (
                      <td key={idx} className={`p-2 border text-center ${idx === 1 ? 'bg-blue-50/30' : ''}`}>
                        <div className="flex flex-col items-center gap-1">
                          <div className="w-16 h-16 border rounded bg-white shadow-inner">
                            <MachineCrossSection
                              geometry={modifiedGeo}
                              winding={specs.winding}
                              showLabels={false}
                            />
                          </div>
                          <span className={idx === 1 ? 'font-bold text-blue-600' : 'text-slate-500'}>
                            {v.toFixed(2)}{param.unit || 'mm'}
                          </span>
                        </div>
                      </td>
                    );
                  })}
                </tr>
              );
            })}
          </tbody>
        </table>
      </CardContent>
    </Card>
  );
});
GeometryMatrix.displayName = "GeometryMatrix";
