"use client";

import React, { useEffect, useState, useMemo, useRef, useCallback } from "react";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/alert";
import { Loader2, AlertTriangle, RefreshCw } from "lucide-react";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";

/** Sections and param keys aligned with backend/codes4/user_minitureMachine.py (GeometrySpecs, WindingSpecs, MaterialSpecs, PerformanceTargets) */
const PARAM_SECTIONS: {
  section: "geometry" | "winding" | "materials" | "targets";
  title: string;
  keys: { key: string; label: string; type?: "number" | "text" }[];
}[] = [
  {
    section: "geometry",
    title: "Geometry",
    keys: [
      { key: "d_stator_outer", label: "Stator OD (mm)" },
      { key: "d_rotor_outer", label: "Rotor OD (mm)" },
      { key: "d_shaft", label: "Shaft (mm)" },
      { key: "air_gap", label: "Air gap (mm)" },
      { key: "magnet_thickness", label: "Magnet thick. (mm)" },
      { key: "tooth_width", label: "Tooth width (mm)" },
      { key: "tooth_depth", label: "Tooth depth (mm)" },
      { key: "stack_length", label: "Stack length (mm)" },
      { key: "nominal_eccentricity", label: "Eccentricity (mm)" },
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
      { key: "stator_steel", label: "Stator steel", type: "text" },
      { key: "steel_thickness", label: "Steel thick. (mm)" },
      { key: "steel_stack_factor", label: "Stack factor" },
      { key: "steel_max_flux_density", label: "B_sat (T)" },
    ],
  },
  {
    section: "targets",
    title: "Performance",
    keys: [
      { key: "magnetic_loading_target", label: "B_g target (T)" },
      { key: "maxwell_stress_radial", label: "Maxwell stress (kPa)" },
      { key: "torque_constant_kt", label: "K_t (mNm/A)" },
      { key: "ump_threshold_max", label: "UMP max (N)" },
      { key: "temp_limit", label: "Temp limit (°C)" },
    ],
  },
];

const PARAM_KEYS = PARAM_SECTIONS.flatMap(({ section, keys }) =>
  keys.map((k) => ({ section, key: k.key, label: k.label, type: k.type ?? "number" }))
);

function paramKey(section: string, key: string) {
  return `${section}.${key}`;
}

function specsToInputValues(specs: MachineSpecs): Record<string, string> {
  const out: Record<string, string> = {};
  for (const { section, key } of PARAM_KEYS) {
    const sec = specs[section as keyof MachineSpecs];
    if (!sec || typeof sec !== "object") continue;
    const v = (sec as Record<string, unknown>)[key];
    if (v !== undefined && v !== null) out[paramKey(section, key)] = String(v);
  }
  return out;
}

/** Matches backend MotorSpecs (GeometrySpecs, WindingSpecs, MaterialSpecs, PerformanceTargets) */
interface MachineSpecs {
  geometry: {
    d_stator_outer: number;
    d_rotor_outer: number;
    d_shaft: number;
    air_gap: number;
    magnet_thickness: number;
    tooth_width: number;
    tooth_depth: number;
    stack_length_options?: number[];
    stack_length?: number;
    nominal_eccentricity?: number;
  };
  winding: {
    num_slots: number;
    num_poles: number;
    coil_pitch_y?: number;
    conductors_per_slot: number;
    wire_diameter_with_insulation: number;
    rated_current_density?: number;
  };
  materials: {
    magnet_grade: string;
    magnet_br: number;
    magnet_h_cj: number;
    magnet_temp_max: number;
    stator_steel: string;
    steel_thickness: number;
    steel_stack_factor: number;
    steel_max_flux_density: number;
  };
  targets: {
    magnetic_loading_target: number;
    maxwell_stress_radial: number;
    torque_constant_kt: number;
    ump_threshold_max: number;
    temp_limit: number;
  };
}

const COUNTDOWN_SEC = 3;

export default function MachineVisualizer() {
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
        // #region agent log
        fetch("http://127.0.0.1:7242/ingest/e5770a5a-cc34-4592-9d22-bc6eef1eb00c", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ hypothesisId: "H2-H3-H5", location: "MachineVisualizer.tsx:fetchSpecs.then", message: "API response geometry", data: { api_tooth_width: data?.geometry?.tooth_width, api_tooth_depth: data?.geometry?.tooth_depth, geometry_keys: data?.geometry ? Object.keys(data.geometry) : [] }, timestamp: Date.now() }) }).catch(() => {});
        // #endregion
        const normalized = {
          geometry: { ...data?.geometry },
          winding: { ...data?.winding },
          materials: { ...data?.materials },
          targets: { ...data?.targets },
        } as MachineSpecs;
        // #region agent log
        const inputV = specsToInputValues(normalized);
        fetch("http://127.0.0.1:7242/ingest/e5770a5a-cc34-4592-9d22-bc6eef1eb00c", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ hypothesisId: "H2-H3-H5", location: "MachineVisualizer.tsx:fetchSpecs.then", message: "normalized and inputValues", data: { norm_tooth_width: normalized.geometry?.tooth_width, norm_tooth_depth: normalized.geometry?.tooth_depth, input_geometry_tooth_width: inputV["geometry.tooth_width"], input_geometry_tooth_depth: inputV["geometry.tooth_depth"] }, timestamp: Date.now() }) }).catch(() => {});
        // #endregion
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
    if (specs && Object.keys(inputValues).length === 0) {
      setInputValues(specsToInputValues(specs));
    }
  }, [specs]);

  const visualizationData = useMemo(() => {
    if (!specs) return null;

    const { geometry, winding } = specs;
    
    const R_so = geometry.d_stator_outer / 2;
    const R_ro = geometry.d_rotor_outer / 2;
    const gap = geometry.air_gap;
    const R_si = R_ro + gap;
    const toothDepth = geometry.tooth_depth;
    const R_sb = R_si + toothDepth; // Slot bottom (outer radius of slot)
    const toothWidth = geometry.tooth_width;
    const numSlots = winding.num_slots;
    const numPoles = winding.num_poles;
    const conductorsPerSlot = winding.conductors_per_slot;
    const wireDiam = winding.wire_diameter_with_insulation;
    const wireRad = wireDiam / 2;
    const magnetThickness = geometry.magnet_thickness;
    const R_shaft = geometry.d_shaft / 2;

    // Helper to get point at radius r, angle theta (degrees)
    const getPt = (r: number, theta: number) => {
      const rad = (theta * Math.PI) / 180;
      return { x: r * Math.cos(rad), y: r * Math.sin(rad) };
    };

    // --- Magnets & Back Iron ---
    const R_mag_in = R_ro - magnetThickness;
    const hasBackIron = R_mag_in > R_shaft + 0.01;
    
    const magnetPaths: {d: string, color: string}[] = [];
    if (numPoles > 0) {
        const poleAngle = 360 / numPoles;
        for (let i = 0; i < numPoles; i++) {
            const startAngle = i * poleAngle;
            const endAngle = (i + 1) * poleAngle;
            
            const p1 = getPt(R_ro, startAngle);
            const p2 = getPt(R_ro, endAngle);
            const p3 = getPt(R_mag_in, endAngle);
            const p4 = getPt(R_mag_in, startAngle);
            
            let d = `M ${p1.x} ${p1.y} `;
            d += `A ${R_ro} ${R_ro} 0 0 1 ${p2.x} ${p2.y} `;
            d += `L ${p3.x} ${p3.y} `;
            d += `A ${R_mag_in} ${R_mag_in} 0 0 0 ${p4.x} ${p4.y} `;
            d += `Z`;
            
            magnetPaths.push({
                d,
                color: i % 2 === 0 ? "#fca5a5" : "#93c5fd"
            });
        }
    }

    // --- Stator Geometry ---
    const getToothPoints = (angle: number, r: number, w: number) => {
        const rad = (angle * Math.PI) / 180;
        const u = { x: Math.cos(rad), y: Math.sin(rad) };
        const n = { x: -Math.sin(rad), y: Math.cos(rad) }; // Normal pointing CCW
        
        const h = w / 2;
        if (r < h) return null;
        
        const d = Math.sqrt(r*r - h*h);
        
        const p1 = { x: d * u.x + h * n.x, y: d * u.y + h * n.y };
        const p2 = { x: d * u.x - h * n.x, y: d * u.y - h * n.y };
        return { cw: p2, ccw: p1 };
    };

    let d_path = "";
    for (let i = 0; i < numSlots; i++) {
        const angle = i * 360 / numSlots;
        const tip = getToothPoints(angle, R_si, toothWidth);
        const root = getToothPoints(angle, R_sb, toothWidth);
        
        if (!tip || !root) continue;
        
        const nextAngle = (i + 1) * 360 / numSlots;
        const nextRoot = getToothPoints(nextAngle, R_sb, toothWidth);
        const nextTip = getToothPoints(nextAngle, R_si, toothWidth);
        
        if (!nextRoot) continue;

        if (i === 0) {
            d_path += `M ${tip.cw.x} ${tip.cw.y} `;
        }
        
        d_path += `A ${R_si} ${R_si} 0 0 1 ${tip.ccw.x} ${tip.ccw.y} `;
        d_path += `L ${root.ccw.x} ${root.ccw.y} `;
        d_path += `A ${R_sb} ${R_sb} 0 0 1 ${nextRoot.cw.x} ${nextRoot.cw.y} `;
        d_path += `L ${nextTip!.cw.x} ${nextTip!.cw.y} `;
    }
    d_path += "Z";
    
    const fullStatorPath = `M ${R_so} 0 A ${R_so} ${R_so} 0 1 1 ${-R_so} 0 A ${R_so} ${R_so} 0 1 1 ${R_so} 0 Z ` + d_path;

    // --- Conductors: layer-by-layer winding (贴齿一层一层向外绕线) ---
    const halfTooth = toothWidth / 2;
    const slotAngle = 360 / numSlots;
    const halfSlotAngle = slotAngle / 2;
    const halfSlotRad = (halfSlotAngle * Math.PI) / 180;
    const angleStep = 360 / numSlots;
    const radStep = (angleStep * Math.PI) / 180;
    const n2 = { x: -Math.sin(radStep), y: Math.cos(radStep) }; // Normal to Tooth 1 ray

    const padding = 0.01;
    const effR = wireRad + padding;
    const conductorsPerSide = Math.ceil(conductorsPerSlot / 2);
    const hexRowDy = effR * Math.sqrt(3);

    /** 绕线工艺：第一层导体相切排列，第二层及以后与上一层交错嵌套（六边形密排），体现机械臂紧密绕制。 */
    const packSideLayers = (isRightSide: boolean): {x: number, y: number}[] => {
        const placed: {x: number, y: number}[] = [];
        const splitRad = halfSlotRad;
        const tanSplit = Math.tan(splitRad);

        if (isRightSide) {
            // 左侧半槽贴 Tooth 0。第一层 y = halfTooth + effR（贴齿相切），层间距 sqrt(3)*effR，奇偶层 x 错位 effR 形成交错嵌套。
            let layerIndex = 0;
            while (placed.length < conductorsPerSide) {
                const yLayer = halfTooth + effR + layerIndex * hexRowDy;
                if (yLayer > R_sb - effR) break;
                const rInnerSq = (R_si + effR) ** 2 - yLayer * yLayer;
                const rOuterSq = (R_sb - effR) ** 2 - yLayer * yLayer;
                if (rInnerSq < 0 || rOuterSq < 0) { layerIndex++; continue; }
                const xInner = Math.sqrt(rInnerSq);
                const xOuter = Math.sqrt(rOuterSq);
                const xSplit = yLayer / tanSplit;
                const xMin = Math.max(xInner, xSplit + effR);
                const xMax = xOuter;
                if (xMin >= xMax - effR) { layerIndex++; continue; }
                const rowOffset = (layerIndex % 2) * effR;
                let x = xMin + effR + rowOffset;
                const xStep = 2 * effR;
                while (x <= xMax - effR && placed.length < conductorsPerSide) {
                    const r = Math.sqrt(x * x + yLayer * yLayer);
                    if (r >= R_si + effR && r <= R_sb - effR && x > xSplit + effR) {
                        placed.push({ x, y: yLayer });
                    }
                    x += xStep;
                }
                layerIndex++;
            }
        } else {
            // 右侧半槽贴 Tooth 1。第一层线 p·n2 = -halfTooth - effR，之后每层向槽内移 sqrt(3)*effR，奇偶层沿齿向错位 effR。
            const u1 = { x: Math.cos(radStep), y: Math.sin(radStep) };
            let layerIndex = 0;
            while (placed.length < conductorsPerSide) {
                const lineVal = -halfTooth - effR - layerIndex * hexRowDy;
                const d = lineVal;
                const tMin = R_si + effR;
                const tMax = R_sb - effR;
                const tRange = tMax - tMin;
                if (tRange < 2 * effR) { layerIndex++; continue; }
                const rowOffset = (layerIndex % 2) * effR;
                for (let j = 0; placed.length < conductorsPerSide; j++) {
                    const t = tMin + effR + rowOffset + j * 2 * effR;
                    if (t > tMax - effR) break;
                    const px = t * u1.x + d * n2.x;
                    const py = t * u1.y + d * n2.y;
                    const r = Math.sqrt(px * px + py * py);
                    if (r < R_si + effR || r > R_sb - effR) continue;
                    const angleP = Math.atan2(py, px);
                    if (angleP <= splitRad + 0.001) continue;
                    if (angleP >= radStep - 0.001) continue;
                    let overlap = false;
                    for (const ex of placed) {
                        const dx = px - ex.x, dy = py - ex.y;
                        if (dx * dx + dy * dy < (2 * effR) ** 2 - 0.0001) { overlap = true; break; }
                    }
                    if (!overlap) placed.push({ x: px, y: py });
                }
                layerIndex++;
                if (Math.abs(lineVal) > R_sb + 1) break;
            }
        }
        return placed;
    };

    const packSide = (isRightSide: boolean): {x: number, y: number}[] => {
        const byLayer = packSideLayers(isRightSide);
        if (byLayer.length >= conductorsPerSide) return byLayer.slice(0, conductorsPerSide);
        return byLayer;
    };

    const leftConductors = packSide(true);
    const rightConductors = packSide(false);
    const slot0Conductors = [...leftConductors, ...rightConductors];
    
    const totalFit = leftConductors.length + rightConductors.length;
    let currentErrorMsg = null;
    let currentMaxFit = totalFit;
    
    if (leftConductors.length < conductorsPerSide || rightConductors.length < conductorsPerSide) {
        currentErrorMsg = `Cannot fit ${conductorsPerSlot} conductors (Need ${conductorsPerSide} per side).`;
    }

    const finalAllConductors: {x: number, y: number}[] = [];
    const slotSplitLines: {x1: number, y1: number, x2: number, y2: number}[] = [];
    
    for (let i = 0; i < numSlots; i++) {
        const angle = i * 360 / numSlots;
        const rad = angle * Math.PI / 180;
        const cos = Math.cos(rad);
        const sin = Math.sin(rad);
        
        slot0Conductors.forEach(c => {
            finalAllConductors.push({
                x: c.x * cos - c.y * sin,
                y: c.x * sin + c.y * cos
            });
        });
        
        const splitA = angle + halfSlotAngle;
        const pStart = getPt(R_si, splitA);
        const pEnd = getPt(R_sb, splitA);
        slotSplitLines.push({
            x1: pStart.x, y1: pStart.y,
            x2: pEnd.x, y2: pEnd.y
        });
    }

    return {
        R_so, R_si, R_ro, R_sb, R_shaft, R_mag_in, hasBackIron,
        fullStatorPath,
        allConductors: finalAllConductors,
        magnetPaths,
        slotSplitLines,
        errorMsg: currentErrorMsg,
        maxFit: currentMaxFit,
        conductorsPerSlot
    };
  }, [specs]);

  const apiStatus =
    error
      ? `请求失败: ${error}（请确认后端已启动且为 http://localhost:8000）`
      : apiDebug
        ? `API: tooth_width_from_py=${String(apiDebug.tooth_width_from_py)} path=${apiDebug.path_used ?? ""}`
        : specs
          ? `API 响应无 _debug，当前 geometry.tooth_width=${specs.geometry?.tooth_width}（可能未打到 FastAPI 或后端未返回 _debug）`
          : loading
            ? "请求中..."
            : "无数据";

  if (loading) return <div className="flex flex-col items-center gap-4 p-8"><Loader2 className="animate-spin" /><p className="text-sm text-muted-foreground">API 状态: {apiStatus}</p></div>;
  if (error) return (
    <div className="space-y-2 p-4">
      <Alert variant="destructive"><AlertTitle>Error</AlertTitle><AlertDescription>{error}</AlertDescription></Alert>
      <p className="text-sm text-muted-foreground">请求 URL: /api/machine-specs（代理到后端 8000）— 请确认后端已启动。</p>
    </div>
  );
  if (!specs || !visualizationData) return <div className="p-4">No data. API 状态: {apiStatus}</div>;

  const { R_so, fullStatorPath, allConductors, magnetPaths, slotSplitLines, errorMsg, maxFit, conductorsPerSlot, R_ro, R_shaft, R_mag_in, hasBackIron } = visualizationData;
  const viewBoxSize = R_so * 2.4;
  
  return (
    <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
      {/* Controls */}
      <Card className="md:col-span-1">
        <CardHeader>
            <CardTitle className="flex items-center justify-between">
                Parameters
                <Button variant="outline" size="icon" onClick={fetchSpecs} title="Reload from File">
                    <RefreshCw className="h-4 w-4" />
                </Button>
            </CardTitle>
        </CardHeader>
        <CardContent className="space-y-3">
            <div className="space-y-3 max-h-[600px] overflow-y-auto pr-1">
              {PARAM_SECTIONS.map(({ section, title, keys }) => {
                const secSpecs = specs?.[section as keyof MachineSpecs] as Record<string, unknown> | undefined;
                const secInitial = initialSpecs?.[section as keyof MachineSpecs] as Record<string, unknown> | undefined;
                return (
                  <div key={section} className="space-y-1.5">
                    <p className="text-xs font-medium text-muted-foreground uppercase tracking-wide sticky top-0 bg-card py-0.5">{title}</p>
                    <div className="grid gap-x-2 gap-y-1.5 grid-cols-[1fr,auto] items-center">
                      {keys.map(({ key, label, type }) => {
                        const k = paramKey(section, key);
                        const displayValue = inputValues[k] ?? (secSpecs ? String(secSpecs[key] ?? "") : "");
                        const defaultVal = secInitial?.[key];
                        const isPending = pendingUpdate?.section === section && pendingUpdate?.key === key;
                        return (
                          <React.Fragment key={k}>
                            <div className="flex items-center justify-between gap-1 min-w-0">
                              <Label className="text-xs shrink-0">{label}</Label>
                              {defaultVal !== undefined && defaultVal !== null && (
                                <span className="text-[10px] text-muted-foreground truncate">默认: {String(defaultVal)}</span>
                              )}
                            </div>
                            <div className="flex items-center gap-1">
                              <Input
                                type={type === "text" ? "text" : "text"}
                                inputMode={type === "text" ? "text" : "decimal"}
                                value={displayValue}
                                onChange={(e) => handleInputChange(section, key, e.target.value)}
                                className="h-7 text-xs w-20 shrink-0"
                              />
                              {isPending && countdownSec > 0 && (
                                <span className="text-[10px] text-muted-foreground whitespace-nowrap">{countdownSec}s</span>
                              )}
                            </div>
                          </React.Fragment>
                        );
                      })}
                    </div>
                  </div>
                );
              })}
            </div>
        </CardContent>
      </Card>

      {/* Visualization */}
      <Card className="md:col-span-2">
        <CardHeader>
            <CardTitle>Cross-Section View</CardTitle>
        </CardHeader>
        <CardContent className="flex flex-col items-center">
            <div className="relative border rounded p-4 bg-white">
                <svg width="600" height="600" viewBox={`${-viewBoxSize/2} ${-viewBoxSize/2} ${viewBoxSize} ${viewBoxSize}`}>
                    {/* Scale Bar - 1mm */}
                    <line x1={-viewBoxSize/2 + 1} y1={viewBoxSize/2 - 2} x2={-viewBoxSize/2 + 1 + 1} y2={viewBoxSize/2 - 2} stroke="black" strokeWidth="0.1" />
                    <text x={-viewBoxSize/2 + 1} y={viewBoxSize/2 - 2.5} fontSize="0.5" fill="black">1mm</text>

                    {/* Stator Core - Grey */}
                    <path d={fullStatorPath} fill="#9ca3af" stroke="#4b5563" strokeWidth="0.05" fillRule="evenodd" />
                    
                    {/* Rotor Back Iron (if exists) - Grey */}
                    {hasBackIron && (
                        <circle cx="0" cy="0" r={R_mag_in} fill="#9ca3af" stroke="#4b5563" strokeWidth="0.05" />
                    )}

                    {/* Magnets */}
                    {magnetPaths.map((m, i) => (
                        <path key={i} d={m.d} fill={m.color} stroke="none" />
                    ))}
                    
                    {/* Shaft - White */}
                    <circle cx="0" cy="0" r={R_shaft} fill="white" stroke="#4b5563" strokeWidth="0.05" />
                    
                    {/* Slot Split Lines (Dash-Dot) */}
                    {slotSplitLines.map((l, i) => (
                        <line 
                            key={`split-${i}`} 
                            x1={l.x1} y1={l.y1} 
                            x2={l.x2} y2={l.y2} 
                            stroke="#94a3b8" 
                            strokeWidth="0.05" 
                            strokeDasharray="0.5 0.2 0.1 0.2" 
                        />
                    ))}

                    {/* Conductors */}
                    {allConductors.map((c, i) => (
                        <circle key={i} cx={c.x} cy={c.y} r={specs.winding.wire_diameter_with_insulation/2} fill="#b91c1c" stroke="none" />
                    ))}
                </svg>
            </div>
            
            {errorMsg && (
                <Alert variant="destructive" className="mt-4">
                    <AlertTriangle className="h-4 w-4" />
                    <AlertTitle>Fitting Error</AlertTitle>
                    <AlertDescription>{errorMsg}</AlertDescription>
                </Alert>
            )}
            
            <div className="mt-4 text-sm text-muted-foreground">
                <p className="text-xs mt-2 border-b pb-2 text-amber-700 bg-amber-50/80 rounded px-2 py-1" title="后端调试信息">
                    API 状态: {apiStatus}
                </p>
                <p>Stator OD: {specs.geometry.d_stator_outer}mm</p>
                <p>Conductors/Slot: {conductorsPerSlot}</p>
                <p>Wire Diameter: {specs.winding.wire_diameter_with_insulation}mm</p>
                {errorMsg && <p className="font-bold text-red-500">Max Fit: {maxFit}</p>}
            </div>
        </CardContent>
      </Card>
    </div>
  );
}
