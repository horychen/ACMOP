"use client";

import React, { useState, useMemo } from 'react';

const WIRE_DATA = {
    30: { bare: 0.254, coated: 0.290 },
    31: { bare: 0.226, coated: 0.260 },
    32: { bare: 0.203, coated: 0.235 },
};

interface StatorParams {
    OD: number;
    ID: number;
    toothDepth: number;
    toothWidth: number;
    yoke: number;
    liner: number;
}
interface RotorParams {
    OD: number;
    ID: number;
    magnetDepth: number;
    airGap: number;
}
interface WindingParams {
    awg: number;
    J: number;
}
interface MotorParams {
    stator: StatorParams;
    rotor: RotorParams;
    winding: WindingParams;
}

export default function StatorValidator() {
    const [awg, setAwg] = useState<number>(31);
    const [motorParams, setMotorParams] = useState<MotorParams | null>(null);
    const [loading, setLoading] = useState(true);

    React.useEffect(() => {
        // Fetch parameters from backend
        fetch('http://localhost:8001/api/motor-parameters')
            .then(res => res.json())
            .then(data => {
                setMotorParams(data);
                setAwg(data.winding.awg);
                setLoading(false);
            })
            .catch(err => {
                console.error("Failed to fetch motor parameters:", err);
                // Fallback hardcoded values for development if backend is not running
                setMotorParams({
                    stator: { OD: 13.0, ID: 8.3, toothDepth: 2.1, toothWidth: 1.2, yoke: 0.25, liner: 0.1 },
                    rotor: { OD: 8.0, ID: 2.0, magnetDepth: 3.0, airGap: 0.15 },
                    winding: { awg: 31, J: 14 }
                });
                setLoading(false);
            });
    }, []);

    const scale = 200; // 1mm = 200px

    // 2. Derived Calculations
    const metrics = useMemo(() => {
        if (!motorParams) return null;

        // 1. Core Geometric Parameters (mm) from Backend
        const OD = motorParams.stator.OD;
        const ID = motorParams.stator.ID;
        const toothDepth = motorParams.stator.toothDepth;
        const toothWidth = motorParams.stator.toothWidth;
        const yoke = motorParams.stator.yoke;
        const liner = motorParams.stator.liner;
        const J = motorParams.winding.J; // A/mm²

        const rInner = ID / 2; // 4.15 mm
        const rOuter = OD / 2 - yoke; // 6.25 mm

        // Slot width calculations (Parallel teeth)
        const slotOpeningNet = ((Math.PI * ID) / 12) - toothWidth - (2 * liner); // ~0.77 mm
        const slotBottomNet = ((Math.PI * (rOuter * 2)) / 12) - toothWidth - (2 * liner); // ~1.87 mm

        // Trapezoidal Slot Area Approximation (Net of liner)
        const netArea = ((slotOpeningNet + slotBottomNet) / 2) * (toothDepth - liner); // ~2.77 mm^2 Net space

        // Wire properties
        const wireData = WIRE_DATA[awg as keyof typeof WIRE_DATA];
        const dCoated = wireData?.coated || 0.260;
        const aBare = Math.PI * Math.pow((wireData?.bare || 0.226) / 2, 2);

        // Orthocyclic packing estimation (Max physical wires per half-slot)
        // Assuming double layer, we pack one side of the tooth
        const effectiveWidth = slotOpeningNet / 2; // Split for two coils
        const layers = Math.floor((toothDepth - liner * 2) / (dCoated * 0.866)); // Hexagonal stacking height
        const wiresPerLayer = Math.max(1, Math.floor(effectiveWidth / dCoated));

        // Adjusting Z_slot based on realistic orthocyclic packing factor (0.88)
        const rawCapacity = (netArea * 0.88) / (Math.PI * Math.pow(dCoated / 2, 2));
        const zSlot = Math.floor(rawCapacity * 0.7); // Heuristic for realistic winding limits and double layer

        const totalCopperArea = zSlot * aBare;
        const fillFactor = totalCopperArea / 3.20; // Using geometric gross area ~3.20 mm^2

        const I = aBare * J;
        const NI = zSlot * I;

        // Thermal Load (Electric Loading A * Current Density J)
        const totalConductors = 12 * zSlot;
        const A_loading = (totalConductors * I) / (Math.PI * ID); // A/mm
        const AJ = (A_loading * 10) * J; // Converted to A^2/(cm*mm^2)

        return { slotOpeningNet, slotBottomNet, netArea, zSlot, fillFactor, I, NI, AJ, dCoated };
    }, [awg, motorParams]);

    // Handle loading state
    if (loading || !motorParams || !metrics) {
        return <div className="p-6 text-center text-slate-500">Loading Motor Parameters...</div>;
    }

    const { OD, ID, toothDepth, toothWidth, yoke, liner } = motorParams.stator;
    const { magnetDepth, airGap } = motorParams.rotor;
    const { J } = motorParams.winding;

    // 3. SVG Rendering Helpers
    const renderWires = () => {
        const wires: React.ReactNode[] = [];
        const { dCoated } = metrics;
        const slots = 12; // 12S
        const rIn = ID / 2;
        const slotAngle = (2 * Math.PI) / slots;

        // Two adjacent teeth forming the slot
        const toothAngles = [-slotAngle / 2, slotAngle / 2];

        toothAngles.forEach((tAngle, idx) => {
            const side = idx === 0 ? 1 : -1;
            const fillC = idx === 0 ? "#B87333" : "#CD7F32";
            const strokeC = "#8B4513";

            for (let layer = 0; layer < 6; layer++) {
                // Horizontal offset from tooth center 
                const hOffset = (toothWidth / 2 + liner + dCoated / 2 + layer * dCoated * 0.88);

                for (let row = 0; row < 15; row++) {
                    const vOffset = liner + dCoated / 2 + row * dCoated;
                    if (vOffset > toothDepth - liner) break; // Reached bottom of slot

                    const r = rIn + vOffset;
                    if (hOffset >= r) continue;

                    const thetaOffset = Math.asin(hOffset / r) * side;
                    const finalAngle = tAngle + thetaOffset;

                    // Slot center collision check
                    if (side === 1 && finalAngle > -0.002) continue;
                    if (side === -1 && finalAngle < 0.002) continue;

                    // Convert to StatorValidator's Cartesian mapping
                    // y axis is depth (r - rIn), which is exactly vOffset
                    // x axis is distance from slot center, which is r * sin(finalAngle)
                    const cx = r * Math.sin(finalAngle) * scale;
                    const cy = vOffset * scale;

                    wires.push(
                        <circle
                            key={`W-${idx}-${layer}-${row}`} cx={cx} cy={cy} r={dCoated / 2 * scale}
                            fill={fillC} stroke={strokeC} strokeWidth="1"
                        />
                    );
                }
            }
        });

        return wires;
    };

    return (
        <div className="p-6 max-w-4xl mx-auto bg-slate-50 rounded-xl shadow-lg font-sans">
            <h2 className="text-2xl font-bold mb-6 text-slate-800 border-b pb-2">12S10P Micro Motor Validation Dashboard</h2>

            <div className="grid grid-cols-1 md:grid-cols-2 gap-8">
                {/* Visualizer Section */}
                <div className="bg-white p-4 rounded-lg shadow border border-slate-200 flex flex-col items-center overflow-hidden">
                    <h3 className="text-lg font-semibold mb-4 text-slate-700">Slot Cross-Section (1mm = 200px)</h3>

                    <svg width="400" height="500" viewBox="-200 -200 400 600" className="bg-slate-100 rounded">
                        {/* Rotor Core (Background) */}
                        <rect x="-200" y="-200" width="400" height={(200 / scale - airGap - magnetDepth) * scale} fill="#cbd5e1" />
                        <text x="-190" y={- (airGap + magnetDepth) * scale - 20} fontSize="12" fill="#475569">Rotor Core</text>

                        {/* Magnet */}
                        <rect x="-200" y={-(airGap + magnetDepth) * scale} width="400" height={magnetDepth * scale} fill="#fca5a5" stroke="#ef4444" strokeWidth="2" />
                        <text x="-190" y={-(airGap + magnetDepth / 2) * scale + 5} fontSize="12" fill="#7f1d1d" fontWeight="bold">Magnet ({magnetDepth}mm)</text>

                        {/* Airgap */}
                        <rect x="-200" y={-airGap * scale} width="400" height={airGap * scale} fill="#e2e8f0" />
                        <text x="-190" y={-airGap * scale / 2 + 5} fontSize="12" fill="#64748b">Airgap ({airGap}mm)</text>

                        {/* Stator Iron (Teeth & Yoke) - Approximated Trapezoidal View */}
                        <path
                            d={`M -200,0 L -${metrics.slotOpeningNet * scale / 2 + (liner * scale)},0 
                  L -${metrics.slotBottomNet * scale / 2 + (liner * scale)},${toothDepth * scale} 
                  L -200,${toothDepth * scale} Z`}
                            fill="#94a3b8"
                        />
                        <path
                            d={`M 200,0 L ${metrics.slotOpeningNet * scale / 2 + (liner * scale)},0 
                  L ${metrics.slotBottomNet * scale / 2 + (liner * scale)},${toothDepth * scale} 
                  L 200,${toothDepth * scale} Z`}
                            fill="#94a3b8"
                        />
                        <rect x="-200" y={toothDepth * scale} width="400" height={yoke * scale} fill="#64748b" />

                        {/* Slot Liner */}
                        <path
                            d={`M -${metrics.slotOpeningNet * scale / 2},0 
                  L -${metrics.slotBottomNet * scale / 2},${toothDepth * scale} 
                  L ${metrics.slotBottomNet * scale / 2},${toothDepth * scale} 
                  L ${metrics.slotOpeningNet * scale / 2},0`}
                            fill="none" stroke="#fef08a" strokeWidth={liner * scale}
                        />

                        {/* Centerline Divider (Double Layer Separation) */}
                        <line x1="0" y1="0" x2="0" y2={toothDepth * scale} stroke="#cbd5e1" strokeWidth="2" strokeDasharray="5,5" />

                        {/* Winding Rendering */}
                        {renderWires()}
                    </svg>
                </div>

                {/* Engineering Metrics Section */}
                <div className="space-y-6">
                    <div className="bg-white p-4 rounded-lg shadow border border-slate-200">
                        <h3 className="text-lg font-semibold mb-3 text-slate-700">Wire Selection</h3>
                        <select
                            value={awg}
                            onChange={(e) => setAwg(Number(e.target.value))}
                            className="w-full p-2 border rounded bg-slate-50 font-mono text-lg"
                        >
                            <option value="30">AWG 30 (0.254 mm bare)</option>
                            <option value="31">AWG 31 (0.226 mm bare)</option>
                            <option value="32">AWG 32 (0.203 mm bare)</option>
                        </select>
                    </div>

                    <div className="bg-white p-4 rounded-lg shadow border border-slate-200">
                        <h3 className="text-lg font-semibold mb-4 text-slate-700">Engineering Metrics (J = 14 A/mm²)</h3>
                        <div className="space-y-3 font-mono text-sm">
                            <div className="flex justify-between border-b pb-1">
                                <span className="text-slate-600">Total Conductors (Z slot):</span>
                                <span className="font-bold">{metrics.zSlot} wires</span>
                            </div>
                            <div className="flex justify-between border-b pb-1">
                                <span className="text-slate-600">Copper Fill Factor (Cu):</span>
                                <span className={`font-bold ${(metrics.fillFactor * 100) > 45 ? 'text-red-500' : 'text-green-600'}`}>
                                    {(metrics.fillFactor * 100).toFixed(1)}%
                                </span>
                            </div>
                            <div className="flex justify-between border-b pb-1">
                                <span className="text-slate-600">Current per Wire ($I$):</span>
                                <span className="font-bold">{metrics.I.toFixed(3)} A</span>
                            </div>
                            <div className="flex justify-between border-b pb-1">
                                <span className="text-slate-600">Total Ampere-Turns ($NI$):</span>
                                <span className="font-bold">{metrics.NI.toFixed(1)} AT</span>
                            </div>
                            <div className="flex justify-between border-b pb-1 bg-red-50 p-2 rounded">
                                <span className="text-red-700 font-bold">Thermal Load ($AJ$):</span>
                                <span className="font-bold text-red-700">{metrics.AJ.toFixed(0)} A²/(cm·mm²)</span>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    );
}