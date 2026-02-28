import React, { useState, useMemo } from 'react';

const WIRE_DATA = {
    30: { bare: 0.254, coated: 0.290 },
    31: { bare: 0.226, coated: 0.260 },
    32: { bare: 0.203, coated: 0.235 },
};

export default function StatorValidator() {
    const [awg, setAwg] = useState(31);
    const [motorParams, setMotorParams] = useState(null);
    const [loading, setLoading] = useState(true);

    React.useEffect(() => {
        // Fetch parameters from backend
        fetch('http://localhost:8000/api/motor-parameters')
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
                    stator: { OD: 13.0, ID: 8.3, toothDepth: 2.1, toothWidth: 1.2, yoke: 0.25 },
                    winding: { awg: 31, J: 14, liner: 0.1 }
                });
                setLoading(false);
            });
    }, []);

    const scale = 200; // 1mm = 200px

    // Handle loading state
    if (loading || !motorParams) {
        return <div className="p-6 text-center text-slate-500">Loading Motor Parameters...</div>;
    }

    // 1. Core Geometric Parameters (mm) from Backend
    const OD = motorParams.stator.OD;
    const ID = motorParams.stator.ID;
    const toothDepth = motorParams.stator.toothDepth;
    const toothShoe = motorParams.stator.toothShoe || 0.15;
    const toothWidth = motorParams.stator.toothWidth;
    const yoke = motorParams.stator.yoke;
    const liner = motorParams.winding.liner;
    const J = motorParams.winding.J; // A/mm²

    // 2. Derived Calculations
    const metrics = useMemo(() => {
        const rInner = ID / 2; // 4.15 mm
        const rOuter = OD / 2 - yoke; // 6.25 mm

        // Slot width calculations (Parallel teeth)
        const slotOpeningNet = ((Math.PI * ID) / 12) - toothWidth - (2 * liner); // ~0.77 mm
        const slotBottomNet = ((Math.PI * (rOuter * 2)) / 12) - toothWidth - (2 * liner); // ~1.87 mm

        // Trapezoidal Slot Area Approximation (Net of liner)
        const dSlot = toothDepth - toothShoe;
        const netArea = ((slotOpeningNet + slotBottomNet) / 2) * (dSlot - liner); // ~2.77 mm^2 Net space

        // Wire properties
        const dCoated = WIRE_DATA[awg].coated;
        const aBare = Math.PI * Math.pow(WIRE_DATA[awg].bare / 2, 2);

        // Orthocyclic packing estimation (Max physical wires per half-slot)
        // Assuming double layer, we pack one side of the tooth
        const effectiveWidth = slotOpeningNet / 2; // Split for two coils
        const layers = Math.floor((dSlot - liner * 2) / (dCoated * 0.866)); // Hexagonal stacking height
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

        return { slotOpeningNet, netArea, zSlot, fillFactor, I, NI, AJ, dCoated };
    }, [awg]);

    // 3. SVG Rendering Helpers
    const renderWires = () => {
        const wires = [];
        const { dCoated, zSlot } = metrics;
        const r = (dCoated * scale) / 2;

        // Simulate wire packing (Right side of the slot - single coil side)
        let currentX = r + (liner * scale);
        let currentY = toothDepth * scale - r - (liner * scale);
        let rowCount = 0;

        for (let i = 0; i < zSlot / 2; i++) {
            wires.push(
                <circle
                    key={`R-${i}`} cx={currentX} cy={currentY} r={r}
                    fill="#B87333" stroke="#8B4513" strokeWidth="2"
                />
            );
            // Mirror for the left side (second coil)
            wires.push(
                <circle
                    key={`L-${i}`} cx={-currentX} cy={currentY} r={r}
                    fill="#CD7F32" stroke="#8B4513" strokeWidth="2"
                />
            );

            currentX += dCoated * scale;
            // Stacking logic with collision threshold (simplified)
            if (currentX > (metrics.slotBottomNet * scale / 2) - r || rowCount >= 4) {
                currentY -= (dCoated * 0.866 * scale); // Orthocyclic offset
                rowCount = 0;
                currentX = r + (liner * scale) + (currentY % 2 === 0 ? r : 0);
            } else {
                rowCount++;
            }
        }
        return wires;
    };

    return (
        <div className="p-6 max-w-4xl mx-auto bg-slate-50 rounded-xl shadow-lg font-sans">
            <h2 className="text-2xl font-bold mb-6 text-slate-800 border-b pb-2">12S10P Micro Motor Validation Dashboard</h2>

            <div className="grid grid-cols-1 md:grid-cols-2 gap-8">
                {/* Visualizer Section */}
                <div className="bg-white p-4 rounded-lg shadow border border-slate-200 flex flex-col items-center overflow-hidden">
                    <h3 className="text-lg font-semibold mb-4 text-slate-700">Slot Cross-Section (1mm = 200px)</h3>

                    <svg width="400" height="500" viewBox="-200 -50 400 500" className="bg-slate-100 rounded">
                        {/* Airgap */}
                        <rect x="-200" y="-30" width="400" height={0.15 * scale} fill="#e2e8f0" />
                        <text x="-190" y="-10" fontSize="12" fill="#64748b">Rotor Airgap (0.15mm)</text>

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
                                <span className="text-slate-600">Total Conductors ($Z_{slot}$):</span>
                                <span className="font-bold">{metrics.zSlot} wires</span>
                            </div>
                            <div className="flex justify-between border-b pb-1">
                                <span className="text-slate-600">Copper Fill Factor ($\eta_{cu}$):</span>
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