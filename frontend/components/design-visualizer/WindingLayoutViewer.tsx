'use client';

import React from 'react';
import { WindingLayout } from '@/lib/DesignData';
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from '@/components/ui/card';

interface WindingLayoutViewerProps {
    exUserData: any; // ExUser type
    Qs: number;
    p?: number; // Pole pairs
}

export default function WindingLayoutViewer({ exUserData, Qs, p = 1 }: WindingLayoutViewerProps) {
    const data = exUserData.wily;
    const phases = ['U', 'V', 'W'];
    const colors: Record<string, string> = {
        'U': '#f43f5e', // rose-500 - Modern coral/pink
        'V': '#06b6d4', // cyan-500 - Bright teal
        'W': '#6366f1', // indigo-500 - Deep purple-blue
        '': '#e5e7eb'   // gray-200
    };

    // Calculate SPP (slots per phase per pole)
    const SPP = Qs / (3 * 2 * p);
    const kd1 = data.kd1 || 0;
    const kp1 = data.kp1 || 0;
    const kw1 = kd1 * kp1; // Winding factor
    const DriveW_zQ = exUserData.DriveW_zQ || 0;
    const BeariW_zQ = exUserData.BeariW_zQ || 0;
    const DriveW_CurrentAmp = exUserData.DriveW_CurrentAmp || 0;
    const BeariW_CurrentAmp = exUserData.BeariW_CurrentAmp || 0;
    const grouping_AC = data.grouping_AC || [];

    const coilPitch = data.coil_pitch_y || 1;
    const slotWidth = 40;
    const slotHeight = 100;
    const slotGap = 10;
    const startX = 50;
    const startY = 100; // Space for top connections

    const totalWidth = startX + Qs * (slotWidth + slotGap) + 50;
    const totalHeight = 400; // Increased height for rotor and axes

    // Rotor parameters
    const numPoles = 2 * p;
    const totalLength = Qs * (slotWidth + slotGap);
    const poleWidth = totalLength / numPoles;

    const renderConnections = () => {
        const connections = [];
        // Draw connections from Layer X (Top) to Layer Y (Bottom) with pitch y
        // Assuming Top layer is at index i, connects to Bottom layer at index i + y
        for (let i = 0; i < Qs; i++) {
            const x1 = startX + i * (slotWidth + slotGap) + slotWidth / 2;
            const targetSlot = (i + coilPitch) % Qs;

            const x2 = startX + (i + coilPitch) * (slotWidth + slotGap) + slotWidth / 2;

            const phase = data.layer_X_phases[i];
            const color = colors[phase] || '#ccc';

            const arcHeight = 40 + (i % 2) * 10; // Stagger slightly

            // If it wraps around the linear diagram
            if (i + coilPitch >= Qs) {
                // Draw to end and from start
                const xEnd = startX + Qs * (slotWidth + slotGap);
                const xStart = startX;
                const xTarget = startX + (targetSlot) * (slotWidth + slotGap) + slotWidth / 2;

                // Line to right edge
                connections.push(
                    <path
                        key={`conn-wrap-1-${i}`}
                        d={`M ${x1} ${startY} Q ${x1 + 20} ${startY - arcHeight} ${xEnd + 20} ${startY - arcHeight * 0.5}`}
                        fill="none"
                        stroke={color}
                        strokeWidth="1.5"
                        strokeDasharray="4 2"
                        opacity="0.5"
                    />
                );
                // Line from left edge
                connections.push(
                    <path
                        key={`conn-wrap-2-${i}`}
                        d={`M ${xStart - 20} ${startY - arcHeight * 0.5} Q ${xTarget - 20} ${startY - arcHeight} ${xTarget} ${startY + slotHeight * 0.75}`} // Connect to center of Layer Y
                        fill="none"
                        stroke={color}
                        strokeWidth="1.5"
                        strokeDasharray="4 2"
                        opacity="0.5"
                    />
                );

            } else {
                connections.push(
                    <path
                        key={`conn-${i}`}
                        d={`M ${x1} ${startY} C ${x1} ${startY - arcHeight}, ${x2} ${startY - arcHeight}, ${x2} ${startY + slotHeight * 0.75}`} // Connects to center of Layer Y
                        fill="none"
                        stroke={color}
                        strokeWidth="2"
                        opacity="0.7"
                    />
                );
            }
        }

        // Add a label for coil pitch on the first valid connection
        if (Qs > coilPitch) {
            const i = 0;
            const x1 = startX + i * (slotWidth + slotGap) + slotWidth / 2;
            const x2 = startX + (i + coilPitch) * (slotWidth + slotGap) + slotWidth / 2;
            const midX = (x1 + x2) / 2;
            const arcHeight = 60;

            connections.push(
                <g key="pitch-label">
                    <path
                        d={`M ${x1} ${startY - 10} L ${x1} ${startY - arcHeight - 10} M ${x2} ${startY - 10} L ${x2} ${startY - arcHeight - 10} M ${x1} ${startY - arcHeight} L ${x2} ${startY - arcHeight}`}
                        stroke="black"
                        strokeWidth="1"
                        markerEnd="url(#arrow)"
                        markerStart="url(#arrow)"
                    />
                    <text x={midX} y={startY - arcHeight - 5} textAnchor="middle" className="text-xs font-bold">
                        y = {coilPitch}
                    </text>
                </g>
            );
        }

        return connections;
    };

    return (
        <div className="space-y-6">
            <Card>
                <CardHeader>
                    <CardTitle>Winding Visualization</CardTitle>
                    <CardDescription>
                        Stator slots (Layer X/Y), Rotor poles, and Coil connections.
                    </CardDescription>
                </CardHeader>
                <CardContent className="overflow-x-auto">
                    <svg width={totalWidth} height={totalHeight} className="min-w-full">
                        <defs>
                            <marker id="arrow" markerWidth="10" markerHeight="10" refX="9" refY="3" orient="auto" markerUnits="strokeWidth">
                                <path d="M0,0 L0,6 L9,3 z" fill="black" />
                            </marker>
                        </defs>

                        {/* Stator Slots */}
                        {Array.from({ length: Qs }).map((_, i) => {
                            const x = startX + i * (slotWidth + slotGap);
                            const phaseX = data.layer_X_phases[i];
                            const signX = data.layer_X_signs[i];
                            const phaseY = data.layer_Y_phases[i];
                            const signY = data.layer_Y_signs[i];
                            const groupAC = grouping_AC[i];

                            return (
                                <g key={`slot-${i}`}>
                                    {/* Slot Number */}
                                    <text x={x + slotWidth / 2} y={startY - 5} textAnchor="middle" className="text-xs fill-muted-foreground">
                                        {i + 1}
                                    </text>

                                    {/* Slot Outline (Teeth) */}
                                    <rect
                                        x={x}
                                        y={startY}
                                        width={slotWidth}
                                        height={slotHeight}
                                        fill="none"
                                        stroke="currentColor"
                                        strokeWidth="1"
                                        className="text-slate-300"
                                    />

                                    {/* Layer X (Top) */}
                                    <rect
                                        x={x + 2}
                                        y={startY + 2}
                                        width={slotWidth - 4}
                                        height={slotHeight / 2 - 4}
                                        fill={colors[phaseX] || '#ccc'}
                                        rx="2"
                                    />
                                    <text
                                        x={x + slotWidth / 2}
                                        y={startY + slotHeight / 4 - 5}
                                        textAnchor="middle"
                                        dominantBaseline="middle"
                                        className={`text-xs font-bold ${phaseX ? 'fill-white' : 'fill-slate-700'}`}
                                    >
                                        {phaseX}{signX}
                                    </text>
                                    {/* Group AC indicator for Layer X */}
                                    {groupAC !== undefined && (
                                        <text
                                            x={x + slotWidth / 2}
                                            y={startY + slotHeight / 4 + 8}
                                            textAnchor="middle"
                                            dominantBaseline="middle"
                                            className={`text-[10px] font-semibold ${phaseX ? 'fill-white' : 'fill-slate-700'}`}
                                        >
                                            {groupAC === 0 ? 'TI' : 'SI'}
                                        </text>
                                    )}

                                    {/* Layer Y (Bottom) */}
                                    <rect
                                        x={x + 2}
                                        y={startY + slotHeight / 2 + 2}
                                        width={slotWidth - 4}
                                        height={slotHeight / 2 - 4}
                                        fill={colors[phaseY] || '#ccc'}
                                        rx="2"
                                    />
                                    <text
                                        x={x + slotWidth / 2}
                                        y={startY + 3 * slotHeight / 4}
                                        textAnchor="middle"
                                        dominantBaseline="middle"
                                        className={`text-xs font-bold ${phaseY ? 'fill-white' : 'fill-slate-700'}`}
                                    >
                                        {phaseY}{signY}
                                    </text>
                                </g>
                            );
                        })}

                        {/* Coil Connections */}
                        {renderConnections()}

                        {/* Rotor Poles */}
                        <g transform={`translate(${startX}, ${startY + slotHeight + 20})`}>
                            {Array.from({ length: numPoles }).map((_, i) => {
                                const x = i * poleWidth;
                                const isN = i % 2 === 0;
                                return (
                                    <g key={`pole-${i}`}>
                                        <rect
                                            x={x}
                                            y={0}
                                            width={poleWidth}
                                            height={40}
                                            fill={isN ? '#ef4444' : '#3b82f6'} // Red for N, Blue for S
                                            opacity="0.8"
                                            stroke="white"
                                            strokeWidth="1"
                                        />
                                        <text x={x + poleWidth / 2} y={25} textAnchor="middle" className="text-sm font-bold fill-white">
                                            {isN ? 'N' : 'S'}
                                        </text>
                                    </g>
                                );
                            })}
                            <text x={-10} y={25} textAnchor="end" className="text-xs font-semibold fill-muted-foreground">Rotor</text>
                        </g>

                        <text x={startX - 10} y={startY + slotHeight / 4} textAnchor="end" className="text-xs font-semibold fill-muted-foreground">Layer X</text>
                        <text x={startX - 10} y={startY + 3 * slotHeight / 4} textAnchor="end" className="text-xs font-semibold fill-muted-foreground">Layer Y</text>

                        {/* Angular Axes */}
                        <g transform={`translate(${startX}, ${startY + slotHeight + 80})`}>
                            {/* Axis Line */}
                            <line x1={0} y1={0} x2={totalLength} y2={0} stroke="currentColor" strokeWidth="1" className="text-slate-400" />

                            {/* Ticks and Labels */}
                            {Array.from({ length: Qs + 1 }).map((_, i) => {
                                const x = i * (slotWidth + slotGap);
                                // Mechanical Angle
                                const mechAngle = (i * 360 / Qs).toFixed(0);
                                // Electrical Angle
                                const elecAngle = (i * 360 * p / Qs).toFixed(0);

                                return (
                                    <g key={`axis-${i}`}>
                                        <line x1={x} y1={0} x2={x} y2={5} stroke="currentColor" strokeWidth="1" className="text-slate-400" />

                                        {/* Mechanical Angle Label */}
                                        <text x={x} y={20} textAnchor="middle" className="text-[10px] fill-muted-foreground">
                                            {mechAngle}°
                                        </text>

                                        {/* Electrical Angle Label */}
                                        <text x={x} y={35} textAnchor="middle" className="text-[10px] fill-muted-foreground font-medium">
                                            {elecAngle}°
                                        </text>
                                    </g>
                                );
                            })}

                            {/* Axis Titles */}
                            <text x={-10} y={20} textAnchor="end" className="text-xs font-semibold fill-muted-foreground">
                                Mech. Θ
                            </text>
                            <text x={-10} y={35} textAnchor="end" className="text-xs font-semibold fill-muted-foreground">
                                Elec. α
                            </text>
                        </g>

                    </svg>
                </CardContent>
            </Card>

            {/* Winding Parameters */}
            <div className="grid grid-cols-2 md:grid-cols-4 gap-4 text-sm">
                <div className="bg-slate-100 dark:bg-slate-800 p-3 rounded-lg">
                    <div className="text-muted-foreground text-xs">SPP (slots/phase/pole)</div>
                    <div className="text-lg font-semibold">{SPP.toFixed(2)}</div>
                </div>
                <div className="bg-slate-100 dark:bg-slate-800 p-3 rounded-lg">
                    <div className="text-muted-foreground text-xs">Distribution Factor kd1</div>
                    <div className="text-lg font-semibold">{kd1.toFixed(3)}</div>
                </div>
                <div className="bg-slate-100 dark:bg-slate-800 p-3 rounded-lg">
                    <div className="text-muted-foreground text-xs">Pitch Factor kp1</div>
                    <div className="text-lg font-semibold">{kp1.toFixed(3)}</div>
                </div>
                <div className="bg-slate-100 dark:bg-slate-800 p-3 rounded-lg">
                    <div className="text-muted-foreground text-xs">Winding Factor kw1</div>
                    <div className="text-lg font-semibold">{kw1.toFixed(3)}</div>
                </div>
                <div className="bg-slate-100 dark:bg-slate-800 p-3 rounded-lg">
                    <div className="text-muted-foreground text-xs">Conductors/Slot (Drive)</div>
                    <div className="text-lg font-semibold">{DriveW_zQ.toFixed(0)}</div>
                </div>
                <div className="bg-slate-100 dark:bg-slate-800 p-3 rounded-lg">
                    <div className="text-muted-foreground text-xs">Conductors/Slot (Bearing)</div>
                    <div className="text-lg font-semibold">{BeariW_zQ.toFixed(0)}</div>
                </div>
                <div className="bg-slate-100 dark:bg-slate-800 p-3 rounded-lg">
                    <div className="text-muted-foreground text-xs">Current Amp (Drive)</div>
                    <div className="text-lg font-semibold">{DriveW_CurrentAmp.toFixed(2)} A</div>
                </div>
                <div className="bg-slate-100 dark:bg-slate-800 p-3 rounded-lg">
                    <div className="text-muted-foreground text-xs">Current Amp (Bearing)</div>
                    <div className="text-lg font-semibold">{BeariW_CurrentAmp.toFixed(2)} A</div>
                </div>
            </div>

            <div className="space-y-2">
                <div className="flex gap-4 text-sm justify-center">
                    <div className="flex items-center gap-2">
                        <div className="w-4 h-4 bg-rose-500 rounded"></div> Phase U
                    </div>
                    <div className="flex items-center gap-2">
                        <div className="w-4 h-4 bg-cyan-500 rounded"></div> Phase V
                    </div>
                    <div className="flex items-center gap-2">
                        <div className="w-4 h-4 bg-indigo-500 rounded"></div> Phase W
                    </div>
                </div>
                {grouping_AC.length > 0 && (
                    <div className="text-xs text-center text-muted-foreground">
                        <span className="font-semibold">Layer X Grouping:</span> TI = Torque Inverter, SI = Suspension Inverter
                    </div>
                )}
            </div>
        </div>
    );
}
