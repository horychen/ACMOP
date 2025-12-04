import React, { useMemo } from 'react';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';

interface WindingDiagramsProps {
    Qs: number;
    p: number;
    ps?: number; // Suspension pole pairs
    m: number;
    layer_X_phases: (string | null)[];
    layer_X_signs: (string | null)[];
}

// Helper to convert polar to cartesian
const polarToCartesian = (radius: number, angleDegrees: number) => {
    const angleRadians = (angleDegrees * Math.PI) / 180;
    return {
        x: radius * Math.cos(angleRadians),
        y: radius * Math.sin(angleRadians)
    };
};

export default function WindingDiagrams({ Qs, p, ps, m, layer_X_phases, layer_X_signs }: WindingDiagramsProps) {
    const RADIUS = 120;
    const CENTER = { x: 150, y: 150 };
    const SVG_SIZE = 300;

    // Phase Colors
    const colors: Record<string, string> = {
        'U': '#f43f5e', // rose-500
        'V': '#06b6d4', // cyan-500
        'W': '#6366f1', // indigo-500
        '': '#e5e7eb'   // gray-200
    };

    const getPhaseColor = (phase: string | null) => {
        if (!phase) return colors[''];
        return colors[phase] || colors[''];
    };

    // Calculate slot data for Torque (p)
    const torqueSlotData = useMemo(() => {
        const slots = [];
        const electricalAnglePerSlot = (360 * p) / Qs;
        const radiusBase = 80;
        const spiralFactor = 1.5;

        for (let i = 0; i < Qs; i++) {
            const slotIndex = i + 1;
            const electricalAngle = (i * electricalAnglePerSlot) % 360;

            const wrapCount = Math.floor((i * electricalAnglePerSlot) / 360);
            const radius = radiusBase - (wrapCount * spiralFactor);

            const pos = polarToCartesian(radius, electricalAngle);

            let phase = 'A';
            let sign = '+';
            if (layer_X_phases && layer_X_phases[i]) {
                phase = layer_X_phases[i] as string;
            }
            if (layer_X_signs && layer_X_signs[i]) {
                sign = layer_X_signs[i] as string;
            }

            slots.push({
                id: slotIndex,
                angle: electricalAngle,
                radius,
                x: pos.x,
                y: pos.y,
                phase,
                sign
            });
        }
        return slots;
    }, [Qs, p, layer_X_phases, layer_X_signs]);

    // Calculate slot data for Suspension (ps)
    const suspensionSlotData = useMemo(() => {
        if (!ps) return [];
        const slots = [];
        const electricalAnglePerSlot = (360 * ps) / Qs;
        const radiusBase = 80;
        const spiralFactor = 1.5;

        for (let i = 0; i < Qs; i++) {
            const slotIndex = i + 1;
            const electricalAngle = (i * electricalAnglePerSlot) % 360;

            const wrapCount = Math.floor((i * electricalAnglePerSlot) / 360);
            const radius = radiusBase - (wrapCount * spiralFactor);

            const pos = polarToCartesian(radius, electricalAngle);

            let phase = 'A';
            if (layer_X_phases && layer_X_phases[i]) {
                phase = layer_X_phases[i] as string;
            }

            slots.push({
                id: slotIndex,
                angle: electricalAngle,
                radius,
                x: pos.x,
                y: pos.y,
                phase,
                sign: '+'
            });
        }
        return slots;
    }, [Qs, ps, layer_X_phases]);

    const renderStarOfSlots = (data: typeof torqueSlotData, title: string) => (
        <div className="flex flex-col items-center">
            <h4 className="text-sm font-medium mb-2">{title}</h4>
            <svg width="300" height="300" viewBox="-120 -120 240 240" className="border rounded bg-white">
                {/* Grid circles */}
                <circle cx="0" cy="0" r="80" fill="none" stroke="#e5e7eb" strokeWidth="1" />
                <circle cx="0" cy="0" r="60" fill="none" stroke="#e5e7eb" strokeWidth="1" strokeDasharray="4 4" />

                {/* Axes */}
                <line x1="-100" y1="0" x2="100" y2="0" stroke="#e5e7eb" strokeWidth="1" />
                <line x1="0" y1="-100" x2="0" y2="100" stroke="#e5e7eb" strokeWidth="1" />

                {/* Slots */}
                {data.map((slot) => (
                    <g key={slot.id}>
                        <line
                            x1="0" y1="0"
                            x2={slot.x} y2={slot.y}
                            stroke={getPhaseColor(slot.phase)}
                            strokeWidth="1.5"
                            markerEnd="url(#arrowhead)"
                        />
                        <text
                            x={slot.x * 1.15}
                            y={slot.y * 1.15}
                            textAnchor="middle"
                            dominantBaseline="middle"
                            fontSize="10"
                            fill="#374151"
                        >
                            {slot.id}
                        </text>
                    </g>
                ))}

                {/* Arrow Marker Definition */}
                <defs>
                    <marker id="arrowhead" markerWidth="10" markerHeight="7" refX="9" refY="3.5" orient="auto">
                        <polygon points="0 0, 10 3.5, 0 7" fill="#9ca3af" />
                    </marker>
                </defs>
            </svg>
        </div>
    );

    const renderConnectionStar = () => (
        <div className="flex flex-col items-center">
            <h4 className="text-sm font-medium mb-2">Torque MMF (Connection Star)</h4>
            <svg width="300" height="300" viewBox="-120 -120 240 240" className="border rounded bg-white">
                {/* Grid circles */}
                <circle cx="0" cy="0" r="80" fill="none" stroke="#e5e7eb" strokeWidth="1" />

                {/* Phase Sectors (Simplified visualization) */}
                {/* Ideally we should draw sectors based on m and phase_belt */}

                {/* Slots with phase shift for negative signs */}
                {torqueSlotData.map((slot) => {
                    let angle = slot.angle;
                    let label = `${slot.sign === '-' ? '-' : '+'}${slot.id}`;

                    // Apply 180 degree shift for negative connections to align MMF
                    if (slot.sign === '-') {
                        angle = (angle + 180) % 360;
                    }

                    const pos = polarToCartesian(slot.radius, angle);

                    return (
                        <g key={slot.id}>
                            <line
                                x1="0" y1="0"
                                x2={pos.x} y2={pos.y}
                                stroke={getPhaseColor(slot.phase)}
                                strokeWidth="1.5"
                                markerEnd="url(#arrowhead)"
                            />
                            <text
                                x={pos.x * 1.15}
                                y={pos.y * 1.15}
                                textAnchor="middle"
                                dominantBaseline="middle"
                                fontSize="10"
                                fill="#374151"
                                fontWeight="bold"
                            >
                                {label}
                            </text>
                        </g>
                    );
                })}
            </svg>
        </div>
    );

    return (
        <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            <Card>
                <CardContent className="pt-6">
                    {renderStarOfSlots(torqueSlotData, `Torque Star of Slots (p=${p})`)}
                </CardContent>
            </Card>

            {ps && (
                <Card>
                    <CardContent className="pt-6">
                        {renderStarOfSlots(suspensionSlotData, `Suspension Star of Slots (ps=${ps})`)}
                    </CardContent>
                </Card>
            )}

            <Card>
                <CardContent className="pt-6">
                    {renderConnectionStar()}
                </CardContent>
            </Card>
        </div>
    );
}
