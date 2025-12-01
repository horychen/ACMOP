import React, { useMemo } from 'react';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';

interface WindingDiagramsProps {
    Qs: number;
    p: number;
    m: number;
    layer_X_phases: string[];
    layer_X_signs: string[];
}

// Helper to convert polar to cartesian
const polarToCartesian = (radius: number, angleDegrees: number) => {
    const angleRadians = (angleDegrees * Math.PI) / 180;
    return {
        x: radius * Math.cos(angleRadians),
        y: radius * Math.sin(angleRadians)
    };
};

export default function WindingDiagrams({ Qs, p, m, layer_X_phases, layer_X_signs }: WindingDiagramsProps) {
    const RADIUS = 120;
    const CENTER = { x: 150, y: 150 };
    const SVG_SIZE = 300;

    // Calculate electrical angle for each slot
    const slotData = useMemo(() => {
        const anglePerSlot = (360 * p) / Qs;
        return Array.from({ length: Qs }).map((_, i) => {
            const mechAngle = (i * 360) / Qs;
            const elecAngle = i * anglePerSlot;
            // Normalize to 0-360 for some calculations, but keep cumulative for spiral
            const normalizedElecAngle = elecAngle % 360;

            // Determine radius based on spiral (wrapping)
            // Backend logic: RADIUS - (PHI // 360) * 1.0
            // We scale it up for pixels. Say 10px per wrap.
            const wrapCount = Math.floor(elecAngle / 360);
            const radius = RADIUS - wrapCount * 10;

            return {
                id: i + 1,
                mechAngle,
                elecAngle,
                normalizedElecAngle,
                radius,
                phase: layer_X_phases[i],
                sign: layer_X_signs[i]
            };
        });
    }, [Qs, p, layer_X_phases, layer_X_signs]);

    // Phase Colors (matching WindingLayoutViewer)
    const colors: Record<string, string> = {
        'U': '#f43f5e', // rose-500
        'V': '#06b6d4', // cyan-500
        'W': '#6366f1', // indigo-500
        '': '#e5e7eb'   // gray-200
    };

    const renderStarOfSlots = () => {
        return (
            <svg width={SVG_SIZE} height={SVG_SIZE} viewBox={`0 0 ${SVG_SIZE} ${SVG_SIZE}`}>
                {/* Background Sectors for Phase Belts (Assuming 60 degree belts for m=3) */}
                {/* This is a simplification; ideally we calculate belts exactly like backend. 
                    For m=3, belts are 60 deg. 
                    +U: -30 to 30
                    -W: 30 to 90
                    +V: 90 to 150
                    -U: 150 to 210
                    +W: 210 to 270
                    -V: 270 to 330
                */}
                <g transform={`translate(${CENTER.x}, ${CENTER.y})`}>
                    {/* Axes/Grid */}
                    <circle cx={0} cy={0} r={RADIUS + 20} fill="none" stroke="#e2e8f0" />

                    {/* Phase Sectors Labels */}
                    {/* Using simple text for now at approximate locations */}
                    <text x={RADIUS + 10} y={0} textAnchor="start" dominantBaseline="middle" className="text-[10px] fill-slate-400 font-mono">+U</text>
                    <text x={0} y={-(RADIUS + 10)} textAnchor="middle" dominantBaseline="auto" className="text-[10px] fill-slate-400 font-mono">-V</text>
                    <text x={-(RADIUS + 10)} y={0} textAnchor="end" dominantBaseline="middle" className="text-[10px] fill-slate-400 font-mono">-U</text>
                    <text x={0} y={RADIUS + 10} textAnchor="middle" dominantBaseline="hanging" className="text-[10px] fill-slate-400 font-mono">+V</text>

                    {/* Slot Phasors */}
                    {slotData.map((slot) => {
                        const { x, y } = polarToCartesian(slot.radius, slot.elecAngle);
                        const color = colors[slot.phase] || '#94a3b8';

                        return (
                            <g key={slot.id}>
                                {/* Arrow Line */}
                                <line
                                    x1={0} y1={0} x2={x} y2={y}
                                    stroke={color}
                                    strokeWidth="1.5"
                                    opacity="0.6"
                                />
                                {/* Arrow Head (Circle for now) */}
                                <circle cx={x} cy={y} r={2} fill={color} />

                                {/* Label */}
                                {/* Push label out a bit */}
                                {(() => {
                                    const labelPos = polarToCartesian(slot.radius + 15, slot.elecAngle);
                                    return (
                                        <text
                                            x={labelPos.x}
                                            y={labelPos.y}
                                            textAnchor="middle"
                                            dominantBaseline="middle"
                                            className="text-[8px] font-bold fill-slate-700"
                                        >
                                            {slot.id}
                                        </text>
                                    );
                                })()}
                            </g>
                        );
                    })}
                </g>
            </svg>
        );
    };

    const renderConnectionStar = () => {
        // Group slots by phase and align them
        // For Connection Star, we rotate phasors so that:
        // U is at 0
        // V is at 120 (or -120 depending on convention)
        // W is at 240
        // And negative phases are flipped 180

        // In backend:
        // if key in 'abc' (negative): phase_shift = 180
        // if key in 'ABC' (positive): phase_shift = 0

        // We can just use the slot's assigned phase to determine the target sector
        // But to replicate the "Star", we usually plot the *actual* electrical phasors 
        // but grouped/colored. 
        // Actually, the Connection Star in the paper/backend often shows the phasors *after* 
        // being referred to the fundamental phase axis? 
        // Let's look at backend: `PHI, label = phase_shift+PHI_ori`
        // It shifts negative belts by 180 to align with positive.

        return (
            <svg width={SVG_SIZE} height={SVG_SIZE} viewBox={`0 0 ${SVG_SIZE} ${SVG_SIZE}`}>
                <g transform={`translate(${CENTER.x}, ${CENTER.y})`}>
                    <circle cx={0} cy={0} r={RADIUS + 20} fill="none" stroke="#e2e8f0" />

                    {/* Main Phase Axes */}
                    <line x1={0} y1={0} x2={RADIUS} y2={0} stroke={colors['U']} strokeWidth="1" strokeDasharray="4 4" />
                    <text x={RADIUS + 5} y={0} className="text-xs fill-rose-500 font-bold">U</text>

                    {/* V at 120? Backend says: U=0, W=120, V=240 (Note V/W transposed in comments) */}
                    {/* Let's stick to standard U=0, V=-120(240), W=120 for now or follow backend comments */}
                    {/* Backend: U at 0. W at 120. V at 240. */}
                    {(() => {
                        const wPos = polarToCartesian(RADIUS, 120);
                        const vPos = polarToCartesian(RADIUS, 240);
                        return (
                            <>
                                <line x1={0} y1={0} x2={wPos.x} y2={wPos.y} stroke={colors['W']} strokeWidth="1" strokeDasharray="4 4" />
                                <text x={wPos.x + 5} y={wPos.y} className="text-xs fill-indigo-500 font-bold">W</text>

                                <line x1={0} y1={0} x2={vPos.x} y2={vPos.y} stroke={colors['V']} strokeWidth="1" strokeDasharray="4 4" />
                                <text x={vPos.x + 5} y={vPos.y} className="text-xs fill-cyan-500 font-bold">V</text>
                            </>
                        );
                    })()}

                    {slotData.map((slot) => {
                        if (!slot.phase) return null;

                        // Determine shift based on phase
                        let shift = 0;
                        let isNegative = slot.sign === '-';

                        // We want to align everything to the "Positive" axis of its phase
                        // If it's U-, we add 180 to bring it to U+ (or vice versa? Backend says phase_shift=180 for negative)
                        // Wait, if it's U- (e.g. at 180), adding 180 brings it to 360 (0). Yes.

                        // However, we also want to visualize them *clustered*.
                        // The backend `draw_connection_star` iterates by phase belt.
                        // It shifts the angle so they all point roughly in the same direction (the resultant MMF direction).

                        // Let's just plot the phasors as they are, but maybe color coded?
                        // No, the connection star specifically shows how they add up.
                        // So we should apply the 180 shift for negative signs.

                        let displayAngle = slot.elecAngle;
                        if (isNegative) {
                            displayAngle += 180;
                        }

                        const { x, y } = polarToCartesian(slot.radius, displayAngle);
                        const color = colors[slot.phase] || '#94a3b8';

                        return (
                            <g key={`conn-${slot.id}`}>
                                <line
                                    x1={0} y1={0} x2={x} y2={y}
                                    stroke={color}
                                    strokeWidth="1.5"
                                    opacity="0.6"
                                />
                                <circle cx={x} cy={y} r={2} fill={color} />
                                {(() => {
                                    const labelPos = polarToCartesian(slot.radius + 15, displayAngle);
                                    return (
                                        <text
                                            x={labelPos.x}
                                            y={labelPos.y}
                                            textAnchor="middle"
                                            dominantBaseline="middle"
                                            className="text-[8px] font-bold fill-slate-700"
                                        >
                                            {isNegative ? `-${slot.id}` : slot.id}
                                        </text>
                                    );
                                })()}
                            </g>
                        );
                    })}
                </g>
            </svg>
        );
    }

    return (
        <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            <Card>
                <CardHeader className="pb-2">
                    <CardTitle className="text-sm font-medium">Star of Slots (Phasor Diagram)</CardTitle>
                </CardHeader>
                <CardContent className="flex justify-center p-4">
                    {renderStarOfSlots()}
                </CardContent>
            </Card>
            <Card>
                <CardHeader className="pb-2">
                    <CardTitle className="text-sm font-medium">Connection Star (Mmf Phasors)</CardTitle>
                </CardHeader>
                <CardContent className="flex justify-center p-4">
                    {renderConnectionStar()}
                </CardContent>
            </Card>
        </div>
    );
}
