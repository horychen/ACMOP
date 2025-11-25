'use client';

import React from 'react';
import { ResponsiveContainer, RadarChart, PolarGrid, PolarAngleAxis, PolarRadiusAxis, Radar, Tooltip } from 'recharts';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import { FEA_Performance } from '@/lib/DesignData';

interface PerformanceRadarProps {
    data: FEA_Performance;
}

export default function PerformanceRadar({ data }: PerformanceRadarProps) {
    // Prepare data for Radar Chart
    // We need to normalize or select comparable metrics.
    // For now, let's plot the 3 objectives (f1, f2, f3) and maybe a few others if they fit.
    // Note: f2 is usually negative efficiency, so we might want to invert it for visualization or show absolute.

    const chartData = [
        { subject: 'Cost (f1)', A: data.f1, fullMark: 100 },
        { subject: 'Efficiency (f2)', A: Math.abs(data.f2) * 100, fullMark: 100 }, // Scaled for visibility
        { subject: 'Torque Ripple (f3)', A: data.f3, fullMark: 100 },
        { subject: 'Torque Avg', A: data.torque_average, fullMark: 100 },
        { subject: 'Force Err', A: data.normalized_force_error_magnitude * 100, fullMark: 100 },
    ];

    return (
        <Card className="h-full flex flex-col" style={{ minHeight: '400px' }}>
            <CardHeader className="pb-2 flex-shrink-0">
                <CardTitle className="text-lg">性能概览</CardTitle>
            </CardHeader>
            <CardContent className="flex-1" style={{ minHeight: '300px', position: 'relative' }}>
                <div className="w-full h-full" style={{ minHeight: '300px', minWidth: '0' }}>
                    <ResponsiveContainer width="100%" height="100%" minHeight={300}>
                        <RadarChart cx="50%" cy="50%" outerRadius="80%" data={chartData}>
                            <PolarGrid />
                            <PolarAngleAxis dataKey="subject" tick={{ fontSize: 12 }} />
                            <PolarRadiusAxis angle={30} domain={[0, 'auto']} tick={{ fontSize: 10 }} />
                            <Radar
                                name="设计"
                                dataKey="A"
                                stroke="#8884d8"
                                fill="#8884d8"
                                fillOpacity={0.6}
                            />
                            <Tooltip />
                        </RadarChart>
                    </ResponsiveContainer>
                </div>
            </CardContent>
        </Card>
    );
}
