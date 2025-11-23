"use client";

import React from 'react';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, AreaChart, Area } from 'recharts';

interface DataPoint {
    speed: number;
    efficiency: number;
}

interface EfficiencyChartProps {
    data: DataPoint[];
}

export const EfficiencyChart: React.FC<EfficiencyChartProps> = ({ data }) => {
    return (
        <div className="h-64 w-full bg-slate-800 rounded-lg p-4 border border-slate-700">
            <h3 className="text-sm font-semibold text-slate-300 mb-4">Efficiency Map (%) vs Speed (RPM)</h3>
            <ResponsiveContainer width="100%" height="100%">
                <AreaChart data={data}>
                    <defs>
                        <linearGradient id="colorEff" x1="0" y1="0" x2="0" y2="1">
                            <stop offset="5%" stopColor="#0ea5e9" stopOpacity={0.3} />
                            <stop offset="95%" stopColor="#0ea5e9" stopOpacity={0} />
                        </linearGradient>
                    </defs>
                    <CartesianGrid strokeDasharray="3 3" stroke="#334155" />
                    <XAxis
                        dataKey="speed"
                        stroke="#94a3b8"
                        tick={{ fontSize: 12 }}
                        label={{ value: 'RPM', position: 'insideBottomRight', offset: -5, fill: '#94a3b8' }}
                    />
                    <YAxis
                        stroke="#94a3b8"
                        tick={{ fontSize: 12 }}
                        domain={[0, 100]}
                    />
                    <Tooltip
                        contentStyle={{ backgroundColor: '#1e293b', borderColor: '#475569', color: '#f1f5f9' }}
                        itemStyle={{ color: '#38bdf8' }}
                    />
                    <Area
                        type="monotone"
                        dataKey="efficiency"
                        stroke="#0ea5e9"
                        fillOpacity={1}
                        fill="url(#colorEff)"
                        strokeWidth={2}
                    />
                </AreaChart>
            </ResponsiveContainer>
        </div>
    );
};
