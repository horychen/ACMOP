"use client";

import React from 'react';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, AreaChart, Area } from 'recharts';
import { useTheme } from '@/context/ThemeContext';

interface DataPoint {
    speed: number;
    efficiency: number;
}

interface EfficiencyChartProps {
    data: DataPoint[];
}

export const EfficiencyChart: React.FC<EfficiencyChartProps> = ({ data }) => {
    const { theme } = useTheme();
    
    const isDark = theme === 'dark';
    
    return (
        <div className={`h-64 w-full rounded-lg p-4 border ${
            isDark 
                ? 'bg-slate-800 border-slate-700' 
                : 'bg-white border-slate-200'
        }`}>
            <h3 className={`text-sm font-semibold mb-4 ${
                isDark ? 'text-slate-300' : 'text-slate-700'
            }`}>Efficiency Map (%) vs Speed (RPM)</h3>
            <ResponsiveContainer width="100%" height="100%">
                <AreaChart data={data}>
                    <defs>
                        <linearGradient id="colorEff" x1="0" y1="0" x2="0" y2="1">
                            <stop offset="5%" stopColor="#0ea5e9" stopOpacity={isDark ? 0.3 : 0.2} />
                            <stop offset="95%" stopColor="#0ea5e9" stopOpacity={0} />
                        </linearGradient>
                    </defs>
                    <CartesianGrid 
                        strokeDasharray="3 3" 
                        stroke={isDark ? "#334155" : "#cbd5e1"} 
                    />
                    <XAxis
                        dataKey="speed"
                        stroke={isDark ? "#94a3b8" : "#64748b"}
                        tick={{ fontSize: 12, fill: isDark ? "#94a3b8" : "#64748b" }}
                        label={{ 
                            value: 'RPM', 
                            position: 'insideBottomRight', 
                            offset: -5, 
                            fill: isDark ? "#94a3b8" : "#64748b" 
                        }}
                    />
                    <YAxis
                        stroke={isDark ? "#94a3b8" : "#64748b"}
                        tick={{ fontSize: 12, fill: isDark ? "#94a3b8" : "#64748b" }}
                        domain={[0, 100]}
                    />
                    <Tooltip
                        contentStyle={{ 
                            backgroundColor: isDark ? '#1e293b' : '#f8fafc', 
                            borderColor: isDark ? '#475569' : '#cbd5e1', 
                            color: isDark ? '#f1f5f9' : '#1e293b' 
                        }}
                        itemStyle={{ color: isDark ? '#38bdf8' : '#0284c7' }}
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
