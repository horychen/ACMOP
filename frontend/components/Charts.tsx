"use client";

import React from 'react';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, AreaChart, Area, PieChart, Pie, Cell, Legend } from 'recharts';
import { useTheme } from '@/context/ThemeContext';

interface DataPoint {
    time: number;
    value: number;
}

interface EfficiencyChartProps {
    data: DataPoint[];
    title?: string;
}

export const EfficiencyChart: React.FC<EfficiencyChartProps> = ({ data, title }) => {
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
            }`}>{title || 'Chart'}</h3>
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
                        stroke={isDark ? "#334155" : "#e5e7eb"} 
                    />
                    <XAxis
                        dataKey="time"
                        stroke={isDark ? "#94a3b8" : "#6b7280"}
                        tick={{ fontSize: 12, fill: isDark ? "#94a3b8" : "#6b7280" }}
                        tickFormatter={(value) => value.toFixed(4)}
                        label={{ 
                            value: 'Time (s)', 
                            position: 'insideBottomRight', 
                            offset: -5, 
                            fill: isDark ? "#94a3b8" : "#6b7280" 
                        }}
                    />
                    <YAxis
                        stroke={isDark ? "#94a3b8" : "#6b7280"}
                        tick={{ fontSize: 12, fill: isDark ? "#94a3b8" : "#6b7280" }}
                        tickFormatter={(value) => value.toFixed(2)}
                        domain={['auto', 'auto']}
                    />
                    <Tooltip
                        contentStyle={{ 
                            backgroundColor: isDark ? '#1e293b' : '#ffffff', 
                            borderColor: isDark ? '#475569' : '#d1d5db', 
                            color: isDark ? '#f1f5f9' : '#111827',
                            borderRadius: '6px'
                        }}
                        itemStyle={{ color: isDark ? '#38bdf8' : '#0284c7' }}
                    />
                    <Area
                        type="monotone"
                        dataKey="value"
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

interface LossDataPoint {
    name: string;
    value: number;
}

interface DonutChartProps {
    data: LossDataPoint[];
    title?: string;
}

export const DonutChart: React.FC<DonutChartProps> = ({ data, title }) => {
    const { theme } = useTheme();
    const isDark = theme === 'dark';
    
    // Color palette for different loss types
    const COLORS = [
        '#3b82f6', // Blue - Stator Copper Loss
        '#8b5cf6', // Purple - Stator Copper Loss (End Turn)
        '#ef4444', // Red - Iron Loss
        '#f59e0b', // Amber - Magnet Joule Loss
        '#10b981', // Green - Windage Loss
        '#06b6d4', // Cyan - Other
    ];
    
    // Calculate total for percentage display
    const total = data.reduce((sum, item) => sum + item.value, 0);
    
    // Format tooltip to show both value and percentage
    const CustomTooltip = ({ active, payload }: any) => {
        if (active && payload && payload.length) {
            const data = payload[0];
            const percentage = total > 0 ? ((data.value / total) * 100).toFixed(2) : 0;
            return (
                <div className={`rounded-lg border p-3 ${
                    isDark 
                        ? 'bg-slate-800 border-slate-700 text-slate-200' 
                        : 'bg-white border-slate-200 text-slate-800'
                }`}>
                    <p className="font-semibold">{data.name}</p>
                    <p className="text-sm">
                        {data.value.toFixed(2)} W ({percentage}%)
                    </p>
                </div>
            );
        }
        return null;
    };
    
    return (
        <div className={`h-full w-full rounded-lg p-4 border ${
            isDark 
                ? 'bg-slate-800 border-slate-700' 
                : 'bg-white border-slate-200'
        }`}>
            {title && (
                <h3 className={`text-sm font-semibold mb-4 ${
                    isDark ? 'text-slate-300' : 'text-slate-700'
                }`}>{title}</h3>
            )}
            <ResponsiveContainer width="100%" height="100%">
                <PieChart>
                    <Pie
                        data={data}
                        cx="50%"
                        cy="50%"
                        labelLine={false}
                        label={({ name, percent }) => 
                            percent > 0.05 ? `${(percent * 100).toFixed(0)}%` : ''
                        }
                        outerRadius={80}
                        innerRadius={50}
                        fill="#8884d8"
                        dataKey="value"
                    >
                        {data.map((entry, index) => (
                            <Cell key={`cell-${index}`} fill={COLORS[index % COLORS.length]} />
                        ))}
                    </Pie>
                    <Tooltip content={<CustomTooltip />} />
                    <Legend 
                        verticalAlign="bottom" 
                        height={36}
                        formatter={(value, entry: any) => (
                            <span style={{ 
                                color: isDark ? '#cbd5e1' : '#475569',
                                fontSize: '12px'
                            }}>
                                {value}
                            </span>
                        )}
                    />
                </PieChart>
            </ResponsiveContainer>
        </div>
    );
};
