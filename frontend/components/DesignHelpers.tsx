"use client";

import React from 'react';
import { useTheme } from '@/context/ThemeContext';

export const InputGroup = ({ label, value, onChange, step = 1, icon }: { label: string, value: number, onChange: (v: string) => void, step?: number, icon?: React.ReactNode }) => {
    const { theme } = useTheme();
    
    return (
        <div className="relative group">
            <label className={`block text-[10px] font-medium mb-1 uppercase ${
                theme === 'dark' ? 'text-slate-500' : 'text-slate-600'
            }`}>{label}</label>
            <div className="relative">
                <input
                    type="number"
                    step={step}
                    value={value}
                    onChange={(e) => onChange(e.target.value)}
                    className={`w-full rounded px-3 py-2 text-sm font-mono outline-none transition-all border focus:border-brand-500 focus:ring-1 focus:ring-brand-500 ${
                        theme === 'dark' 
                            ? 'bg-slate-950 border-slate-700 text-slate-200' 
                            : 'bg-white border-slate-300 text-slate-900'
                    }`}
                />
                {icon && (
                    <div className={`absolute right-3 top-2.5 pointer-events-none opacity-50 ${
                        theme === 'dark' ? 'text-slate-600' : 'text-slate-500'
                    }`}>{icon}</div>
                )}
            </div>
        </div>
    );
};

export const KpiCard = ({ label, value, unit, icon }: { label: string, value: string | undefined, unit: string, icon: React.ReactNode }) => {
    const { theme } = useTheme();
    
    return (
        <div className={`p-4 rounded-lg border flex flex-col justify-between ${
            theme === 'dark' 
                ? 'bg-slate-800/50 border-slate-700' 
                : 'bg-slate-50 border-slate-200'
        }`}>
            <div className="flex items-start justify-between mb-2">
                <span className={`text-xs font-medium ${
                    theme === 'dark' ? 'text-slate-400' : 'text-slate-600'
                }`}>{label}</span>
                {icon}
            </div>
            <div className="flex items-baseline">
                <span className={`text-2xl font-bold tracking-tight mr-1 ${
                    theme === 'dark' ? 'text-white' : 'text-slate-900'
                }`}>{value || '-'}</span>
                <span className={`text-xs font-mono ${
                    theme === 'dark' ? 'text-slate-500' : 'text-slate-600'
                }`}>{unit}</span>
            </div>
        </div>
    );
};

export const ResultRow = ({ label, value, unit }: { label: string, value: string | undefined, unit: string }) => {
    const { theme } = useTheme();
    
    return (
        <tr className={`transition-colors ${
            theme === 'dark' ? 'hover:bg-slate-800/50' : 'hover:bg-slate-100/50'
        }`}>
            <td className={`px-6 py-3 font-medium ${
                theme === 'dark' ? 'text-slate-300' : 'text-slate-700'
            }`}>{label}</td>
            <td className={`px-6 py-3 font-mono ${
                theme === 'dark' ? 'text-brand-100' : 'text-brand-600'
            }`}>{value}</td>
            <td className={`px-6 py-3 ${
                theme === 'dark' ? 'text-slate-500' : 'text-slate-600'
            }`}>{unit}</td>
        </tr>
    );
};
