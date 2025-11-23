import React from 'react';

export const InputGroup = ({ label, value, onChange, step = 1, icon }: { label: string, value: number, onChange: (v: string) => void, step?: number, icon?: React.ReactNode }) => (
    <div className="relative group">
        <label className="block text-[10px] text-slate-500 font-medium mb-1 uppercase">{label}</label>
        <div className="relative">
            <input
                type="number"
                step={step}
                value={value}
                onChange={(e) => onChange(e.target.value)}
                className="w-full bg-slate-950 border border-slate-700 rounded px-3 py-2 text-sm text-slate-200 focus:border-brand-500 focus:ring-1 focus:ring-brand-500 outline-none transition-all font-mono"
            />
            {icon && <div className="absolute right-3 top-2.5 text-slate-600 pointer-events-none opacity-50">{icon}</div>}
        </div>
    </div>
);

export const KpiCard = ({ label, value, unit, icon }: { label: string, value: string | undefined, unit: string, icon: React.ReactNode }) => (
    <div className="bg-slate-800/50 p-4 rounded-lg border border-slate-700 flex flex-col justify-between">
        <div className="flex items-start justify-between mb-2">
            <span className="text-xs text-slate-400 font-medium">{label}</span>
            {icon}
        </div>
        <div className="flex items-baseline">
            <span className="text-2xl font-bold text-white tracking-tight mr-1">{value || '-'}</span>
            <span className="text-xs text-slate-500 font-mono">{unit}</span>
        </div>
    </div>
);

export const ResultRow = ({ label, value, unit }: { label: string, value: string | undefined, unit: string }) => (
    <tr className="hover:bg-slate-800/50 transition-colors">
        <td className="px-6 py-3 font-medium text-slate-300">{label}</td>
        <td className="px-6 py-3 font-mono text-brand-100">{value}</td>
        <td className="px-6 py-3 text-slate-500">{unit}</td>
    </tr>
);
