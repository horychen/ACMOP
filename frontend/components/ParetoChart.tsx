"use client"

import { ScatterChart, Scatter, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, Label } from 'recharts';

interface ParetoChartProps {
    data: any[];
    xKey: string;
    yKey: string;
    xLabel?: string;
    yLabel?: string;
}

export function ParetoChart({ data, xKey, yKey, xLabel, yLabel }: ParetoChartProps) {
    return (
        <ResponsiveContainer width="100%" height={400}>
            <ScatterChart margin={{ top: 20, right: 20, bottom: 20, left: 20 }}>
                <CartesianGrid />
                <XAxis type="number" dataKey={xKey} name={xLabel || xKey}>
                    <Label value={xLabel || xKey} offset={0} position="insideBottom" />
                </XAxis>
                <YAxis type="number" dataKey={yKey} name={yLabel || yKey}>
                    <Label value={yLabel || yKey} angle={-90} position="insideLeft" />
                </YAxis>
                <Tooltip cursor={{ strokeDasharray: '3 3' }} />
                <Scatter name="Designs" data={data} fill="#8884d8" />
            </ScatterChart>
        </ResponsiveContainer>
    );
}
