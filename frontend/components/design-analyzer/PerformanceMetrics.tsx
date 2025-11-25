'use client';

import React from 'react';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from '@/components/ui/table';
import { FEA_Performance } from '@/lib/DesignData';
import { PieChart, Pie, Cell, ResponsiveContainer, Tooltip, Legend } from 'recharts';

interface PerformanceMetricsProps {
    data: FEA_Performance;
}

export default function PerformanceMetrics({ data }: PerformanceMetricsProps) {
    // Group metrics for better readability
    const generalMetrics = [
        { label: 'Torque Average', value: data.torque_average, unit: 'Nm' },
        { label: 'Torque Ripple (Norm)', value: data.normalized_torque_ripple, unit: '' },
        { label: 'Force Error (Norm)', value: data.normalized_force_error_magnitude, unit: '' },
        { label: 'Force Error Angle', value: data.force_error_angle, unit: 'deg' },
        { label: 'Power Factor', value: data.power_factor, unit: '' },
        { label: 'Efficiency (f2)', value: Math.abs(data.f2) * 100, unit: '%' },
    ];

    const costMetrics = [
        { label: 'Total Cost (f1)', value: data.Cost, unit: '$' },
        { label: 'Iron Cost', value: data.Cost_Fe, unit: '$' },
        { label: 'Copper Cost', value: data.Cost_Cu, unit: '$' },
        { label: 'Magnet Cost', value: data.Cost_PM, unit: '$' },
    ];

    const lossMetrics = [
        { label: 'Total Loss', value: data.rated_total_loss, unit: 'W' },
        { label: 'Iron Loss', value: data.rated_iron_loss, unit: 'W' },
        { label: 'Windage Loss', value: data.rated_windage_loss, unit: 'W' },
        { label: 'Magnet Joule Loss', value: data.rated_magnet_Joule_loss, unit: 'W' },
        { label: 'Stator Cu Loss', value: data.rated_stator_copper_loss_along_stack, unit: 'W' },
    ];

    // Data for Donut Chart
    const lossChartData = [
        { name: 'Iron Loss', value: data.rated_iron_loss },
        { name: 'Windage Loss', value: data.rated_windage_loss },
        { name: 'Magnet Loss', value: data.rated_magnet_Joule_loss },
        { name: 'Stator Cu Loss', value: data.rated_stator_copper_loss_along_stack },
    ].filter(item => item.value > 0);

    const COLORS = ['#0088FE', '#00C49F', '#FFBB28', '#FF8042'];

    const renderTable = (title: string, metrics: { label: string, value: number, unit: string }[]) => (
        <div className="space-y-2">
            <h4 className="font-semibold text-sm text-muted-foreground">{title}</h4>
            <Table>
                <TableBody>
                    {metrics.map((m) => (
                        <TableRow key={m.label} className="h-8">
                            <TableCell className="py-1 font-medium text-xs">{m.label}</TableCell>
                            <TableCell className="py-1 text-right font-mono text-xs">
                                {m.value?.toFixed(4)} {m.unit}
                            </TableCell>
                        </TableRow>
                    ))}
                </TableBody>
            </Table>
        </div>
    );

    return (
        <Card className="h-full overflow-hidden flex flex-col">
            <CardHeader className="pb-2 flex-shrink-0">
                <CardTitle className="text-lg">详细指标</CardTitle>
            </CardHeader>
            <CardContent className="flex-1 overflow-auto space-y-6">
                <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                    <div className="space-y-6">
                        {renderTable("一般性能", generalMetrics)}
                        {renderTable("成本分解", costMetrics)}
                    </div>
                    <div className="space-y-6">
                        <div className="w-full flex flex-col" style={{ minHeight: '200px', height: '200px' }}>
                            <h4 className="font-semibold text-sm text-muted-foreground mb-2 text-center">损耗分解</h4>
                            <div className="flex-1" style={{ minHeight: '150px', minWidth: '0' }}>
                                <ResponsiveContainer width="100%" height="100%" minHeight={150}>
                                    <PieChart>
                                        <Pie
                                            data={lossChartData}
                                            cx="50%"
                                            cy="50%"
                                            innerRadius={40}
                                            outerRadius={70}
                                            fill="#8884d8"
                                            paddingAngle={5}
                                            dataKey="value"
                                        >
                                            {lossChartData.map((entry, index) => (
                                                <Cell key={`cell-${index}`} fill={COLORS[index % COLORS.length]} />
                                            ))}
                                        </Pie>
                                        <Tooltip formatter={(value: number) => `${value.toFixed(2)} W`} />
                                        <Legend verticalAlign="bottom" height={36} iconSize={8} wrapperStyle={{ fontSize: '10px' }} />
                                    </PieChart>
                                </ResponsiveContainer>
                            </div>
                        </div>
                        {renderTable("损耗值", lossMetrics)}
                    </div>
                </div>
            </CardContent>
        </Card>
    );
}
