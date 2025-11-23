'use client';

import React from 'react';
import { DesignData } from '@/lib/DesignData';
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from '@/components/ui/table';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';

interface ParameterTableProps {
    data: DesignData;
}

export default function ParameterTable({ data }: ParameterTableProps) {
    const formatValue = (val: any) => {
        if (val === null || val === undefined) return '-';
        if (typeof val === 'number') return val.toFixed(4);
        if (typeof val === 'object' && val.value !== undefined) return formatValue(val.value);
        if (Array.isArray(val)) return `[${val.map(v => typeof v === 'number' ? v.toFixed(2) : v).join(', ')}]`;
        return String(val);
    };

    const renderSection = (title: string, obj: any, filterComplex = true) => {
        const entries = Object.entries(obj).filter(([key, val]) => {
            if (filterComplex) {
                if (key === 'wily' || key === 'GeometricComponentsObjects') return false;
                if (typeof val === 'object' && val !== null && !('value' in val) && !Array.isArray(val)) return false;
            }
            return true;
        });

        if (entries.length === 0) return null;

        return (
            <Card className="h-fit">
                <CardHeader className="pb-2">
                    <CardTitle className="text-lg">{title}</CardTitle>
                </CardHeader>
                <CardContent>
                    <Table>
                        <TableHeader>
                            <TableRow>
                                <TableHead>Parameter</TableHead>
                                <TableHead className="text-right">Value</TableHead>
                            </TableRow>
                        </TableHeader>
                        <TableBody>
                            {entries.map(([key, val]) => (
                                <TableRow key={key}>
                                    <TableCell className="font-medium max-w-[200px] truncate" title={key}>{key}</TableCell>
                                    <TableCell className="text-right font-mono text-xs">{formatValue(val)}</TableCell>
                                </TableRow>
                            ))}
                        </TableBody>
                    </Table>
                </CardContent>
            </Card>
        );
    };

    // Extract geometric details
    const renderGeometricDetails = () => {
        return Object.entries(data.GeometricComponentsObjects).map(([compName, compData]) => {
            // Filter out points (P1, P2...) and objects to show only scalar parameters
            const scalarData: Record<string, any> = {};
            Object.entries(compData).forEach(([k, v]) => {
                if (k.startsWith('P') && /\d/.test(k)) return; // Skip points like P1, P2
                if (k === 'location' || k === 'py/object') return;
                if (typeof v !== 'object' || v === null) {
                    scalarData[k] = v;
                }
            });

            return renderSection(`Geometry: ${compName}`, scalarData, false);
        });
    };

    return (
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
            {renderSection("General Parameters", {
                machine_type: data.machine_type,
                m: data.m,
                Qs: data.Qs,
                p: data.p,
                ps: data.ps,
                mec_power: data.mec_power,
                efficiency: data.guess_efficiency,
                power_factor: data.guess_power_factor
            })}

            {renderSection("Excitation (User)", data["EX-user"])}

            {renderSection("Geometric Parameters (User)", data["GP-user"])}

            {/* Render detailed geometric parameters for each component */}
            {renderGeometricDetails()}
        </div>
    );
}
