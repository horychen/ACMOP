'use client';

import React from 'react';
import { GP, GPParameter } from '@/lib/DesignData';
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from '@/components/ui/card';
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from '@/components/ui/table';
import { Badge } from '@/components/ui/badge';
import { Input } from '@/components/ui/input';
import { Ruler, Cog, Zap } from 'lucide-react';

interface GPViewerProps {
    data?: GP;
    onParameterChange?: (key: string, value: number) => void;
}

const GPViewer: React.FC<GPViewerProps> = ({ data, onParameterChange }) => {
    if (!data) {
        return (
            <div className="flex items-center justify-center h-64 text-muted-foreground">
                No GP data available
            </div>
        );
    }

    // Categorize parameters
    const rotorParams = [
        'mm_r_ro', 'mm_r_ri', 'mm_d_pm', 'mm_d_ri', 'mm_d_rp', 'mm_d_rs',
        'deg_alpha_rm', 'deg_alpha_rs'
    ];
    const statorParams = [
        'mm_r_si', 'mm_r_so', 'mm_w_st', 'mm_d_st', 'mm_d_sy', 'mm_d_sts',
        'mm_d_sto', 'deg_alpha_st', 'deg_alpha_sto'
    ];
    const mechanicalParams = [
        'mm_d_mech_air_gap', 'mm_d_sleeve', 'split_ratio'
    ];

    const getTypeBadge = (type: string) => {
        switch (type) {
            case 'fixed':
                return <Badge variant="default" className="bg-blue-500">Fixed</Badge>;
            case 'free':
                return <Badge variant="default" className="bg-green-500">Free</Badge>;
            case 'derived':
                return <Badge variant="secondary">Derived</Badge>;
            default:
                return <Badge variant="outline">{type}</Badge>;
        }
    };

    const renderParameterTable = (paramKeys: string[], title: string, icon: React.ReactNode, colorClass: string) => {
        const params = paramKeys.map(key => ({ key, ...data[key] })).filter(p => p.value !== undefined);

        if (params.length === 0) return null;

        return (
            <Card>
                <CardHeader>
                    <div className="flex items-center gap-2">
                        <div className={`p-2 rounded-md ${colorClass}`}>
                            {icon}
                        </div>
                        <div>
                            <CardTitle className="text-lg">{title}</CardTitle>
                            <CardDescription>{params.length} parameters</CardDescription>
                        </div>
                    </div>
                </CardHeader>
                <CardContent>
                    <Table>
                        <TableHeader>
                            <TableRow>
                                <TableHead>Parameter</TableHead>
                                <TableHead>Type</TableHead>
                                <TableHead className="text-right">Value</TableHead>
                                <TableHead>Bounds</TableHead>
                            </TableRow>
                        </TableHeader>
                        <TableBody>
                            {params.map(param => (
                                <TableRow key={param.key}>
                                    <TableCell className="font-medium">
                                        <div className="flex flex-col">
                                            <span className="text-sm font-mono">{param.key}</span>
                                            <span className="text-xs text-muted-foreground">{param.description}</span>
                                        </div>
                                    </TableCell>
                                    <TableCell>
                                        {getTypeBadge(param.type)}
                                    </TableCell>
                                    <TableCell className="text-right font-mono">
                                        {onParameterChange && param.type === 'free' ? (
                                            <Input
                                                type="number"
                                                value={param.value}
                                                onChange={(e) => onParameterChange(param.key, parseFloat(e.target.value))}
                                                className="h-8 w-24 text-right ml-auto"
                                                step={0.1}
                                            />
                                        ) : (
                                            param.value?.toFixed(4)
                                        )}
                                    </TableCell>
                                    <TableCell className="text-xs text-muted-foreground">
                                        {param.bounds ? (
                                            <span className="font-mono">
                                                [{param.bounds[0].toFixed(2)}, {param.bounds[1].toFixed(2)}]
                                            </span>
                                        ) : (
                                            <span className="text-muted-foreground/50">—</span>
                                        )}
                                    </TableCell>
                                </TableRow>
                            ))}
                        </TableBody>
                    </Table>
                </CardContent>
            </Card>
        );
    };

    return (
        <div className="space-y-6">
            <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
                {renderParameterTable(
                    rotorParams,
                    'Rotor Parameters',
                    <Cog className="h-5 w-5 text-orange-600" />,
                    'bg-orange-100 dark:bg-orange-900/20'
                )}
                {renderParameterTable(
                    statorParams,
                    'Stator Parameters',
                    <Zap className="h-5 w-5 text-blue-600" />,
                    'bg-blue-100 dark:bg-blue-900/20'
                )}
            </div>

            <div className="grid grid-cols-1">
                {renderParameterTable(
                    mechanicalParams,
                    'Mechanical & Air Gap',
                    <Ruler className="h-5 w-5 text-green-600" />,
                    'bg-green-100 dark:bg-green-900/20'
                )}
            </div>
        </div>
    );
};

export default GPViewer;
