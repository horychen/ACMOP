'use client';

import React, { useState, useMemo } from 'react';
import { DesignData } from '@/lib/DesignData';
import { generateGeometryFromGP } from '@/lib/geometryGenerator';
import GPViewer from '@/components/design-visualizer/GPViewer';
import CrossSectionViewer from '@/components/design-visualizer/CrossSectionViewer';
import PerformanceRadar from './PerformanceRadar';
import PerformanceMetrics from './PerformanceMetrics';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { Alert, AlertDescription, AlertTitle } from '@/components/ui/alert';
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from '@/components/ui/table';
import { Info, Settings, BarChart } from 'lucide-react';

interface AnalyzerDashboardProps {
    initialData: DesignData;
    fileName: string;
}

export default function AnalyzerDashboard({ initialData, fileName }: AnalyzerDashboardProps) {
    const [data, setData] = useState<DesignData>(initialData);

    const handleParameterChange = (key: string, value: number) => {
        setData(prev => {
            if (!prev.GP) return prev;
            return {
                ...prev,
                GP: {
                    ...prev.GP,
                    [key]: {
                        ...prev.GP[key],
                        value: value
                    }
                }
            };
        });
    };

    // Check if geometry data is available
    const hasOriginalGeometry = data.GeometricComponentsObjects && 
        Object.values(data.GeometricComponentsObjects).some(comp => comp !== null);

    // Generate geometry from GP if original is not available
    const generatedGeometry = useMemo(() => {
        if (hasOriginalGeometry) return null;
        return generateGeometryFromGP(data);
    }, [data, hasOriginalGeometry]);

    // Use original geometry if available, otherwise use generated
    const displayGeometry = hasOriginalGeometry ? data.GeometricComponentsObjects : generatedGeometry;
    const hasGeometry = hasOriginalGeometry || generatedGeometry !== null;
    const isGenerated = !hasOriginalGeometry && generatedGeometry !== null;

    // Format input parameters for display
    const formatInputValue = (value: any): string => {
        if (value === null || value === undefined) return 'N/A';
        if (typeof value === 'number') {
            if (Math.abs(value) < 0.001) return value.toExponential(2);
            if (Math.abs(value) > 1000) return value.toExponential(2);
            return value.toFixed(4);
        }
        return String(value);
    };

    return (
        <div className="space-y-4 h-[calc(100vh-8rem)] flex flex-col overflow-hidden">
            <div className="flex items-center justify-between flex-shrink-0">
                <div>
                    <h2 className="text-2xl font-bold tracking-tight">{fileName}</h2>
                    <p className="text-muted-foreground">
                        {data.machine_type} | {data.mec_power / 1000}kW | {data.p} Pole Pairs | Qs: {data.Qs}
                    </p>
                </div>
            </div>

            <div className="flex-1 overflow-hidden">
                <Tabs defaultValue="analysis" className="h-full flex flex-col">
                    <TabsList className="flex-shrink-0">
                        <TabsTrigger value="analysis">分析视图</TabsTrigger>
                        <TabsTrigger value="inputs">输入参数</TabsTrigger>
                        <TabsTrigger value="outputs">FEA输出</TabsTrigger>
                    </TabsList>

                    {/* Main Analysis View */}
                    <TabsContent value="analysis" className="flex-1 mt-4 overflow-hidden">
                        <div className="grid grid-cols-1 lg:grid-cols-12 gap-4 h-full">
                            {/* Left Column: Parameters (3 cols) */}
                            <div className="lg:col-span-3 flex flex-col gap-4 overflow-hidden">
                                <Card className="flex-1 overflow-hidden flex flex-col">
                                    <CardHeader className="py-3 px-4 bg-muted/50">
                                        <CardTitle className="text-sm font-medium">几何参数 (GP)</CardTitle>
                                    </CardHeader>
                                    <CardContent className="flex-1 overflow-y-auto p-0">
                                        <div className="p-4">
                                            {data.GP ? (
                                                <GPViewer data={data.GP} onParameterChange={handleParameterChange} />
                                            ) : (
                                                <p className="text-sm text-muted-foreground">无几何参数数据</p>
                                            )}
                                        </div>
                                    </CardContent>
                                </Card>
                            </div>

                            {/* Middle Column: Visualization (5 cols) */}
                            <div className="lg:col-span-5 flex flex-col gap-4 overflow-hidden">
                                <Card className="flex-1 overflow-hidden flex flex-col">
                                    <CardHeader className="py-3 px-4 bg-muted/50 flex flex-row items-center justify-between">
                                        <CardTitle className="text-sm font-medium">横截面视图</CardTitle>
                                        <span className="text-xs text-muted-foreground">
                                            {hasGeometry ? (isGenerated ? "✓ 从参数生成" : "✓ 几何数据已加载") : "✗ 无几何数据"}
                                        </span>
                                    </CardHeader>
                                    <CardContent className="flex-1 p-0 relative bg-slate-50 dark:bg-slate-950/50 min-h-0">
                                        {hasGeometry && displayGeometry ? (
                                            <>
                                                <CrossSectionViewer
                                                    geometry={displayGeometry}
                                                    showParameters={isGenerated}
                                                    gpData={isGenerated ? data.GP : undefined}
                                                />
                                                <div className="absolute bottom-2 left-2 right-2 z-10">
                                                    <Alert variant="default" className="bg-background/80 backdrop-blur-sm border-primary/20">
                                                        <Info className="h-4 w-4" />
                                                        <AlertTitle className="text-xs">提示</AlertTitle>
                                                        <AlertDescription className="text-xs">
                                                            {isGenerated 
                                                                ? "此几何图形是根据GP参数自动生成的近似视图。编辑参数会更新数值，但不会实时更新几何图形。"
                                                                : "编辑参数会更新数值，但不会在本地重新生成几何图形。"}
                                                        </AlertDescription>
                                                    </Alert>
                                                </div>
                                            </>
                                        ) : (
                                            <div className="h-full flex items-center justify-center">
                                                <div className="text-center space-y-2">
                                                    <p className="text-muted-foreground text-sm">无几何数据可用</p>
                                                    <p className="text-xs text-muted-foreground">
                                                        JSON文件中的GeometricComponentsObjects为null，且无法从GP参数生成几何图形
                                                    </p>
                                                </div>
                                            </div>
                                        )}
                                    </CardContent>
                                </Card>
                            </div>

                            {/* Right Column: Performance (4 cols) */}
                            <div className="lg:col-span-4 flex flex-col gap-4 overflow-hidden">
                                <Tabs defaultValue="radar" className="flex-1 flex flex-col overflow-hidden">
                                    <div className="flex items-center justify-between mb-2 flex-shrink-0">
                                        <h3 className="font-semibold text-sm">FEA性能</h3>
                                        <TabsList className="h-8">
                                            <TabsTrigger value="radar" className="text-xs">雷达图</TabsTrigger>
                                            <TabsTrigger value="metrics" className="text-xs">详细指标</TabsTrigger>
                                        </TabsList>
                                    </div>

                                    <TabsContent value="radar" className="flex-1 mt-0 overflow-hidden" style={{ minHeight: '400px' }}>
                                        {data["FEA_Evaluated_Performance--1-Initial"] ? (
                                            <div className="h-full w-full" style={{ minHeight: '400px' }}>
                                                <PerformanceRadar data={data["FEA_Evaluated_Performance--1-Initial"]} />
                                            </div>
                                        ) : (
                                            <Card className="h-full flex items-center justify-center">
                                                <p className="text-muted-foreground text-sm">无FEA数据可用</p>
                                            </Card>
                                        )}
                                    </TabsContent>

                                    <TabsContent value="metrics" className="flex-1 mt-0 overflow-hidden" style={{ minHeight: '400px' }}>
                                        {data["FEA_Evaluated_Performance--1-Initial"] ? (
                                            <div className="h-full w-full" style={{ minHeight: '400px' }}>
                                                <PerformanceMetrics data={data["FEA_Evaluated_Performance--1-Initial"]} />
                                            </div>
                                        ) : (
                                            <Card className="h-full flex items-center justify-center">
                                                <p className="text-muted-foreground text-sm">无FEA数据可用</p>
                                            </Card>
                                        )}
                                    </TabsContent>
                                </Tabs>
                            </div>
                        </div>
                    </TabsContent>

                    {/* Inputs Tab */}
                    <TabsContent value="inputs" className="flex-1 mt-4 overflow-auto">
                        <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
                            {/* EX-user Inputs */}
                            <Card>
                                <CardHeader>
                                    <CardTitle className="flex items-center gap-2">
                                        <Settings className="h-5 w-5" />
                                        激励参数 (EX-user)
                                    </CardTitle>
                                </CardHeader>
                                <CardContent>
                                    <Table>
                                        <TableHeader>
                                            <TableRow>
                                                <TableHead>参数</TableHead>
                                                <TableHead className="text-right">数值</TableHead>
                                            </TableRow>
                                        </TableHeader>
                                        <TableBody>
                                            {data["EX-user"] && Object.entries(data["EX-user"]).map(([key, value]) => {
                                                if (key === 'wily' || typeof value === 'object') return null;
                                                return (
                                                    <TableRow key={key}>
                                                        <TableCell className="font-medium">{key}</TableCell>
                                                        <TableCell className="text-right font-mono text-sm">
                                                            {formatInputValue(value)}
                                                        </TableCell>
                                                    </TableRow>
                                                );
                                            })}
                                        </TableBody>
                                    </Table>
                                </CardContent>
                            </Card>

                            {/* GP-user Inputs */}
                            <Card>
                                <CardHeader>
                                    <CardTitle className="flex items-center gap-2">
                                        <Settings className="h-5 w-5" />
                                        几何参数 (GP-user)
                                    </CardTitle>
                                </CardHeader>
                                <CardContent>
                                    <Table>
                                        <TableHeader>
                                            <TableRow>
                                                <TableHead>参数</TableHead>
                                                <TableHead className="text-right">数值</TableHead>
                                            </TableRow>
                                        </TableHeader>
                                        <TableBody>
                                            {data["GP-user"] && Object.entries(data["GP-user"]).map(([key, value]) => {
                                                const val = typeof value === 'object' && value !== null ? value.value : value;
                                                return (
                                                    <TableRow key={key}>
                                                        <TableCell className="font-medium">{key}</TableCell>
                                                        <TableCell className="text-right font-mono text-sm">
                                                            {formatInputValue(val)}
                                                        </TableCell>
                                                    </TableRow>
                                                );
                                            })}
                                        </TableBody>
                                    </Table>
                                </CardContent>
                            </Card>

                            {/* General Machine Parameters */}
                            <Card className="lg:col-span-2">
                                <CardHeader>
                                    <CardTitle className="flex items-center gap-2">
                                        <Settings className="h-5 w-5" />
                                        机器基本参数
                                    </CardTitle>
                                </CardHeader>
                                <CardContent>
                                    <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
                                        <div>
                                            <p className="text-sm text-muted-foreground">机器类型</p>
                                            <p className="font-semibold">{data.machine_type}</p>
                                        </div>
                                        <div>
                                            <p className="text-sm text-muted-foreground">相数 (m)</p>
                                            <p className="font-semibold">{data.m}</p>
                                        </div>
                                        <div>
                                            <p className="text-sm text-muted-foreground">定子槽数 (Qs)</p>
                                            <p className="font-semibold">{data.Qs}</p>
                                        </div>
                                        <div>
                                            <p className="text-sm text-muted-foreground">极对数 (p)</p>
                                            <p className="font-semibold">{data.p}</p>
                                        </div>
                                        <div>
                                            <p className="text-sm text-muted-foreground">极数 (ps)</p>
                                            <p className="font-semibold">{data.ps}</p>
                                        </div>
                                        <div>
                                            <p className="text-sm text-muted-foreground">机械功率 (kW)</p>
                                            <p className="font-semibold">{(data.mec_power / 1000).toFixed(2)}</p>
                                        </div>
                                        <div>
                                            <p className="text-sm text-muted-foreground">预估效率</p>
                                            <p className="font-semibold">{(data.guess_efficiency * 100).toFixed(2)}%</p>
                                        </div>
                                        <div>
                                            <p className="text-sm text-muted-foreground">功率因数</p>
                                            <p className="font-semibold">{data.guess_power_factor?.toFixed(3) || 'N/A'}</p>
                                        </div>
                                    </div>
                                </CardContent>
                            </Card>
                        </div>
                    </TabsContent>

                    {/* Outputs Tab */}
                    <TabsContent value="outputs" className="flex-1 mt-4 overflow-auto">
                        {data["FEA_Evaluated_Performance--1-Initial"] ? (
                            <div className="space-y-4">
                                <Card>
                                    <CardHeader>
                                        <CardTitle className="flex items-center gap-2">
                                            <BarChart className="h-5 w-5" />
                                            FEA评估性能 - 完整输出
                                        </CardTitle>
                                    </CardHeader>
                                    <CardContent>
                                        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
                                            {Object.entries(data["FEA_Evaluated_Performance--1-Initial"]).map(([key, value]) => {
                                                if (typeof value === 'object' && value !== null) return null;
                                                return (
                                                    <div key={key} className="p-3 border rounded-lg">
                                                        <p className="text-xs text-muted-foreground mb-1">{key}</p>
                                                        <p className="font-mono font-semibold text-sm">
                                                            {formatInputValue(value)}
                                                        </p>
                                                    </div>
                                                );
                                            })}
                                        </div>
                                    </CardContent>
                                </Card>
                                <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
                                    <div className="w-full" style={{ minHeight: '400px', height: '400px' }}>
                                        <PerformanceRadar data={data["FEA_Evaluated_Performance--1-Initial"]} />
                                    </div>
                                    <div className="w-full" style={{ minHeight: '400px', height: '400px' }}>
                                        <PerformanceMetrics data={data["FEA_Evaluated_Performance--1-Initial"]} />
                                    </div>
                                </div>
                            </div>
                        ) : (
                            <Card className="h-full flex items-center justify-center">
                                <div className="text-center space-y-2">
                                    <p className="text-muted-foreground">无FEA输出数据可用</p>
                                    <p className="text-xs text-muted-foreground">
                                        JSON文件中缺少"FEA_Evaluated_Performance--1-Initial"字段
                                    </p>
                                </div>
                            </Card>
                        )}
                    </TabsContent>
                </Tabs>
            </div>
        </div>
    );
}
