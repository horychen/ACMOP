"use client";

import React, { useState, useEffect } from 'react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { Loader2, Code2, Eye, AlertCircle, Download, RefreshCw } from 'lucide-react';
import axios from 'axios';
import DesignVisualizerClient from '@/components/design-visualizer/DesignVisualizerClient';
import { DesignData, parseDesignData, parseGPData } from '@/lib/DesignData';
import { Badge } from '@/components/ui/badge';

interface MachineSpec {
    machine_type: string;
    DPNV_or_SEPA?: boolean;
    m: number;
    Qs: number;
    p: number;
    ps: number;
    coil_pitch_y?: number;
    mec_power: number;
    guess_air_gap_flux_density_Bg?: number;
    guess_stator_tooth_flux_density_Bst?: number;
    guess_stator_yoke_flux_density_Bsy?: number;
    guess_efficiency?: number;
    guess_power_factor?: number;
    bool_skew_stator?: boolean;
    bool_skew_rotor?: boolean;
    no_segmented_magnets?: number;
    GeometricComponentsObjects?: any;
    "EX-user": any;
    "GP-user": any;
    [key: string]: any;
}

export default function DeveloperPage() {
    const [availableSpecs, setAvailableSpecs] = useState<Record<string, MachineSpec>>({});
    const [selectedSpecKey, setSelectedSpecKey] = useState<string>('');
    const [selectedSpec, setSelectedSpec] = useState<MachineSpec | null>(null);
    const [designData, setDesignData] = useState<DesignData | null>(null);
    const [isLoading, setIsLoading] = useState(false);
    const [error, setError] = useState<string | null>(null);
    const [activeView, setActiveView] = useState<'visualizer' | 'json'>('visualizer');
    const [jsonView, setJsonView] = useState<string>('');

    // Load available specifications
    useEffect(() => {
        const loadSpecs = async () => {
            setIsLoading(true);
            setError(null);
            try {
                const response = await axios.get('/api/machine-specs');
                setAvailableSpecs(response.data);
                
                // Auto-select first spec if available
                const keys = Object.keys(response.data);
                if (keys.length > 0 && !selectedSpecKey) {
                    setSelectedSpecKey(keys[0]);
                }
            } catch (err: any) {
                console.error('Failed to load machine specs:', err);
                setError(err.response?.data?.error || 'Failed to load machine specifications');
            } finally {
                setIsLoading(false);
            }
        };

        loadSpecs();
    }, []);

    // Process selected specification
    useEffect(() => {
        if (!selectedSpecKey || !availableSpecs[selectedSpecKey]) {
            setSelectedSpec(null);
            setDesignData(null);
            setJsonView('');
            return;
        }

        const spec = availableSpecs[selectedSpecKey];
        setSelectedSpec(spec);

            // Convert to DesignData format
        try {
            // Parse GP-user data from machine_specifications.json format
            const parseGPUserData = (gpUser: any): any => {
                if (!gpUser) return undefined;
                
                const gp: any = {};
                Object.entries(gpUser).forEach(([key, param]: [string, any]) => {
                    if (param && typeof param === 'object' && param.type) {
                        // Handle null values - only include if value is not null or if it's a derived parameter
                        if (param.value !== null || param.type === 'derived') {
                            gp[key] = {
                                type: param.type,
                                description: key, // Use key as description if not provided
                                value: param.value,
                                bounds: param.bounds || null
                            };
                        }
                    }
                });
                return Object.keys(gp).length > 0 ? gp : undefined;
            };

            const processed: DesignData = {
                machine_type: spec.machine_type,
                m: spec.m || 3,
                Qs: spec.Qs,
                p: spec.p,
                ps: spec.ps,
                coil_pitch_y: spec.coil_pitch_y,
                mec_power: spec.mec_power,
                guess_efficiency: spec.guess_efficiency,
                guess_power_factor: spec.guess_power_factor,
                guess_air_gap_flux_density_Bg: spec.guess_air_gap_flux_density_Bg,
                guess_stator_tooth_flux_density_Bst: spec.guess_stator_tooth_flux_density_Bst,
                guess_stator_yoke_flux_density_Bsy: spec.guess_stator_yoke_flux_density_Bsy,
                bool_skew_stator: spec.bool_skew_stator,
                bool_skew_rotor: spec.bool_skew_rotor,
                no_segmented_magnets: spec.no_segmented_magnets,
                GeometricComponentsObjects: spec.GeometricComponentsObjects || {
                    rotorCore: null,
                    shaft: null,
                    rotorMagnet: null,
                    sleeve: null,
                    statorCore: null,
                    coils: null
                },
                "EX-user": spec["EX-user"] || {},
                "GP-user": spec["GP-user"] || {},
                GP: parseGPUserData(spec["GP-user"])
            };

            setDesignData(processed);
            
            // Format JSON for display
            setJsonView(JSON.stringify(spec, null, 2));
        } catch (err) {
            console.error('Error processing specification:', err);
            setError('Failed to process specification data');
        }
    }, [selectedSpecKey, availableSpecs]);

    const handleDownloadJson = () => {
        if (!jsonView) return;
        const blob = new Blob([jsonView], { type: 'application/json' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = `${selectedSpecKey || 'spec'}.json`;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
    };

    const handleRefresh = () => {
        window.location.reload();
    };

    return (
        <div className="space-y-6">
            {/* Header */}
            <div className="flex items-center justify-between">
                <div>
                    <h1 className="text-3xl font-bold tracking-tight">开发者页面</h1>
                    <p className="text-muted-foreground mt-1">
                        可视化 machine_specifications.json 中的设计规格
                    </p>
                </div>
                <div className="flex items-center space-x-2">
                    <Button
                        variant="outline"
                        size="sm"
                        onClick={handleRefresh}
                        disabled={isLoading}
                    >
                        <RefreshCw className={`h-4 w-4 mr-2 ${isLoading ? 'animate-spin' : ''}`} />
                        刷新
                    </Button>
                    {jsonView && (
                        <Button
                            variant="outline"
                            size="sm"
                            onClick={handleDownloadJson}
                        >
                            <Download className="h-4 w-4 mr-2" />
                            下载 JSON
                        </Button>
                    )}
                </div>
            </div>

            {/* Error Display */}
            {error && (
                <Card className="border-destructive">
                    <CardContent className="pt-6">
                        <div className="flex items-center space-x-2 text-destructive">
                            <AlertCircle className="h-5 w-5" />
                            <span>{error}</span>
                        </div>
                    </CardContent>
                </Card>
            )}

            {/* Specification Selector */}
            <Card>
                <CardHeader>
                    <CardTitle>选择设计规格</CardTitle>
                    <CardDescription>
                        从 machine_specifications.json 中选择要可视化的设计
                    </CardDescription>
                </CardHeader>
                <CardContent>
                    <div className="flex items-center space-x-4">
                        <div className="flex-1">
                            <Select
                                value={selectedSpecKey}
                                onValueChange={setSelectedSpecKey}
                                disabled={isLoading || Object.keys(availableSpecs).length === 0}
                            >
                                <SelectTrigger>
                                    <SelectValue placeholder="选择机器规格..." />
                                </SelectTrigger>
                                <SelectContent>
                                    {Object.keys(availableSpecs).map((key) => (
                                        <SelectItem key={key} value={key}>
                                            {key}
                                        </SelectItem>
                                    ))}
                                </SelectContent>
                            </Select>
                        </div>
                        {isLoading && (
                            <Loader2 className="h-5 w-5 animate-spin text-muted-foreground" />
                        )}
                        {selectedSpec && (
                            <div className="flex items-center space-x-2">
                                <Badge variant="outline">{selectedSpec.machine_type}</Badge>
                                <Badge variant="secondary">Qs={selectedSpec.Qs}</Badge>
                                <Badge variant="secondary">p={selectedSpec.p}</Badge>
                                <Badge variant="secondary">ps={selectedSpec.ps}</Badge>
                            </div>
                        )}
                    </div>
                </CardContent>
            </Card>

            {/* Main Content */}
            {selectedSpec && designData && (
                <Tabs value={activeView} onValueChange={(v) => setActiveView(v as 'visualizer' | 'json')}>
                    <TabsList>
                        <TabsTrigger value="visualizer">
                            <Eye className="h-4 w-4 mr-2" />
                            可视化视图
                        </TabsTrigger>
                        <TabsTrigger value="json">
                            <Code2 className="h-4 w-4 mr-2" />
                            JSON 数据
                        </TabsTrigger>
                    </TabsList>

                    <TabsContent value="visualizer" className="mt-4">
                        <Card>
                            <CardHeader>
                                <CardTitle>设计可视化</CardTitle>
                                <CardDescription>
                                    完整的设计规格可视化：几何形状、绕组布局、参数等
                                </CardDescription>
                            </CardHeader>
                            <CardContent>
                                <DesignVisualizerClient data={designData} />
                            </CardContent>
                        </Card>
                    </TabsContent>

                    <TabsContent value="json" className="mt-4">
                        <Card>
                            <CardHeader>
                                <CardTitle>原始 JSON 数据</CardTitle>
                                <CardDescription>
                                    machine_specifications.json 中的完整数据结构
                                </CardDescription>
                            </CardHeader>
                            <CardContent>
                                <div className="relative">
                                    <pre className="bg-muted p-4 rounded-lg overflow-auto max-h-[800px] text-sm font-mono">
                                        <code>{jsonView}</code>
                                    </pre>
                                </div>
                            </CardContent>
                        </Card>
                    </TabsContent>
                </Tabs>
            )}

            {/* Empty State */}
            {!selectedSpec && !isLoading && !error && (
                <Card>
                    <CardContent className="pt-6">
                        <div className="text-center text-muted-foreground py-8">
                            <Code2 className="h-12 w-12 mx-auto mb-4 opacity-50" />
                            <p>请从上方选择器中选择一个设计规格</p>
                        </div>
                    </CardContent>
                </Card>
            )}
        </div>
    );
}

