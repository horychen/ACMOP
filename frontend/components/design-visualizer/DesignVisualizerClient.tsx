'use client';

import React, { useState } from 'react';
import { DesignData } from '@/lib/DesignData';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import CrossSectionViewer from './CrossSectionViewer';
import WindingLayoutViewer from './WindingLayoutViewer';
import ParameterTable from './ParameterTable';
import ExcitationViewer from './ExcitationViewer';
import GeometryDetails from './GeometryDetails';
import GeometryPlayground from './GeometryPlayground';
import GPViewer from './GPViewer';

interface DesignVisualizerClientProps {
    data: DesignData;
}

type TabValue = 'geometry' | 'winding' | 'performance' | 'excitation' | 'gp';

export default function DesignVisualizerClient({ data }: DesignVisualizerClientProps) {
    const [activeTab, setActiveTab] = useState<TabValue>('geometry');
    const [selectedComponent, setSelectedComponent] = useState<string | null>(null);
    const [visibility, setVisibility] = useState<Record<string, boolean>>({});

    // Initialize visibility state
    React.useEffect(() => {
        if (!data.GeometricComponentsObjects) return;
        const initialVisibility: Record<string, boolean> = {};
        Object.keys(data.GeometricComponentsObjects).forEach(key => {
            // Default all to visible except sleeve
            initialVisibility[key] = key.toLowerCase().includes('sleeve') ? false : true;
        });
        setVisibility(initialVisibility);
    }, [data.GeometricComponentsObjects]);

    return (
        <div className="container mx-auto p-4 space-y-6">
            <div className="flex flex-col space-y-2">
                <h1 className="text-3xl font-bold tracking-tight">Design Visualizer</h1>
                <p className="text-muted-foreground">
                    Visualizing {data.machine_type} Design (Qs={data.Qs}, p={data.p})
                </p>
            </div>

            <div className="space-y-4">
                <div className="flex space-x-2 border-b pb-2 overflow-x-auto">
                    <Button
                        variant={activeTab === 'geometry' ? 'default' : 'ghost'}
                        onClick={() => setActiveTab('geometry')}
                    >
                        Geometry & Playground
                    </Button>
                    <Button
                        variant={activeTab === 'winding' ? 'default' : 'ghost'}
                        onClick={() => setActiveTab('winding')}
                    >
                        Winding Layout
                    </Button>
                    <Button
                        variant={activeTab === 'performance' ? 'default' : 'ghost'}
                        onClick={() => setActiveTab('performance')}
                    >
                        Performance & Parameters
                    </Button>
                    <Button
                        variant={activeTab === 'excitation' ? 'default' : 'ghost'}
                        onClick={() => setActiveTab('excitation')}
                    >
                        Circuit Excitation
                    </Button>
                    <Button
                        variant={activeTab === 'gp' ? 'default' : 'ghost'}
                        onClick={() => setActiveTab('gp')}
                    >
                        GP Parameters
                    </Button>
                </div>

                <div className="mt-4">
                    {activeTab === 'geometry' && (
                        <div className="grid grid-cols-1 lg:grid-cols-3 gap-4">
                            {/* Top row: Cross Section spanning 2 columns */}
                            <div className="lg:col-span-2 h-[600px]">
                                <Card className="h-full flex flex-col">
                                    <CardHeader>
                                        <CardTitle>Cross Section</CardTitle>
                                        <CardDescription>
                                            Geometric representation of the rotor and stator.
                                        </CardDescription>
                                    </CardHeader>
                                    <CardContent className="flex-1 bg-slate-50 dark:bg-slate-900 rounded-md overflow-hidden relative m-4 mt-0">
                                        <CrossSectionViewer
                                            geometry={data.GeometricComponentsObjects}
                                            selectedComponent={selectedComponent}
                                            visibility={visibility}
                                            onVisibilityChange={setVisibility}
                                        />
                                    </CardContent>
                                </Card>
                            </div>

                            {/* Top right: Geometry Details */}
                            <div className="lg:col-span-1 h-[600px]">
                                <GeometryDetails
                                    data={data.GeometricComponentsObjects}
                                    onComponentSelect={setSelectedComponent}
                                    visibility={visibility}
                                    onVisibilityChange={setVisibility}
                                />
                            </div>

                            {/* Bottom: Playground Editor spanning all 3 columns */}
                            <div className="lg:col-span-3 h-[600px]">
                                <GeometryPlayground data={data.GeometricComponentsObjects} />
                            </div>
                        </div>
                    )}

                    {activeTab === 'winding' && (
                        <Card>
                            <CardHeader>
                                <CardTitle>Winding Layout</CardTitle>
                                <CardDescription>
                                    Stator winding distribution and phases.
                                </CardDescription>
                            </CardHeader>
                            <CardContent>
                                <WindingLayoutViewer exUserData={data["EX-user"]} Qs={data.Qs} p={data.p} />
                            </CardContent>
                        </Card>
                    )}

                    {activeTab === 'performance' && (
                        <Card>
                            <CardHeader>
                                <CardTitle>Design Parameters</CardTitle>
                                <CardDescription>
                                    Key specifications and performance metrics.
                                </CardDescription>
                            </CardHeader>
                            <CardContent>
                                <ParameterTable data={data} />
                            </CardContent>
                        </Card>
                    )}

                    {activeTab === 'excitation' && (
                        <Card>
                            <CardHeader>
                                <CardTitle>Circuit Excitation</CardTitle>
                                <CardDescription>
                                    Voltage, current, and frequency settings.
                                </CardDescription>
                            </CardHeader>
                            <CardContent>
                                <ExcitationViewer data={data["EX-user"]} />
                            </CardContent>
                        </Card>
                    )}

                    {activeTab === 'gp' && (
                        <Card>
                            <CardHeader>
                                <CardTitle>Geometric Parameters (GP)</CardTitle>
                                <CardDescription>
                                    All machine design parameters with optimization bounds.
                                </CardDescription>
                            </CardHeader>
                            <CardContent>
                                <GPViewer data={data.GP} />
                            </CardContent>
                        </Card>
                    )}
                </div>
            </div>
        </div>
    );
}
