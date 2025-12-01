'use client';

import React, { useState } from 'react';
import { DesignData } from '@/lib/DesignData';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { ScrollArea } from "@/components/ui/scroll-area";
import CrossSectionViewer from './CrossSectionViewer';
import WindingLayoutViewer from './WindingLayoutViewer';
import ParameterTable from './ParameterTable';
import ExcitationViewer from './ExcitationViewer';
import GeometryDetails from './GeometryDetails';
import GeometryPlayground from './GeometryPlayground';
import GPViewer from './GPViewer';
import WindingDiagrams from './WindingDiagrams';
import { Maximize2, Minimize2, Layout, Activity, Settings, Zap, Database, CircleDashed } from 'lucide-react';

interface DesignVisualizerClientProps {
    data: DesignData;
}

type MainViewMode = 'cross-section' | 'winding';
type SidebarTab = 'details' | 'performance' | 'excitation' | 'gp' | 'diagrams';

export default function DesignVisualizerClient({ data }: DesignVisualizerClientProps) {
    const [mainView, setMainView] = useState<MainViewMode>('cross-section');
    const [sidebarTab, setSidebarTab] = useState<SidebarTab>('details');
    const [selectedComponent, setSelectedComponent] = useState<string | null>(null);
    const [visibility, setVisibility] = useState<Record<string, boolean>>({});
    const [isPlaygroundExpanded, setIsPlaygroundExpanded] = useState(false);

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
        <div className="h-[calc(100vh-4rem)] flex flex-col overflow-hidden bg-background">
            {/* Header */}
            <div className="flex items-center justify-between px-4 py-2 border-b shrink-0">
                <div className="flex items-center space-x-4">
                    <h1 className="text-xl font-bold tracking-tight">Design Visualizer</h1>
                    <div className="text-sm text-muted-foreground border-l pl-4">
                        {data.machine_type} (Qs={data.Qs}, p={data.p})
                    </div>
                </div>
                <div className="flex items-center space-x-2">
                    <div className="flex bg-muted rounded-lg p-1">
                        <Button
                            variant={mainView === 'cross-section' ? 'secondary' : 'ghost'}
                            size="sm"
                            onClick={() => setMainView('cross-section')}
                            className="h-7 text-xs"
                        >
                            <Layout className="w-3 h-3 mr-1" />
                            Cross Section
                        </Button>
                        <Button
                            variant={mainView === 'winding' ? 'secondary' : 'ghost'}
                            size="sm"
                            onClick={() => setMainView('winding')}
                            className="h-7 text-xs"
                        >
                            <Activity className="w-3 h-3 mr-1" />
                            Winding
                        </Button>
                    </div>
                </div>
            </div>

            {/* Main Content Area */}
            <div className="flex-1 flex overflow-hidden">
                {/* Left/Center Panel: Visualizer & Playground */}
                <div className="flex-1 flex flex-col min-w-0">
                    {/* Visualizer Area */}
                    <div className={`flex-1 relative bg-slate-50 dark:bg-slate-950 overflow-hidden transition-all duration-300 ${isPlaygroundExpanded ? 'h-1/3' : 'h-2/3'}`}>
                        <div className="absolute inset-0 p-4">
                            {mainView === 'cross-section' ? (
                                <CrossSectionViewer
                                    geometry={data.GeometricComponentsObjects}
                                    selectedComponent={selectedComponent}
                                    visibility={visibility}
                                    onVisibilityChange={setVisibility}
                                />
                            ) : (
                                <WindingLayoutViewer
                                    exUserData={data["EX-user"]}
                                    Qs={data.Qs}
                                    p={data.p}
                                    m={data.m}
                                />
                            )}
                        </div>
                    </div>

                    {/* Playground Area (Bottom) */}
                    <div className={`border-t bg-background flex flex-col transition-all duration-300 ${isPlaygroundExpanded ? 'flex-1' : 'h-[300px]'}`}>
                        <div className="flex items-center justify-between px-4 py-2 border-b bg-muted/30 shrink-0">
                            <span className="text-sm font-medium flex items-center">
                                <Database className="w-4 h-4 mr-2" />
                                Geometry Playground
                            </span>
                            <Button
                                variant="ghost"
                                size="icon"
                                className="h-6 w-6"
                                onClick={() => setIsPlaygroundExpanded(!isPlaygroundExpanded)}
                            >
                                {isPlaygroundExpanded ? <Minimize2 className="w-3 h-3" /> : <Maximize2 className="w-3 h-3" />}
                            </Button>
                        </div>
                        <div className="flex-1 overflow-hidden">
                            <GeometryPlayground data={data.GeometricComponentsObjects} />
                        </div>
                    </div>
                </div>

                {/* Right Sidebar: Details & Parameters */}
                <div className="w-[400px] border-l bg-background flex flex-col shrink-0">
                    <Tabs value={sidebarTab} onValueChange={(v) => setSidebarTab(v as SidebarTab)} className="flex-1 flex flex-col">
                        <div className="px-2 pt-2 border-b">
                            <TabsList className="w-full grid grid-cols-5">
                                <TabsTrigger value="details" title="Geometry Details"><Layout className="w-4 h-4" /></TabsTrigger>
                                <TabsTrigger value="performance" title="Parameters"><Settings className="w-4 h-4" /></TabsTrigger>
                                <TabsTrigger value="excitation" title="Excitation"><Zap className="w-4 h-4" /></TabsTrigger>
                                <TabsTrigger value="gp" title="GP"><Database className="w-4 h-4" /></TabsTrigger>
                                <TabsTrigger value="diagrams" title="Diagrams"><CircleDashed className="w-4 h-4" /></TabsTrigger>
                            </TabsList>
                        </div>

                        <div className="flex-1 overflow-hidden">
                            <ScrollArea className="h-full">
                                <div className="p-4">
                                    <TabsContent value="details" className="mt-0 space-y-4">
                                        <div className="space-y-1">
                                            <h3 className="font-semibold">Geometry Details</h3>
                                            <p className="text-xs text-muted-foreground">Component dimensions and properties</p>
                                        </div>
                                        <GeometryDetails
                                            data={data.GeometricComponentsObjects}
                                            onComponentSelect={setSelectedComponent}
                                            visibility={visibility}
                                            onVisibilityChange={setVisibility}
                                        />
                                    </TabsContent>

                                    <TabsContent value="performance" className="mt-0 space-y-4">
                                        <div className="space-y-1">
                                            <h3 className="font-semibold">Design Parameters</h3>
                                            <p className="text-xs text-muted-foreground">Key specifications and metrics</p>
                                        </div>
                                        <ParameterTable data={data} />
                                    </TabsContent>

                                    <TabsContent value="excitation" className="mt-0 space-y-4">
                                        <div className="space-y-1">
                                            <h3 className="font-semibold">Circuit Excitation</h3>
                                            <p className="text-xs text-muted-foreground">Voltage, current, and frequency</p>
                                        </div>
                                        <ExcitationViewer data={data["EX-user"]} />
                                    </TabsContent>

                                    <TabsContent value="gp" className="mt-0 space-y-4">
                                        <div className="space-y-1">
                                            <h3 className="font-semibold">Geometric Parameters</h3>
                                            <p className="text-xs text-muted-foreground">Optimization bounds and values</p>
                                        </div>
                                        <GPViewer data={data.GP} />
                                    </TabsContent>

                                    <TabsContent value="diagrams" className="mt-0 space-y-4">
                                        <div className="space-y-1">
                                            <h3 className="font-semibold">Winding Diagrams</h3>
                                            <p className="text-xs text-muted-foreground">Star of Slots and Connection Star</p>
                                        </div>
                                        <WindingDiagrams
                                            Qs={data.Qs}
                                            p={data.p}
                                            m={data.m || 3}
                                            layer_X_phases={data["EX-user"]?.wily?.layer_X_phases || []}
                                            layer_X_signs={data["EX-user"]?.wily?.layer_X_signs || []}
                                        />
                                    </TabsContent>
                                </div>
                            </ScrollArea>
                        </div>
                    </Tabs>
                </div>
            </div>
        </div>
    );
}
