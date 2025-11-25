'use client';

import React from 'react';
import { GeometricComponentsObjects, GeometricComponent } from '@/lib/DesignData';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from '@/components/ui/table';
import { ScrollArea } from '@/components/ui/scroll-area';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { Checkbox } from '@/components/ui/checkbox';

interface GeometryDetailsProps {
    data: GeometricComponentsObjects;
    onComponentSelect?: (componentKey: string | null) => void;
    visibility?: Record<string, boolean>;
    onVisibilityChange?: (visibility: Record<string, boolean>) => void;
}

export default function GeometryDetails({ data, onComponentSelect, visibility = {}, onVisibilityChange }: GeometryDetailsProps) {
    const renderComponentDetails = (name: string, component: GeometricComponent | null) => {
        if (!component) {
            return (
                <div className="text-center text-muted-foreground py-8">
                    <p>No geometry data available for {name}</p>
                </div>
            );
        }

        // Separate points and other parameters
        const parameters: Record<string, any> = {};
        const points: Record<string, [number, number]> = {};

        Object.entries(component).forEach(([key, val]) => {
            if (key === 'py/object' || key === 'location' || key === 'color' || key === 'name') return;

            // Check if it's a point (starts with P and is array of 2 numbers)
            // Also handle P1p5 etc.
            if (key.startsWith('P') && Array.isArray(val) && val.length === 2 && typeof val[0] === 'number') {
                points[key] = val as [number, number];
            } else if (typeof val !== 'object' || val === null) {
                parameters[key] = val;
            }
        });

        return (
            <div className="space-y-4">
                <div>
                    <h4 className="text-sm font-semibold mb-2">Parameters</h4>
                    <Table>
                        <TableBody>
                            {Object.entries(parameters).map(([k, v]) => (
                                <TableRow key={k} className="h-8">
                                    <TableCell className="py-1 font-medium text-xs">{k}</TableCell>
                                    <TableCell className="py-1 text-right text-xs font-mono">
                                        {typeof v === 'number' ? v.toFixed(4) : String(v)}
                                    </TableCell>
                                </TableRow>
                            ))}
                        </TableBody>
                    </Table>
                </div>

                <div>
                    <h4 className="text-sm font-semibold mb-2">Points</h4>
                    <Table>
                        <TableHeader>
                            <TableRow className="h-8">
                                <TableHead className="h-8 py-1 text-xs">Point</TableHead>
                                <TableHead className="h-8 py-1 text-right text-xs">X</TableHead>
                                <TableHead className="h-8 py-1 text-right text-xs">Y</TableHead>
                            </TableRow>
                        </TableHeader>
                        <TableBody>
                            {Object.entries(points).sort((a, b) => a[0].localeCompare(b[0], undefined, { numeric: true })).map(([k, [x, y]]) => (
                                <TableRow key={k} className="h-8">
                                    <TableCell className="py-1 font-medium text-xs">{k}</TableCell>
                                    <TableCell className="py-1 text-right text-xs font-mono">{x.toFixed(3)}</TableCell>
                                    <TableCell className="py-1 text-right text-xs font-mono">{y.toFixed(3)}</TableCell>
                                </TableRow>
                            ))}
                        </TableBody>
                    </Table>
                </div>
            </div>
        );
    };

    return (
        <Card className="h-full flex flex-col">
            <CardHeader className="pb-2">
                <CardTitle className="text-lg">Geometry Details</CardTitle>
            </CardHeader>
            <CardContent className="flex-1 overflow-hidden p-0">
                <Tabs
                    defaultValue={Object.keys(data)[0]}
                    className="h-full flex flex-col"
                    onValueChange={(value) => onComponentSelect?.(value)}
                >
                    <div className="px-4 pt-2">
                        <TabsList className="w-full justify-start overflow-x-auto h-auto flex-wrap gap-1 bg-transparent p-0">
                            {Object.keys(data).map(key => {
                                const component = data[key as keyof GeometricComponentsObjects];
                                return (
                                    <div key={key} className="flex items-center gap-1">
                                        {onVisibilityChange && (
                                            <Checkbox
                                                checked={visibility[key] !== false}
                                                onCheckedChange={(checked: boolean) => {
                                                    onVisibilityChange({
                                                        ...visibility,
                                                        [key]: checked === true
                                                    });
                                                }}
                                                className="h-3 w-3"
                                            />
                                        )}
                                        <TabsTrigger
                                            value={key}
                                            className="data-[state=active]:bg-primary data-[state=active]:text-primary-foreground border text-xs px-2 py-1 h-auto"
                                        >
                                            {component?.name || key}
                                        </TabsTrigger>
                                    </div>
                                );
                            })}
                        </TabsList>
                    </div>

                    <ScrollArea className="flex-1 p-4">
                        {Object.entries(data).map(([key, component]) => {
                            if (!component) {
                                return (
                                    <TabsContent key={key} value={key} className="mt-0">
                                        <div className="text-center text-muted-foreground py-8">
                                            <p>No geometry data available for {key}</p>
                                        </div>
                                    </TabsContent>
                                );
                            }
                            return (
                                <TabsContent key={key} value={key} className="mt-0">
                                    {renderComponentDetails(key, component)}
                                </TabsContent>
                            );
                        })}
                    </ScrollArea>
                </Tabs>
            </CardContent>
        </Card>
    );
}
