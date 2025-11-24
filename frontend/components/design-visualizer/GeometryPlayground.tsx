'use client';

import React, { useState, useEffect, useCallback, useRef } from 'react';
import Editor from '@monaco-editor/react';
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { Play, RotateCcw, ZoomIn, ZoomOut, Maximize } from "lucide-react";
import CrossSectionViewer from './CrossSectionViewer';
import { GeometricComponentsObjects, GeometricComponent, getPoints } from '@/lib/DesignData';

const DEFAULT_CODE = `// Define points
const P1 = [10, 0];
const P2 = [15, 5];
const P3 = [5, 10];

// Draw a shape with lines
MoveTo(P1);
LineTo(P2);
LineTo(P3);

// Draw an arc (independent of previous point)
// ArcTo(center, radius, startAngle, endAngle)
// Angles are in radians, automatically moves to arc start
ArcTo([0, 0], 10, Math.PI/2, Math.PI);

// Continue drawing
LineTo(P1);
`;

interface GeometryPlaygroundProps {
    data?: GeometricComponentsObjects;
}

export default function GeometryPlayground({ data }: GeometryPlaygroundProps) {
    const [code, setCode] = useState(DEFAULT_CODE);
    const [geometry, setGeometry] = useState<GeometricComponentsObjects | null>(null);
    const [error, setError] = useState<string | null>(null);
    const [selectedComponent, setSelectedComponent] = useState<string>("custom");
    const [variableNames, setVariableNames] = useState<Record<string, string>>({}); // Track point variable names

    // Reverse engineer list_region to DSL
    const generateCodeFromComponent = (component: GeometricComponent) => {
        let generatedCode = `// Generated code for ${component.name}\n\n`;

        // 1. Extract predefined Points (P1, P2, etc.)
        const pointMap = new Map<string, string>(); // "x,y" -> "P1"
        const points: { name: string, val: [number, number] }[] = [];

        Object.keys(component).forEach(key => {
            if ((/^P\d+/.test(key) || /^P[A-Za-z]+/.test(key)) && Array.isArray(component[key])) {
                const p = component[key] as [number, number];
                if (p.length === 2 && typeof p[0] === 'number') {
                    const name = key;
                    points.push({ name, val: p });
                    pointMap.set(`${p[0]},${p[1]}`, name);
                }
            }
        });

        // 2. Collect all unique coordinates from list_region
        const allCoords = new Set<string>();
        if (component.list_region) {
            component.list_region.forEach((region: any) => {
                if (!Array.isArray(region)) return;

                region.forEach((segment: any) => {
                    let moveTo = segment.move_to;
                    if (moveTo && moveTo['py/tuple']) moveTo = moveTo['py/tuple'];
                    else if (Array.isArray(moveTo)) moveTo = moveTo;

                    if (moveTo) {
                        allCoords.add(`${moveTo[0]},${moveTo[1]}`);
                    }

                    if (segment.line_to) {
                        let lineTo = segment.line_to;
                        if (lineTo && lineTo['py/tuple']) lineTo = lineTo['py/tuple'];
                        else if (Array.isArray(lineTo)) lineTo = lineTo;
                        if (lineTo) {
                            allCoords.add(`${lineTo[0]},${lineTo[1]}`);
                        }
                    }
                });
            });
        }

        // 3. Create variables for coordinates not already defined
        let coordCounter = 1;
        allCoords.forEach(coordKey => {
            if (!pointMap.has(coordKey)) {
                const [x, y] = coordKey.split(',').map(Number);
                const varName = `p${coordCounter++}`;
                pointMap.set(coordKey, varName);
                points.push({ name: varName, val: [x, y] });
            }
        });

        // Sort: predefined points (P*) first, then generated (p*)
        points.sort((a, b) => {
            const aIsPredefined = /^P[A-Z0-9]/i.test(a.name);
            const bIsPredefined = /^P[A-Z0-9]/i.test(b.name);
            if (aIsPredefined && !bIsPredefined) return -1;
            if (!aIsPredefined && bIsPredefined) return 1;
            return a.name.localeCompare(b.name, undefined, { numeric: true });
        });

        points.forEach(p => {
            generatedCode += `const ${p.name} = [${p.val[0]}, ${p.val[1]}];\n`;
        });
        generatedCode += '\n';

        // Helper to format point
        const formatPoint = (p: number[]) => {
            const key = `${p[0]},${p[1]}`;
            return pointMap.get(key) || `[${p[0]}, ${p[1]}]`;
        };

        // 4. Iterate list_region
        if (component.list_region) {
            component.list_region.forEach((region: any, rIdx: number) => {
                if (!Array.isArray(region)) return;
                if (rIdx > 0) generatedCode += `\n// Region ${rIdx + 1}\n`;

                let currentPen: { x: number, y: number } | null = null;

                region.forEach((segment: any) => {
                    // Extract move_to
                    let moveTo = segment.move_to;
                    if (moveTo && moveTo['py/tuple']) moveTo = moveTo['py/tuple'];
                    else if (Array.isArray(moveTo)) moveTo = moveTo;

                    // Check if we need to move
                    if (moveTo) {
                        const mx = moveTo[0];
                        const my = moveTo[1];

                        if (!currentPen || Math.abs(currentPen.x - mx) > 1e-6 || Math.abs(currentPen.y - my) > 1e-6) {
                            generatedCode += `MoveTo(${formatPoint(moveTo)});\n`;
                            currentPen = { x: mx, y: my };
                        }
                    }

                    if (segment.line_to) {
                        let lineTo = segment.line_to;
                        if (lineTo && lineTo['py/tuple']) lineTo = lineTo['py/tuple'];
                        else if (Array.isArray(lineTo)) lineTo = lineTo;

                        if (lineTo) {
                            generatedCode += `LineTo(${formatPoint(lineTo)});\n`;
                            currentPen = { x: lineTo[0], y: lineTo[1] };
                        }
                    } else if (segment.arc) {
                        let arc = segment.arc;
                        if (arc && arc['py/tuple']) arc = arc['py/tuple'];
                        else if (Array.isArray(arc)) arc = arc;

                        if (arc && moveTo) {
                            generatedCode += `ArcTo(${formatPoint(moveTo)}, ${arc[0]}, ${arc[1]}, ${arc[2]});\n`;

                            const r = arc[0];
                            const endAngle = arc[2];
                            const endX = moveTo[0] + r * Math.cos(endAngle);
                            const endY = moveTo[1] + r * Math.sin(endAngle);
                            currentPen = { x: endX, y: endY };
                        }
                    }
                });
            });
        }

        return generatedCode;
    };

    const handleComponentSelect = (value: string) => {
        setSelectedComponent(value);
        if (value === "custom") {
            setCode(DEFAULT_CODE);
        } else if (data && data[value as keyof GeometricComponentsObjects]) {
            const comp = data[value as keyof GeometricComponentsObjects];
            const newCode = generateCodeFromComponent(comp);
            setCode(newCode);
        }
    };

    const executeCode = useCallback(() => {
        setError(null);
        try {
            const generatedRegions: any[][] = [];
            const activeRegionSegments: any[] = [];
            let currentPen = [0, 0];
            const varNames: Record<string, string> = {}; // Map "x,y" -> "varName"

            // Extract variable definitions from code
            const varRegex = /const\s+(\w+)\s*=\s*\[\s*([\d.-]+)\s*,\s*([\d.-]+)\s*\]/g;
            let match;
            while ((match = varRegex.exec(code)) !== null) {
                const varName = match[1];
                const x = parseFloat(match[2]);
                const y = parseFloat(match[3]);
                varNames[`${x},${y}`] = varName;
            }
            setVariableNames(varNames);

            const dslMoveTo = (p: number[]) => {
                currentPen = p;
            };

            const dslLineTo = (p: number[]) => {
                activeRegionSegments.push({
                    move_to: { 'py/tuple': currentPen },
                    line_to: { 'py/tuple': p }
                });
                currentPen = p;
            };

            const dslArcTo = (center: number[], radius: number, startAngle: number, endAngle: number) => {
                // Calculate start point of the arc
                const startX = center[0] + radius * Math.cos(startAngle);
                const startY = center[1] + radius * Math.sin(startAngle);

                // If current pen is not at arc start, add a line to connect
                if (currentPen[0] !== startX || currentPen[1] !== startY) {
                    // Optional: Add line from current position to arc start
                    // For now, just move the pen without drawing
                }

                activeRegionSegments.push({
                    move_to: { 'py/tuple': center },
                    arc: { 'py/tuple': [radius, startAngle, endAngle] }
                });

                const endX = center[0] + radius * Math.cos(endAngle);
                const endY = center[1] + radius * Math.sin(endAngle);
                currentPen = [endX, endY];
            };

            const funcBody = `
                "use strict";
                ${code}
            `;

            const func = new Function('MoveTo', 'LineTo', 'ArcTo', funcBody);
            func(dslMoveTo, dslLineTo, dslArcTo);

            if (activeRegionSegments.length > 0) {
                generatedRegions.push(activeRegionSegments);
            }

            const mockGeometry: GeometricComponentsObjects = {
                // @ts-ignore
                customShape: {
                    name: selectedComponent === "custom" ? "Custom Shape" : selectedComponent,
                    color: "#4ade80",
                    list_region: generatedRegions,
                    variableNames: varNames // Pass variable names to geometry
                }
            };

            setGeometry(mockGeometry);

        } catch (err: any) {
            setError(err.message);
        }
    }, [code, selectedComponent]);

    useEffect(() => {
        executeCode();
    }, [executeCode]);

    return (
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-4 h-full">
            <Card className="flex flex-col">
                <CardHeader className="flex flex-row items-center justify-between py-2">
                    <div className="flex items-center gap-4">
                        <CardTitle>Script Editor</CardTitle>
                        <Select value={selectedComponent} onValueChange={handleComponentSelect}>
                            <SelectTrigger className="w-[180px]">
                                <SelectValue placeholder="Select Component" />
                            </SelectTrigger>
                            <SelectContent>
                                <SelectItem value="custom">Custom Playground</SelectItem>
                                {data && Object.entries(data).map(([key, component]) => (
                                    <SelectItem key={key} value={key}>
                                        {component.name || key}
                                    </SelectItem>
                                ))}
                            </SelectContent>
                        </Select>
                    </div>
                    <div className="flex gap-2">
                        <Button size="sm" variant="outline" onClick={() => setCode(DEFAULT_CODE)}>
                            <RotateCcw className="mr-2 h-4 w-4" /> Reset
                        </Button>
                        <Button size="sm" onClick={executeCode}>
                            <Play className="mr-2 h-4 w-4" /> Run
                        </Button>
                    </div>
                </CardHeader>
                <CardContent className="flex-1 p-0 overflow-hidden border-t">
                    <Editor
                        height="100%"
                        defaultLanguage="javascript"
                        value={code}
                        onChange={(value) => setCode(value || '')}
                        theme="vs-dark"
                        options={{
                            minimap: { enabled: false },
                            fontSize: 14,
                        }}
                    />
                </CardContent>
            </Card>

            <Card className="flex flex-col">
                <CardHeader>
                    <CardTitle>Preview</CardTitle>
                </CardHeader>
                <CardContent className="flex-1 bg-slate-50 dark:bg-slate-900 p-0 relative overflow-hidden">
                    {error && (
                        <div className="absolute top-0 left-0 right-0 bg-red-100 text-red-800 p-2 text-sm z-10 border-b border-red-200">
                            Error: {error}
                        </div>
                    )}
                    {geometry && <PlaygroundPreview geometry={geometry} variableNames={variableNames} />}
                </CardContent>
            </Card>
        </div>
    );
}

// Custom preview component for playground with variable name labels
function PlaygroundPreview({ geometry, variableNames }: { geometry: GeometricComponentsObjects, variableNames: Record<string, string> }) {
    const [scale, setScale] = useState(10);
    const [offset, setOffset] = useState({ x: 400, y: 300 });
    const [isDragging, setIsDragging] = useState(false);
    const [dragStart, setDragStart] = useState({ x: 0, y: 0 });
    const svgRef = useRef<SVGSVGElement>(null);
    const containerRef = useRef<HTMLDivElement>(null);

    const extractTuple = (obj: any): number[] | null => {
        if (!obj) return null;
        if (Array.isArray(obj)) return obj;
        if (obj['py/tuple']) return obj['py/tuple'];
        return null;
    };

    const renderComponent = (component: GeometricComponent) => {
        if (!component.list_region) return null;

        let pathData = "";
        const pointsToLabel = new Map<string, [number, number]>(); // Track unique points

        component.list_region.forEach(region => {
            if (!Array.isArray(region)) return;

            let currentPen: { x: number, y: number } | null = null;

            region.forEach((segment: any) => {
                const moveTo = extractTuple(segment.move_to);

                if (segment.line_to) {
                    const lineTo = extractTuple(segment.line_to);
                    if (moveTo && lineTo) {
                        const startX = moveTo[0];
                        const startY = -moveTo[1];
                        const endX = lineTo[0];
                        const endY = -lineTo[1];

                        if (!currentPen || Math.abs(currentPen.x - startX) > 1e-6 || Math.abs(currentPen.y - startY) > 1e-6) {
                            pathData += `M ${startX} ${startY} `;
                        }
                        pathData += `L ${endX} ${endY} `;
                        currentPen = { x: endX, y: endY };

                        // Track points for labeling
                        pointsToLabel.set(`${moveTo[0]},${moveTo[1]}`, [moveTo[0], moveTo[1]]);
                        pointsToLabel.set(`${lineTo[0]},${lineTo[1]}`, [lineTo[0], lineTo[1]]);
                    }
                } else if (segment.arc) {
                    const arcParams = extractTuple(segment.arc);
                    if (moveTo && arcParams && arcParams.length >= 3) {
                        const cx = moveTo[0];
                        const cy = -moveTo[1];
                        const r = arcParams[0];
                        const startAngle = arcParams[1];
                        const endAngle = arcParams[2];

                        const startX = cx + r * Math.cos(startAngle);
                        const startY = cy - r * Math.sin(startAngle);
                        const endX = cx + r * Math.cos(endAngle);
                        const endY = cy - r * Math.sin(endAngle);

                        if (!currentPen || Math.abs(currentPen.x - startX) > 1e-6 || Math.abs(currentPen.y - startY) > 1e-6) {
                            pathData += `M ${startX} ${startY} `;
                        }

                        let delta = endAngle - startAngle;
                        const largeArc = Math.abs(delta) > Math.PI ? 1 : 0;
                        const sweep = delta > 0 ? 0 : 1;

                        pathData += `A ${r} ${r} 0 ${largeArc} ${sweep} ${endX} ${endY} `;
                        currentPen = { x: endX, y: endY };

                        // Track arc end points
                        pointsToLabel.set(`${cx},${moveTo[1]}`, [cx, moveTo[1]]);
                    }
                }
            });
        });

        if (pathData) {
            pathData += 'Z';
        }

        return (
            <>
                {pathData && (
                    <path
                        d={pathData}
                        fill={component.color}
                        fillOpacity={0}
                        stroke="black"
                        strokeWidth={0.2 / scale}
                        className="transition-all duration-200"
                    />
                )}
                {/* Render point labels with variable names */}
                {Array.from(pointsToLabel.entries()).map(([key, [x, y]]) => {
                    const varName = variableNames[key];
                    if (!varName) return null;

                    return (
                        <g key={`label-${key}`}>
                            <circle
                                cx={x}
                                cy={-y}
                                r={2 / scale}
                                fill="#ef4444"
                                className="cursor-pointer"
                            />
                            <text
                                x={x + 3 / scale}
                                y={-y - 3 / scale}
                                fontSize={12 / scale}
                                fill="#1f2937"
                                className="font-mono font-bold"
                                style={{ pointerEvents: 'none' }}
                            >
                                {varName}
                            </text>
                        </g>
                    );
                })}
            </>
        );
    };

    // Wheel event handling with non-passive listener
    useEffect(() => {
        const container = containerRef.current;
        if (!container) return;

        const handleWheel = (e: WheelEvent) => {
            e.preventDefault();
            const zoomSensitivity = 0.001;
            setScale(s => Math.max(0.1, s * (1 - e.deltaY * zoomSensitivity)));
        };

        container.addEventListener('wheel', handleWheel, { passive: false });
        return () => {
            container.removeEventListener('wheel', handleWheel);
        };
    }, []);

    const handleMouseDown = (e: React.MouseEvent) => {
        setIsDragging(true);
        setDragStart({ x: e.clientX - offset.x, y: e.clientY - offset.y });
    };

    const handleMouseMove = (e: React.MouseEvent) => {
        if (isDragging) {
            setOffset({ x: e.clientX - dragStart.x, y: e.clientY - dragStart.y });
        }
    };


    const handleMouseUp = () => {
        setIsDragging(false);
    };

    // Auto-scale logic
    const handleAutoScale = useCallback(() => {
        if (!geometry || Object.keys(geometry).length === 0 || !svgRef.current) return;

        const viewWidth = svgRef.current.clientWidth;
        const viewHeight = svgRef.current.clientHeight;

        if (viewWidth === 0 || viewHeight === 0) return;

        let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
        let hasPoints = false;

        Object.values(geometry).forEach(comp => {
            const points = getPoints(comp);
            points.forEach(p => {
                hasPoints = true;
                minX = Math.min(minX, p[0]);
                maxX = Math.max(maxX, p[0]);
                minY = Math.min(minY, p[1]);
                maxY = Math.max(maxY, p[1]);
            });
        });

        if (!hasPoints) return;

        const padding = 0.1;
        const geomWidth = maxX - minX;
        const geomHeight = maxY - minY;
        const safeGeomWidth = geomWidth || 1;
        const safeGeomHeight = geomHeight || 1;

        const scaleX = viewWidth / (safeGeomWidth * (1 + padding));
        const scaleY = viewHeight / (safeGeomHeight * (1 + padding));
        const newScale = Math.min(scaleX, scaleY);

        const centerX = (minX + maxX) / 2;
        const centerY = (minY + maxY) / 2;

        const newOffsetX = viewWidth / 2 - centerX * newScale;
        const newOffsetY = viewHeight / 2 - (-centerY) * newScale;

        setScale(newScale);
        setOffset({ x: newOffsetX, y: newOffsetY });
    }, [geometry]);

    // Trigger auto-scale when geometry changes or container resizes
    useEffect(() => {
        handleAutoScale();

        const currentSvg = svgRef.current;
        if (!currentSvg) return;

        const resizeObserver = new ResizeObserver(() => {
            handleAutoScale();
        });

        resizeObserver.observe(currentSvg);

        return () => {
            resizeObserver.disconnect();
        };
    }, [handleAutoScale]);


    return (
        <div
            ref={containerRef}
            className="w-full h-full relative"
        >
            <div className="absolute top-2 left-2 bg-white/80 dark:bg-slate-800/80 p-2 rounded text-xs font-mono pointer-events-none z-10">
                Scale: {scale.toFixed(2)} px/mm
            </div>
            <div className="absolute top-2 right-2 flex gap-2 z-10">
                <Button variant="outline" size="icon" onClick={(e) => { e.preventDefault(); setScale((s: number) => s * 1.2); }}>
                    <ZoomIn className="h-4 w-4" />
                </Button>
                <Button variant="outline" size="icon" onClick={(e) => { e.preventDefault(); setScale((s: number) => s / 1.2); }}>
                    <ZoomOut className="h-4 w-4" />
                </Button>
                <Button variant="outline" size="icon" onClick={(e) => { e.preventDefault(); handleAutoScale(); }}>
                    <Maximize className="h-4 w-4" />
                </Button>
            </div>
            <svg
                ref={svgRef}
                className="w-full h-full cursor-move bg-slate-50 dark:bg-slate-900"
                onMouseDown={handleMouseDown}
                onMouseMove={handleMouseMove}
                onMouseUp={handleMouseUp}
                onMouseLeave={handleMouseUp}
            >
                <g transform={`translate(${offset.x}, ${offset.y}) scale(${scale})`}>
                    <line x1="-1000" y1="0" x2="1000" y2="0" stroke="#ddd" strokeWidth={1 / scale} />
                    <line x1="0" y1="-1000" x2="0" y2="1000" stroke="#ddd" strokeWidth={1 / scale} />
                    {Object.values(geometry).map((comp, idx) => (
                        <g key={idx}>{renderComponent(comp)}</g>
                    ))}
                </g>
            </svg>
        </div>
    );
}
