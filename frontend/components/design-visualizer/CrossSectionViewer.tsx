'use client';

import React, { useState, useRef, useEffect, useCallback } from 'react';
import { GeometricComponentsObjects, GeometricComponent, getPoints, GP } from '@/lib/DesignData';
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { ZoomIn, ZoomOut, Maximize, Settings } from "lucide-react";
import {
    Tooltip,
    TooltipContent,
    TooltipProvider,
    TooltipTrigger,
} from "@/components/ui/tooltip";
import { Checkbox } from "@/components/ui/checkbox";
import {
    Popover,
    PopoverContent,
    PopoverTrigger,
} from "@/components/ui/popover";

interface CrossSectionViewerProps {
    geometry: GeometricComponentsObjects;
    selectedComponent?: string | null;
    visibility?: Record<string, boolean>;
    onVisibilityChange?: (visibility: Record<string, boolean>) => void;
    showParameters?: boolean;
    gpData?: GP;
}

const CrossSectionViewer: React.FC<CrossSectionViewerProps> = ({
    geometry,
    selectedComponent,
    visibility: externalVisibility,
    onVisibilityChange,
    showParameters = false,
    gpData
}) => {
    const [scale, setScale] = useState(10);
    const [offset, setOffset] = useState({ x: 400, y: 300 });
    const [isDragging, setIsDragging] = useState(false);
    const [dragStart, setDragStart] = useState({ x: 0, y: 0 });
    const [internalVisibility, setInternalVisibility] = useState<Record<string, boolean>>({});
    const [isMounted, setIsMounted] = useState(false);
    const svgRef = useRef<SVGSVGElement>(null);

    // Track mounting to prevent hydration mismatch
    useEffect(() => {
        setIsMounted(true);
    }, []);

    // Use external visibility if provided, otherwise use internal
    const visibility = externalVisibility || internalVisibility;
    const setVisibility = onVisibilityChange || setInternalVisibility;

    // Auto-scale logic
    const handleAutoScale = useCallback(() => {
        if (!geometry || Object.keys(geometry).length === 0 || !svgRef.current) return;

        const viewWidth = svgRef.current.clientWidth;
        const viewHeight = svgRef.current.clientHeight;

        if (viewWidth === 0 || viewHeight === 0) return;

        let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
        let hasPoints = false;

        Object.values(geometry).forEach(comp => {
            if (!comp) return;
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

    // Initialize internal visibility state when geometry changes (only if not controlled by parent)
    useEffect(() => {
        if (externalVisibility || !geometry) return;
        const initialVisibility: Record<string, boolean> = {};
        Object.keys(geometry).forEach(key => {
            // Default all to visible except sleeve
            initialVisibility[key] = key.toLowerCase().includes('sleeve') ? false : true;
        });
        setInternalVisibility(initialVisibility);
    }, [geometry, externalVisibility]);

    if (!geometry) {
        return (
            <div className="h-full flex items-center justify-center p-4">
                <div className="text-center space-y-2">
                    <p className="text-red-500 font-medium">无几何数据可用</p>
                    <p className="text-sm text-muted-foreground">
                        GeometricComponentsObjects 为 null 或未定义
                    </p>
                </div>
            </div>
        );
    }

    // Check if geometry has any valid components
    const hasValidComponents = Object.values(geometry).some(comp => 
        comp !== null && comp !== undefined && 
        (comp.list_region || getPoints(comp).length > 0)
    );

    if (!hasValidComponents) {
        return (
            <div className="h-full flex items-center justify-center p-4">
                <div className="text-center space-y-2">
                    <p className="text-yellow-600 dark:text-yellow-400 font-medium">几何数据为空</p>
                    <p className="text-sm text-muted-foreground">
                        所有几何组件均为 null 或没有有效的点数据
                    </p>
                </div>
            </div>
        );
    }

    // Helper to extract tuple values from python pickle JSON format
    const extractTuple = (obj: any): number[] | null => {
        if (!obj) return null;
        if (Array.isArray(obj)) return obj;
        if (obj['py/tuple']) return obj['py/tuple'];
        return null;
    };

    const renderPathFromRegions = (component: GeometricComponent) => {
        if (!component.list_region) return null;

        // Check if list_region is empty or contains only empty/invalid regions
        const hasValidRegions = component.list_region.some(region =>
            Array.isArray(region) && region.length > 0
        );

        if (!hasValidRegions) return null;

        let pathData = "";

        component.list_region.forEach(region => {
            // Skip if region is not an array (e.g., {"py/id": 232} references)
            if (!Array.isArray(region)) return;

            // Skip empty regions
            if (region.length === 0) return;

            let currentPen: { x: number, y: number } | null = null;
            let regionPath = ""; // Separate path for each region to close individually

            region.forEach((segment: any) => {
                const moveTo = extractTuple(segment.move_to);

                if (segment.line_to) {
                    const lineTo = extractTuple(segment.line_to);
                    if (moveTo && lineTo) {
                        const startX = moveTo[0];
                        const startY = -moveTo[1]; // Flip Y for SVG
                        const endX = lineTo[0];
                        const endY = -lineTo[1];

                        if (!currentPen || Math.abs(currentPen.x - startX) > 1e-6 || Math.abs(currentPen.y - startY) > 1e-6) {
                            regionPath += `M ${startX} ${startY} `;
                        }
                        regionPath += `L ${endX} ${endY} `;
                        currentPen = { x: endX, y: endY };
                    }
                } else if (segment.arc) {
                    const arcParams = extractTuple(segment.arc);
                    if (moveTo && arcParams && arcParams.length >= 3) {
                        const cx = moveTo[0];
                        const cy = -moveTo[1]; // Flip Y
                        const r = arcParams[0];
                        const startAngle = arcParams[1];
                        const endAngle = arcParams[2];

                        const startX = cx + r * Math.cos(startAngle);
                        const startY = cy - r * Math.sin(startAngle);
                        const endX = cx + r * Math.cos(endAngle);
                        const endY = cy - r * Math.sin(endAngle);

                        if (!currentPen || Math.abs(currentPen.x - startX) > 1e-6 || Math.abs(currentPen.y - startY) > 1e-6) {
                            regionPath += `M ${startX} ${startY} `;
                        }

                        let delta = endAngle - startAngle;
                        const largeArc = Math.abs(delta) > Math.PI ? 1 : 0;
                        const sweep = delta > 0 ? 0 : 1;

                        regionPath += `A ${r} ${r} 0 ${largeArc} ${sweep} ${endX} ${endY} `;
                        currentPen = { x: endX, y: endY };
                    }
                }
            });

            // Close each region individually and add to main path
            if (regionPath) {
                pathData += regionPath + 'Z ';
            }
        });

        // Return null if no valid path data was generated
        if (!pathData) return null;

        return pathData;
    };

    const renderComponent = (component: GeometricComponent | null, componentKey: string) => {
        // Skip if component is null or hidden
        if (!component || visibility[componentKey] === false) return null;

        const isSelected = selectedComponent === componentKey;
        const isDimmed = selectedComponent && !isSelected;
        const fillOpacity = 0; // No fill, only outlines
        const strokeOpacity = isDimmed ? 0.3 : 1;
        const strokeWidth = isSelected ? 4.0 / scale : 1.5 / scale; // Thicker stroke for selected component

        // Prefer list_region if available
        if (component.list_region) {
            const d = renderPathFromRegions(component);
            if (d) {
                return (
                    <path
                        key={component.name}
                        d={d}
                        fill={component.color}
                        fillOpacity={fillOpacity}
                        stroke="black"
                        strokeWidth={strokeWidth}
                        strokeOpacity={strokeOpacity}
                        className="transition-all duration-200"
                    />
                );
            }
        }

        // Fallback to points
        const points = getPoints(component);
        if (points.length < 2) return null;
        const pathData = points.map((p, i) => {
            const x = p[0];
            const y = -p[1]; // Flip Y
            return `${i === 0 ? 'M' : 'L'} ${x} ${y}`;
        }).join(' ') + ' Z';
        return (
            <path
                key={component.name}
                d={pathData}
                fill={component.color}
                fillOpacity={fillOpacity}
                stroke="black"
                strokeWidth={strokeWidth}
                strokeOpacity={strokeOpacity}
                className="transition-all duration-200"
            />
        );
    };



    const renderPoints = (component: GeometricComponent | null) => {
        if (!component) return null;
        const points = getPoints(component);
        return points.map((p, i) => (
            <TooltipProvider key={`${component.name}-p${i}`}>
                <Tooltip>
                    <TooltipTrigger asChild>
                        <circle
                            cx={p[0]}
                            cy={-p[1]}
                            r={2 / scale}
                            fill="red"
                            className="cursor-pointer hover:fill-yellow-400"
                        />
                    </TooltipTrigger>
                    <TooltipContent>
                        <p className="font-mono text-xs">
                            {`P${i + 1}: (${p[0].toFixed(2)}, ${p[1].toFixed(2)})`}
                        </p>
                    </TooltipContent>
                </Tooltip>
            </TooltipProvider>
        ));
    };

    const renderParameterAnnotations = () => {
        if (!gpData) return null;

        const annotations: JSX.Element[] = [];
        const fontSize = 10 / scale;
        const strokeWidth = 0.5 / scale;

        // Get key parameters
        const r_so = gpData.mm_r_so?.value ?? 0;
        const r_si = gpData.mm_r_si?.value ?? 0;
        const r_ro = gpData.mm_r_ro?.value ?? 0;
        const r_ri = gpData.mm_r_ri?.value ?? 0;
        const d_pm = gpData.mm_d_pm?.value ?? 0;
        const d_mech_air_gap = gpData.mm_d_mech_air_gap?.value ?? 0;

        // Annotate stator outer radius (r_so)
        if (r_so > 0) {
            const angle = Math.PI / 6;
            const x1 = r_so * Math.cos(angle);
            const y1 = -r_so * Math.sin(angle);
            const x2 = (r_so + 10) * Math.cos(angle);
            const y2 = -(r_so + 10) * Math.sin(angle);
            annotations.push(
                <g key="r_so">
                    <line x1={x1} y1={y1} x2={x2} y2={y2} stroke="#3b82f6" strokeWidth={strokeWidth} />
                    <line x1={x2 - 3} y1={y2 - 3} x2={x2} y2={y2} stroke="#3b82f6" strokeWidth={strokeWidth} />
                    <line x1={x2 - 3} y1={y2 + 3} x2={x2} y2={y2} stroke="#3b82f6" strokeWidth={strokeWidth} />
                    <text
                        x={x2 + 5}
                        y={y2}
                        fill="#3b82f6"
                        fontSize={fontSize}
                        textAnchor="start"
                        dominantBaseline="middle"
                    >
                        r_so = {r_so.toFixed(1)}mm
                    </text>
                </g>
            );
        }

        // Annotate stator inner radius (r_si)
        if (r_si > 0) {
            const angle = Math.PI / 4;
            const x1 = r_si * Math.cos(angle);
            const y1 = -r_si * Math.sin(angle);
            const x2 = (r_si + 8) * Math.cos(angle);
            const y2 = -(r_si + 8) * Math.sin(angle);
            annotations.push(
                <g key="r_si">
                    <line x1={x1} y1={y1} x2={x2} y2={y2} stroke="#10b981" strokeWidth={strokeWidth} />
                    <text
                        x={x2 + 5}
                        y={y2}
                        fill="#10b981"
                        fontSize={fontSize}
                        textAnchor="start"
                        dominantBaseline="middle"
                    >
                        r_si = {r_si.toFixed(1)}mm
                    </text>
                </g>
            );
        }

        // Annotate rotor outer radius (r_ro)
        if (r_ro > 0) {
            const angle = -Math.PI / 3;
            const x1 = r_ro * Math.cos(angle);
            const y1 = -r_ro * Math.sin(angle);
            const x2 = (r_ro - 8) * Math.cos(angle);
            const y2 = -(r_ro - 8) * Math.sin(angle);
            annotations.push(
                <g key="r_ro">
                    <line x1={x1} y1={y1} x2={x2} y2={y2} stroke="#f59e0b" strokeWidth={strokeWidth} />
                    <text
                        x={x2 - 5}
                        y={y2}
                        fill="#f59e0b"
                        fontSize={fontSize}
                        textAnchor="end"
                        dominantBaseline="middle"
                    >
                        r_ro = {r_ro.toFixed(1)}mm
                    </text>
                </g>
            );
        }

        // Annotate rotor inner radius (r_ri)
        if (r_ri > 0) {
            const angle = -Math.PI / 5;
            const x1 = r_ri * Math.cos(angle);
            const y1 = -r_ri * Math.sin(angle);
            const x2 = (r_ri - 6) * Math.cos(angle);
            const y2 = -(r_ri - 6) * Math.sin(angle);
            annotations.push(
                <g key="r_ri">
                    <line x1={x1} y1={y1} x2={x2} y2={y2} stroke="#8b5cf6" strokeWidth={strokeWidth} />
                    <text
                        x={x2 - 5}
                        y={y2}
                        fill="#8b5cf6"
                        fontSize={fontSize}
                        textAnchor="end"
                        dominantBaseline="middle"
                    >
                        r_ri = {r_ri.toFixed(1)}mm
                    </text>
                </g>
            );
        }

        // Annotate air gap
        if (d_mech_air_gap > 0 && r_ro > 0 && r_si > 0) {
            const angle = Math.PI / 2;
            const r_mid = (r_ro + r_si) / 2;
            const x = r_mid * Math.cos(angle);
            const y = -r_mid * Math.sin(angle);
            annotations.push(
                <g key="air_gap">
                    <line x1={r_ro * Math.cos(angle)} y1={-r_ro * Math.sin(angle)} 
                          x2={r_si * Math.cos(angle)} y2={-r_si * Math.sin(angle)} 
                          stroke="#ef4444" strokeWidth={strokeWidth} strokeDasharray={`${2/scale} ${2/scale}`} />
                    <text
                        x={x}
                        y={y - 8}
                        fill="#ef4444"
                        fontSize={fontSize}
                        textAnchor="middle"
                        dominantBaseline="middle"
                    >
                        gap = {d_mech_air_gap.toFixed(2)}mm
                    </text>
                </g>
            );
        }

        // Annotate magnet thickness
        if (d_pm > 0 && r_ro > 0) {
            const angle = -Math.PI / 6;
            const r_mid = r_ro - d_pm / 2;
            const x1 = (r_ro - d_pm) * Math.cos(angle);
            const y1 = -(r_ro - d_pm) * Math.sin(angle);
            const x2 = r_ro * Math.cos(angle);
            const y2 = -r_ro * Math.sin(angle);
            annotations.push(
                <g key="d_pm">
                    <line x1={x1} y1={y1} x2={x2} y2={y2} stroke="#ec4899" strokeWidth={strokeWidth} />
                    <line x1={x1 - 2} y1={y1 - 2} x2={x1} y2={y1} stroke="#ec4899" strokeWidth={strokeWidth} />
                    <line x1={x1 - 2} y1={y1 + 2} x2={x1} y2={y1} stroke="#ec4899" strokeWidth={strokeWidth} />
                    <line x1={x2 + 2} y1={y2 - 2} x2={x2} y2={y2} stroke="#ec4899" strokeWidth={strokeWidth} />
                    <line x1={x2 + 2} y1={y2 + 2} x2={x2} y2={y2} stroke="#ec4899" strokeWidth={strokeWidth} />
                    <text
                        x={r_mid * Math.cos(angle) + 5}
                        y={-r_mid * Math.sin(angle)}
                        fill="#ec4899"
                        fontSize={fontSize}
                        textAnchor="start"
                        dominantBaseline="middle"
                    >
                        d_pm = {d_pm.toFixed(1)}mm
                    </text>
                </g>
            );
        }

        return annotations;
    };

    const handleWheel = (e: React.WheelEvent) => {
        e.preventDefault();
        const zoomSensitivity = 0.001;
        const newScale = Math.max(0.1, scale * (1 - e.deltaY * zoomSensitivity));
        setScale(newScale);
    };

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

    const Ruler = () => {
        const targetPx = 100;
        const rawMm = targetPx / scale;

        let niceMm = 1;
        const bases = [1, 2, 5];
        let power = 0;

        while (true) {
            const m = Math.pow(10, power);
            let found = false;
            for (const base of bases) {
                const val = base * m;
                if (val * scale >= 50) {
                    niceMm = val;
                    found = true;
                    break;
                }
            }
            if (found) break;
            power++;
            if (power > 10) break;
        }

        const widthPx = niceMm * scale;

        return (
            <div className="absolute bottom-4 right-4 bg-white/80 dark:bg-slate-800/80 p-2 rounded border shadow-sm pointer-events-none z-10">
                <div className="flex flex-col items-center gap-1">
                    <div className="flex flex-col items-center">
                        <div className="h-2 border-l border-r border-b border-black dark:border-white w-full" style={{ width: widthPx }}></div>
                    </div>
                    <span className="text-xs font-mono text-black dark:text-white">{niceMm} mm</span>
                </div>
            </div>
        );
    };

    return (
        <Card className="h-full flex flex-col">
            <CardHeader className="flex flex-row items-center justify-between py-2">
                <CardTitle className="text-lg">Cross Section</CardTitle>
                <div className="flex gap-2">
                    <Popover>
                        <PopoverTrigger asChild>
                            <Button variant="outline" size="icon">
                                <Settings className="h-4 w-4" />
                            </Button>
                        </PopoverTrigger>
                        <PopoverContent className="w-64">
                            <div className="space-y-2">
                                <h4 className="font-medium text-sm mb-3">Component Visibility</h4>
                                {Object.entries(geometry).map(([key, component]) => (
                                    <div key={key} className="flex items-center space-x-2">
                                        <Checkbox
                                            id={`visibility-${key}`}
                                            checked={visibility[key] !== false}
                                            onCheckedChange={(checked: boolean) => {
                                                const newVisibility = {
                                                    ...visibility,
                                                    [key]: checked === true
                                                };
                                                setVisibility(newVisibility);
                                            }}
                                        />
                                        <label
                                            htmlFor={`visibility-${key}`}
                                            className="text-sm font-normal leading-none peer-disabled:cursor-not-allowed peer-disabled:opacity-70 cursor-pointer"
                                        >
                                            {component?.name || key}
                                        </label>
                                    </div>
                                ))}
                            </div>
                        </PopoverContent>
                    </Popover>
                    <Button variant="outline" size="icon" onClick={() => setScale(s => s * 1.2)}>
                        <ZoomIn className="h-4 w-4" />
                    </Button>
                    <Button variant="outline" size="icon" onClick={() => setScale(s => s / 1.2)}>
                        <ZoomOut className="h-4 w-4" />
                    </Button>
                    <Button variant="outline" size="icon" onClick={handleAutoScale}>
                        <Maximize className="h-4 w-4" />
                    </Button>
                </div>
            </CardHeader>
            <CardContent className="flex-1 p-0 overflow-hidden relative">
                <div className="absolute top-2 left-2 bg-white/80 p-2 rounded text-xs font-mono pointer-events-none z-10">
                    Scale: {scale.toFixed(2)} px/mm
                </div>
                <Ruler />
                <svg
                    ref={svgRef}
                    className="w-full h-full cursor-move bg-slate-50"
                    onWheel={handleWheel}
                    onMouseDown={handleMouseDown}
                    onMouseMove={handleMouseMove}
                    onMouseUp={handleMouseUp}
                    onMouseLeave={handleMouseUp}
                >
                    <g transform={`translate(${offset.x}, ${offset.y}) scale(${scale})`}>
                        {/* Grid lines */}
                        <line x1="-1000" y1="0" x2="1000" y2="0" stroke="#ddd" strokeWidth={1 / scale} />
                        <line x1="0" y1="-1000" x2="0" y2="1000" stroke="#ddd" strokeWidth={1 / scale} />

                        {/* Components - only render after mount to prevent hydration mismatch */}
                        {isMounted && Object.entries(geometry).map(([key, comp]) => comp ? renderComponent(comp, key) : null)}

                        {/* Points Overlay */}
                        {isMounted && Object.values(geometry).map((comp, idx) => renderPoints(comp))}

                        {/* Parameter Annotations */}
                        {isMounted && showParameters && gpData && renderParameterAnnotations()}
                    </g>
                </svg>
            </CardContent>
        </Card>
    );
};

export default CrossSectionViewer;
