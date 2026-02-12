"use client"

import { useEffect, useState, useRef, useMemo } from 'react'
import axios from 'axios'
import * as d3 from 'd3'
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from '@/components/ui/card'
import { Loader2, RefreshCw, ZoomIn, ZoomOut, Maximize, MousePointer2 } from 'lucide-react'
import { Button } from '@/components/ui/button'
import { Badge } from '@/components/ui/badge'
import { cn } from "@/lib/utils"
import { ScrollArea } from '@/components/ui/scroll-area'
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from '@/components/ui/table'
import { Checkbox } from '@/components/ui/checkbox'
import { Label } from '@/components/ui/label'

interface Point {
    0: number; // x
    1: number; // y
}

interface PointData {
    HP: Record<string, [number, number]>;
    HP_mirror: Record<string, [number, number]>;
    RP: Record<string, [number, number]>;
    RP_mirror: Record<string, [number, number]>;
    parameters: Record<string, any>;
    num_slots: number;
    num_poles: number;
}

interface GeometrySegment {
    type: 'line' | 'arc';
    p1: [number, number];
    p2: [number, number];
    center?: [number, number];
}

interface GeometryRegion {
    innerCoord: [number, number];
    list_regions: GeometrySegment[][];
    mirrorAxis: [number, number] | null;
    bMirror: boolean;
    iRotateCopy: number;
    color: string;
}

interface GeometryData {
    regions: GeometryRegion[];
}

export default function DebugPointsPage() {
    const [data, setData] = useState<PointData | null>(null)
    const [geometry, setGeometry] = useState<GeometryData | null>(null)
    const [loading, setLoading] = useState(true)
    const [showGeometry, setShowGeometry] = useState(true)
    const [error, setError] = useState<string | null>(null)
    const svgRef = useRef<SVGSVGElement>(null)
    const gRef = useRef<SVGGElement>(null)
    const [zoomLevel, setZoomLevel] = useState(1)
    const [visibleTypes, setVisibleTypes] = useState<string[]>(['HP', 'RP', 'HP_M', 'RP_M'])

    const fetchData = async () => {
        setLoading(true)
        setError(null)
        try {
            const [pointsRes, geomRes] = await Promise.all([
                axios.get('http://localhost:8000/api/debug/points'),
                axios.get('http://localhost:8000/api/debug/geometry')
            ])
            setData(pointsRes.data)
            setGeometry(geomRes.data)
        } catch (err: any) {
            console.error('Error fetching debug data:', err)
            setError(err.message || 'Failed to fetch data')
        } finally {
            setLoading(false)
        }
    }

    useEffect(() => {
        fetchData()
    }, [])

    useEffect(() => {
        if (!data || !svgRef.current || !gRef.current) return

        const svg = d3.select(svgRef.current)
        const g = d3.select(gRef.current)

        // Clear previous elements
        g.selectAll('*').remove()

        const margin = 50
        const width = 800
        const height = 800

        // Helper functions for geometry
        const rotatePoint = (p: [number, number], deg: number): [number, number] => {
            const rad = (deg * Math.PI) / 180
            const cos = Math.cos(rad)
            const sin = Math.sin(rad)
            return [
                p[0] * cos - p[1] * sin,
                p[0] * sin + p[1] * cos
            ]
        }

        const mirrorPoint = (p: [number, number]): [number, number] => {
            return [p[0], -p[1]]
        }
        const allPoints: { x: number; y: number; name: string; type: string }[] = []

        Object.entries(data.HP).forEach(([id, coords]) => {
            if (visibleTypes.includes('HP')) allPoints.push({ x: coords[0], y: coords[1], name: `HP[${id}]`, type: 'HP' })
        })
        Object.entries(data.HP_mirror).forEach(([id, coords]) => {
            if (visibleTypes.includes('HP_M')) allPoints.push({ x: coords[0], y: coords[1], name: `HP[${id}]_M`, type: 'HP_M' })
        })
        Object.entries(data.RP).forEach(([id, coords]) => {
            if (visibleTypes.includes('RP')) allPoints.push({ x: coords[0], y: coords[1], name: `RP[${id}]`, type: 'RP' })
        })
        Object.entries(data.RP_mirror || {}).forEach(([id, coords]) => {
            if (visibleTypes.includes('RP_M')) allPoints.push({ x: coords[0], y: coords[1], name: `RP[${id}]_M`, type: 'RP_M' })
        })


        const xExtent = d3.extent(allPoints, d => d.x) as [number, number]
        const yExtent = d3.extent(allPoints, d => d.y) as [number, number]

        // Make sure we have a square-ish domain for proportional rendering
        const maxAbs = Math.max(
            Math.abs(xExtent[0] || 0), Math.abs(xExtent[1] || 0),
            Math.abs(yExtent[0] || 0), Math.abs(yExtent[1] || 0)
        ) * 1.2 || 10

        const xScale = d3.scaleLinear()
            .domain([-maxAbs, maxAbs])
            .range([0, width])

        const yScale = d3.scaleLinear()
            .domain([-maxAbs, maxAbs])
            .range([height, 0]) // Invert Y for Cartesian coordinates in SVG

        const rScale = (r: number) => r * (width / (2 * maxAbs))

        // Draw Geometry
        if (showGeometry && geometry) {
            geometry.regions.forEach((region, regionIdx) => {
                const copyCount = region.iRotateCopy || 1
                for (let i = 0; i < copyCount; i++) {
                    const rotationDeg = (i * 360) / copyCount

                    const drawSegments = (segments: GeometrySegment[], isMirrored: boolean) => {
                        segments.forEach((seg, segIdx) => {
                            let p1 = seg.p1
                            let p2 = seg.p2
                            let center = seg.center || [0, 0]

                            if (isMirrored) {
                                p1 = mirrorPoint(p1)
                                p2 = mirrorPoint(p2)
                                center = mirrorPoint(center)
                            }

                            p1 = rotatePoint(p1, rotationDeg)
                            p2 = rotatePoint(p2, rotationDeg)
                            center = rotatePoint(center, rotationDeg)

                            if (seg.type === 'line') {
                                g.append('line')
                                    .attr('x1', xScale(p1[0]))
                                    .attr('y1', yScale(p1[1]))
                                    .attr('x2', xScale(p2[0]))
                                    .attr('y2', yScale(p2[1]))
                                    .attr('stroke', region.color || '#444')
                                    .attr('stroke-width', 1.5)
                                    .attr('opacity', 0.6)
                            } else if (seg.type === 'arc') {
                                const r = Math.sqrt(Math.pow(p1[0] - center[0], 2) + Math.pow(p1[1] - center[1], 2))
                                const x1 = xScale(p1[0])
                                const y1 = yScale(p1[1])
                                const x2 = xScale(p2[0])
                                const y2 = yScale(p2[1])
                                const rx = rScale(r)

                                // Sweep flag logic: if mirrored, we might need to flip the sweep
                                // For now, let's keep it 1.
                                g.append('path')
                                    .attr('d', `M ${x1} ${y1} A ${rx} ${rx} 0 0 1 ${x2} ${y2}`)
                                    .attr('stroke', region.color || '#444')
                                    .attr('fill', 'none')
                                    .attr('stroke-width', 1.5)
                                    .attr('opacity', 0.6)
                            }
                        })
                    }

                    region.list_regions.forEach(segments => drawSegments(segments, false))
                    if (region.bMirror) {
                        region.list_regions.forEach(segments => drawSegments(segments, true))
                    }
                }
            })
        }

        // Add grid lines
        const gridLines = 10
        const step = maxAbs * 2 / gridLines
        for (let i = 0; i <= gridLines; i++) {
            const val = -maxAbs + i * step
            // Vertical
            g.append('line')
                .attr('x1', xScale(val))
                .attr('y1', 0)
                .attr('x2', xScale(val))
                .attr('y2', height)
                .attr('stroke', 'rgba(0,0,0,0.05)')
                .attr('stroke-width', 1)

            // Horizontal
            g.append('line')
                .attr('x1', 0)
                .attr('y1', yScale(val))
                .attr('x2', width)
                .attr('y2', yScale(val))
                .attr('stroke', 'rgba(0,0,0,0.05)')
                .attr('stroke-width', 1)
        }

        // Add Axes
        g.append('line')
            .attr('x1', 0)
            .attr('y1', yScale(0))
            .attr('x2', width)
            .attr('y2', yScale(0))
            .attr('stroke', 'rgba(0,0,0,0.3)')
            .attr('stroke-width', 2)

        g.append('line')
            .attr('x1', xScale(0))
            .attr('y1', 0)
            .attr('x2', xScale(0))
            .attr('y2', height)
            .attr('stroke', 'rgba(0,0,0,0.3)')
            .attr('stroke-width', 2)

        // Draw Points
        const pointGroups = g.selectAll('.point-group')
            .data(allPoints)
            .enter()
            .append('g')
            .attr('class', 'point-group')

        // Point Circle
        pointGroups.append('circle')
            .attr('cx', d => xScale(d.x))
            .attr('cy', d => yScale(d.y))
            .attr('r', 5)
            .attr('fill', d => {
                if (d.type === 'HP') return '#3b82f6' // Blue
                if (d.type === 'RP') return '#ef4444' // Red
                if (d.type === 'RP_M') return '#f59e0b' // Amber
                return '#10b981' // Green for HP Mirror
            })
            .attr('stroke', '#000')
            .attr('stroke-width', 1.5)
            .style('cursor', 'pointer')
            .on('mouseover', function () {
                d3.select(this).transition().duration(200).attr('r', 8)
            })
            .on('mouseout', function () {
                d3.select(this).transition().duration(200).attr('r', 5)
            })

        // Labels
        pointGroups.append('text')
            .attr('x', d => xScale(d.x) + 8)
            .attr('y', (d, i) => {
                // Alternate label position for clustered points (especially along axes)
                return yScale(d.y) + (i % 2 === 0 ? -12 : 18)
            })
            .attr('font-size', '12px')
            .attr('font-weight', '600')
            .attr('fill', '#1e293b')
            .style('pointer-events', 'none')
            .text(d => d.name)

        pointGroups.append('text')
            .attr('x', d => xScale(d.x) + 8)
            .attr('y', (d, i) => {
                return yScale(d.y) + (i % 2 === 0 ? 2 : 32)
            })
            .attr('font-size', '10px')
            .attr('fill', '#64748b')
            .style('pointer-events', 'none')
            .text(d => `(${d.x.toFixed(2)}, ${d.y.toFixed(2)})`)

        // Set up zoom
        const zoom = d3.zoom<SVGSVGElement, unknown>()
            .scaleExtent([0.1, 10])
            .on('zoom', (event) => {
                g.attr('transform', event.transform)
                setZoomLevel(event.transform.k)
            })

        svg.call(zoom as any)

        // Initial zoom to fit
        svg.call(zoom.transform as any, d3.zoomIdentity.translate(0, 0).scale(1))

    }, [data, geometry, visibleTypes, showGeometry])

    const handleResetZoom = () => {
        if (!svgRef.current) return
        d3.select(svgRef.current).call(d3.zoom().transform as any, d3.zoomIdentity)
    }

    return (
        <div className="flex-1 space-y-4 p-8 pt-6">
            <div className="flex items-center justify-between space-y-2">
                <div>
                    <h2 className="text-3xl font-bold tracking-tight">Geometry Debugger</h2>
                    <p className="text-muted-foreground">
                        Visualization of point coordinates from machine_designer_v2.py
                    </p>
                </div>
                <div className="flex items-center space-x-4">
                    <div className="flex items-center gap-4 bg-slate-50 p-1.5 px-3 rounded-full border border-slate-200">
                        {[
                            { id: 'HP', label: 'HP', color: 'bg-blue-500' },
                            { id: 'RP', label: 'RP', color: 'bg-red-500' },
                            { id: 'RP_M', label: 'RP_M', color: 'bg-amber-500' },
                            { id: 'HP_M', label: 'HP_M', color: 'bg-emerald-500' }
                        ].map(type => (
                            <label key={type.id} className="flex items-center gap-1.5 cursor-pointer hover:opacity-80 transition-opacity">
                                <Checkbox
                                    id={`toggle-${type.id}`}
                                    checked={visibleTypes.includes(type.id)}
                                    onCheckedChange={(checked) => {
                                        setVisibleTypes(prev =>
                                            checked
                                                ? [...prev, type.id]
                                                : prev.filter(t => t !== type.id)
                                        )
                                    }}
                                />
                                <div className={cn("w-2 h-2 rounded-full", type.color)} />
                                <span className="text-xs font-semibold text-slate-700">{type.label}</span>
                            </label>
                        ))}
                    </div>
                    <Button
                        onClick={() => setShowGeometry(!showGeometry)}
                        variant={showGeometry ? "default" : "outline"}
                        size="sm"
                    >
                        {showGeometry ? "Hide Geometry" : "Show Geometry"}
                    </Button>
                    <Button onClick={fetchData} disabled={loading} variant="outline" size="sm">
                        {loading ? <Loader2 className="mr-2 h-4 w-4 animate-spin" /> : <RefreshCw className="mr-2 h-4 w-4" />}
                        Refresh Data
                    </Button>
                </div>
            </div>

            <div className="grid gap-4 md:grid-cols-2 lg:grid-cols-7">
                <Card className="col-span-5 bg-white border-slate-200 shadow-sm overflow-hidden relative">
                    <div className="absolute top-4 right-4 z-10 flex flex-col space-y-2">
                        <Button variant="outline" size="icon" onClick={handleResetZoom} title="Reset Zoom">
                            <Maximize className="h-4 w-4" />
                        </Button>
                        <div className="bg-white/90 rounded-md px-2 py-1 text-[10px] text-slate-600 font-medium text-center border border-slate-200 shadow-sm">
                            {Math.round(zoomLevel * 100)}%
                        </div>
                    </div>
                    <CardContent className="p-0 flex items-center justify-center min-h-[800px]">
                        {loading ? (
                            <div className="flex flex-col items-center gap-2">
                                <Loader2 className="h-8 w-8 animate-spin text-primary" />
                                <p className="text-sm text-muted-foreground">Calculating points...</p>
                            </div>
                        ) : error ? (
                            <div className="text-center p-8">
                                <p className="text-destructive font-medium mb-2">Error</p>
                                <p className="text-sm text-muted-foreground mb-4">{error}</p>
                                <Button onClick={fetchData} variant="outline">Try Again</Button>
                            </div>
                        ) : (
                            <svg
                                ref={svgRef}
                                width="800"
                                height="800"
                                viewBox="0 0 800 800"
                                className="cursor-crosshair w-full h-auto max-h-[800px]"
                            >
                                <g ref={gRef} />
                            </svg>
                        )}
                    </CardContent>
                </Card>

                <Card className="col-span-2 border-slate-200 bg-white shadow-sm overflow-hidden flex flex-col h-[800px]">
                    <CardHeader className="pb-2">
                        <CardTitle className="text-lg">Points List</CardTitle>
                        <CardDescription>Coordinate breakdown by point name</CardDescription>
                    </CardHeader>
                    <CardContent className="flex-1 overflow-hidden p-0">
                        <ScrollArea className="h-full">
                            <Table>
                                <TableHeader className="sticky top-0 bg-white z-10">
                                    <TableRow className="sticky top-0 bg-white z-10">
                                        <TableHead className="w-[100px]">Name</TableHead>
                                        <TableHead className="text-right">X (mm)</TableHead>
                                        <TableHead className="text-right">Y (mm)</TableHead>
                                    </TableRow>
                                </TableHeader>
                                <TableBody>
                                    {(() => {
                                        const filteredSortedPoints = [
                                            ...Object.entries(data?.HP || {}).filter(() => visibleTypes.includes('HP')).map(([id, coords]) => ({ name: `HP[${id}]`, x: coords[0], y: coords[1], type: 'HP' })),
                                            ...Object.entries(data?.HP_mirror || {}).filter(() => visibleTypes.includes('HP_M')).map(([id, coords]) => ({ name: `HP[${id}]_M`, x: coords[0], y: coords[1], type: 'HP_M' })),
                                            ...Object.entries(data?.RP || {}).filter(() => visibleTypes.includes('RP')).map(([id, coords]) => ({ name: `RP[${id}]`, x: coords[0], y: coords[1], type: 'RP' })),
                                            ...Object.entries(data?.RP_mirror || {}).filter(() => visibleTypes.includes('RP_M')).map(([id, coords]) => ({ name: `RP[${id}]_M`, x: coords[0], y: coords[1], type: 'RP_M' })),
                                        ].sort((a, b) => a.name.localeCompare(b.name, undefined, { numeric: true }));

                                        return filteredSortedPoints.map((point) => (
                                            <TableRow key={point.name} className="hover:bg-slate-50 transition-colors">
                                                <TableCell className="font-medium py-2">
                                                    <div className="flex items-center gap-2">
                                                        <div className={cn("w-2 h-2 rounded-full",
                                                            point.type === 'HP' ? "bg-blue-500" :
                                                                point.type === 'RP' ? "bg-red-500" :
                                                                    point.type === 'RP_M' ? "bg-amber-500" : "bg-emerald-500"
                                                        )} />
                                                        {point.name}
                                                    </div>
                                                </TableCell>
                                                <TableCell className="text-right font-mono py-2">{point.x.toFixed(3)}</TableCell>
                                                <TableCell className="text-right font-mono py-2">{point.y.toFixed(3)}</TableCell>
                                            </TableRow>
                                        ));
                                    })()}
                                </TableBody>
                            </Table>
                        </ScrollArea>
                    </CardContent>
                </Card>
            </div>

            <div className="grid gap-4 md:grid-cols-2 lg:grid-cols-7 mt-4">
                <Card className="col-span-5 border-slate-200 bg-white shadow-sm">
                    <CardHeader className="py-4">
                        <CardTitle className="text-lg">Legend & Interaction</CardTitle>
                    </CardHeader>
                    <CardContent className="grid grid-cols-4 gap-8 pb-4">
                        <div className="space-y-2">
                            <div className="flex items-center gap-2">
                                <div className="w-3 h-3 rounded-full bg-blue-500" />
                                <span className="text-sm font-semibold">Horizontal (HP)</span>
                            </div>
                            <p className="text-xs text-muted-foreground">Main axis points.</p>
                        </div>
                        <div className="space-y-2">
                            <div className="flex items-center gap-2">
                                <div className="w-3 h-3 rounded-full bg-red-500" />
                                <span className="text-sm font-semibold">Rotated (RP)</span>
                            </div>
                            <p className="text-xs text-muted-foreground">Slot center points.</p>
                        </div>
                        <div className="space-y-2">
                            <div className="flex items-center gap-2">
                                <div className="w-3 h-3 rounded-full bg-amber-500" />
                                <span className="text-sm font-semibold">RP Mirror (RP_M)</span>
                            </div>
                            <p className="text-xs text-muted-foreground">Mirrored rotated.</p>
                        </div>
                        <div className="space-y-2">
                            <div className="flex items-center gap-2">
                                <div className="w-3 h-3 rounded-full bg-emerald-500" />
                                <span className="text-sm font-semibold">HP Mirror (HP_M)</span>
                            </div>
                            <p className="text-xs text-muted-foreground">Mirrored horizontal.</p>
                        </div>
                    </CardContent>
                </Card>

                <Card className="col-span-2 border-slate-200 bg-white shadow-sm">
                    <CardHeader className="py-4">
                        <CardTitle className="text-lg">Parameters</CardTitle>
                    </CardHeader>
                    <CardContent className="pb-4">
                        <div className="grid grid-cols-2 gap-2">
                            {Object.entries(data?.parameters || {}).slice(0, 4).map(([key, val]) => (
                                <div key={key} className="bg-slate-50 p-2 rounded border border-slate-100">
                                    <div className="text-[9px] text-slate-500 uppercase font-bold">{key}</div>
                                    <div className="text-xs font-mono">{typeof val === 'number' ? val.toFixed(2) : val}</div>
                                </div>
                            ))}
                        </div>
                    </CardContent>
                </Card>
            </div>
        </div>
    )
}
