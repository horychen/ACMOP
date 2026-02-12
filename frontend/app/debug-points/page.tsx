"use client"

import { useEffect, useState, useRef, useMemo } from 'react'
import axios from 'axios'
import * as d3 from 'd3'
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from '@/components/ui/card'
import { Loader2, RefreshCw, ZoomIn, ZoomOut, Maximize, MousePointer2 } from 'lucide-react'
import { Button } from '@/components/ui/button'
import { Badge } from '@/components/ui/badge'
import { Separator } from '@/components/ui/separator'

interface Point {
    0: number; // x
    1: number; // y
}

interface PointData {
    HP: Record<string, [number, number]>;
    HP_mirror: Record<string, [number, number]>;
    RP: Record<string, [number, number]>;
    parameters: Record<string, any>;
}

export default function DebugPointsPage() {
    const [data, setData] = useState<PointData | null>(null)
    const [loading, setLoading] = useState(true)
    const [error, setError] = useState<string | null>(null)
    const svgRef = useRef<SVGSVGElement>(null)
    const gRef = useRef<SVGGElement>(null)
    const [zoomLevel, setZoomLevel] = useState(1)

    const fetchData = async () => {
        setLoading(true)
        setError(null)
        try {
            const response = await axios.get('http://localhost:8000/api/debug/points')
            setData(response.data)
        } catch (err: any) {
            console.error('Error fetching debug points:', err)
            setError(err.message || 'Failed to fetch debug points')
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

        // Flatten points for extent calculation
        const allPoints: { x: number; y: number; name: string; type: string }[] = []

        Object.entries(data.HP).forEach(([id, coords]) => {
            allPoints.push({ x: coords[0], y: coords[1], name: `HP[${id}]`, type: 'HP' })
        })
        Object.entries(data.HP_mirror).forEach(([id, coords]) => {
            allPoints.push({ x: coords[0], y: coords[1], name: `HP[${id}]_M`, type: 'HP_M' })
        })
        Object.entries(data.RP).forEach(([id, coords]) => {
            allPoints.push({ x: coords[0], y: coords[1], name: `RP[${id}]`, type: 'RP' })
        })

        const margin = 50
        const width = 800
        const height = 800

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
                .attr('stroke', 'rgba(255,255,255,0.05)')
                .attr('stroke-width', 1)

            // Horizontal
            g.append('line')
                .attr('x1', 0)
                .attr('y1', yScale(val))
                .attr('x2', width)
                .attr('y2', yScale(val))
                .attr('stroke', 'rgba(255,255,255,0.05)')
                .attr('stroke-width', 1)
        }

        // Add Axes
        g.append('line')
            .attr('x1', 0)
            .attr('y1', yScale(0))
            .attr('x2', width)
            .attr('y2', yScale(0))
            .attr('stroke', 'rgba(255,255,255,0.2)')
            .attr('stroke-width', 2)

        g.append('line')
            .attr('x1', xScale(0))
            .attr('y1', 0)
            .attr('x2', xScale(0))
            .attr('y2', height)
            .attr('stroke', 'rgba(255,255,255,0.2)')
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
                return '#10b981' // Green for Mirror
            })
            .attr('stroke', '#fff')
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
            .attr('y', d => yScale(d.y) + 8)
            .attr('font-size', '12px')
            .attr('font-weight', '500')
            .attr('fill', '#e2e8f0')
            .text(d => d.name)

        pointGroups.append('text')
            .attr('x', d => xScale(d.x) + 8)
            .attr('y', d => yScale(d.y) + 22)
            .attr('font-size', '10px')
            .attr('fill', '#94a3b8')
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

    }, [data])

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
                <div className="flex items-center space-x-2">
                    <Button onClick={fetchData} disabled={loading} variant="outline" size="sm">
                        {loading ? <Loader2 className="mr-2 h-4 w-4 animate-spin" /> : <RefreshCw className="mr-2 h-4 w-4" />}
                        Refresh Data
                    </Button>
                </div>
            </div>

            <div className="grid gap-4 md:grid-cols-2 lg:grid-cols-7">
                <Card className="col-span-5 bg-black/40 backdrop-blur-sm border-slate-800 overflow-hidden relative">
                    <div className="absolute top-4 right-4 z-10 flex flex-col space-y-2">
                        <Button variant="secondary" size="icon" onClick={handleResetZoom} title="Reset Zoom">
                            <Maximize className="h-4 w-4" />
                        </Button>
                        <div className="bg-slate-900/80 rounded-md px-2 py-1 text-[10px] text-slate-400 text-center border border-slate-700">
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

                <Card className="col-span-2 border-slate-800 bg-card/50">
                    <CardHeader>
                        <CardTitle className="text-lg">Point Legend</CardTitle>
                        <CardDescription>Types of points defined in the design</CardDescription>
                    </CardHeader>
                    <CardContent className="space-y-6">
                        <div className="space-y-4">
                            <div className="flex items-center justify-between">
                                <div className="flex items-center gap-2">
                                    <div className="w-3 h-3 rounded-full bg-blue-500 shadow-[0_0_8px_rgba(59,130,246,0.5)]" />
                                    <span className="text-sm font-medium">Horizontal (HP)</span>
                                </div>
                                <Badge variant="outline" className="text-[10px] px-1 py-0">{Object.keys(data?.HP || {}).length}</Badge>
                            </div>
                            <p className="text-xs text-muted-foreground pl-5">
                                Points defined on the horizontal axis (Tooth Axis).
                            </p>

                            <Separator className="bg-slate-800" />

                            <div className="flex items-center justify-between">
                                <div className="flex items-center gap-2">
                                    <div className="w-3 h-3 rounded-full bg-red-500 shadow-[0_0_8px_rgba(239,68,68,0.5)]" />
                                    <span className="text-sm font-medium">Rotated (RP)</span>
                                </div>
                                <Badge variant="outline" className="text-[10px] px-1 py-0">{Object.keys(data?.RP || {}).length}</Badge>
                            </div>
                            <p className="text-xs text-muted-foreground pl-5">
                                Points rotated at half slot pitch (Slot Center Axis).
                            </p>

                            <Separator className="bg-slate-800" />

                            <div className="flex items-center justify-between">
                                <div className="flex items-center gap-2">
                                    <div className="w-3 h-3 rounded-full bg-emerald-500 shadow-[0_0_8px_rgba(16,185,129,0.5)]" />
                                    <span className="text-sm font-medium">Mirror (HP_M)</span>
                                </div>
                                <Badge variant="outline" className="text-[10px] px-1 py-0">{Object.keys(data?.HP_mirror || {}).length}</Badge>
                            </div>
                            <p className="text-xs text-muted-foreground pl-5">
                                Coordinates mirrored for symmetric part generation.
                            </p>
                        </div>

                        <div className="mt-8 pt-6 border-t border-slate-800">
                            <h4 className="text-sm font-medium mb-4 flex items-center gap-2">
                                <MousePointer2 className="h-3.5 w-3.5" />
                                Interaction Tips
                            </h4>
                            <ul className="text-xs text-muted-foreground space-y-2 list-disc pl-4">
                                <li>Use mouse wheel to <b>zoom</b> in/out on points.</li>
                                <li>Click and drag to <b>pan</b> the coordinate system.</li>
                                <li>Hover over a point to highlight its location.</li>
                                <li>Scale is in <b>mm</b>. Center is (0,0).</li>
                            </ul>
                        </div>

                        {data?.parameters && (
                            <div className="mt-8 pt-6 border-t border-slate-800">
                                <h4 className="text-sm font-medium mb-3">Active Parameters</h4>
                                <div className="grid grid-cols-2 gap-2">
                                    {Object.entries(data.parameters).slice(0, 6).map(([key, val]) => (
                                        <div key={key} className="bg-slate-900/50 p-2 rounded border border-slate-800">
                                            <div className="text-[10px] text-slate-500 uppercase tracking-wider">{key.replace('r_', 'Radius ').replace('d_', 'Depth ').replace('w_', 'Width ')}</div>
                                            <div className="text-xs font-mono font-semibold">{typeof val === 'number' ? val.toFixed(2) : val}</div>
                                        </div>
                                    ))}
                                </div>
                            </div>
                        )}
                    </CardContent>
                </Card>
            </div>
        </div>
    )
}
