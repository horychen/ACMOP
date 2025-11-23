"use client";

import React, { useEffect, useRef } from 'react';
import * as d3 from 'd3';
import { MachineGeometry } from '../types';

interface MotorVisualizerProps {
    geometry: MachineGeometry;
}

const MotorVisualizer: React.FC<MotorVisualizerProps> = ({ geometry }) => {
    const svgRef = useRef<SVGSVGElement>(null);

    useEffect(() => {
        if (!svgRef.current) return;

        const svg = d3.select(svgRef.current);
        svg.selectAll("*").remove(); // Clear previous

        const width = 400;
        const height = 400;
        const centerX = width / 2;
        const centerY = height / 2;
        const scale = Math.min(width, height) / (2.2 * geometry.statorOuterRadius);

        const g = svg.append("g")
            .attr("transform", `translate(${centerX}, ${centerY}) scale(${scale})`);

        // --- Draw Stator ---

        // Stator Outer Circle (Yoke)
        g.append("circle")
            .attr("r", geometry.statorOuterRadius)
            .attr("fill", "#334155") // Slate 700
            .attr("stroke", "#94a3b8")
            .attr("stroke-width", 1 / scale);

        // Stator Slots & Teeth
        // We draw the negative space (slots) to verify geometry or draw teeth as paths
        const toothAngle = (2 * Math.PI) / geometry.slots;

        // Helper to polar -> cartesian
        const pol2car = (r: number, a: number) => ({ x: r * Math.cos(a), y: r * Math.sin(a) });

        // Inner limit of stator
        g.append("circle")
            .attr("r", geometry.statorInnerRadius)
            .attr("fill", "#0f172a") // Background color (hollow)
            .attr("stroke", "none");

        // Draw Teeth
        for (let i = 0; i < geometry.slots; i++) {
            const theta = i * toothAngle;

            g.append("rect")
                .attr("x", 0)
                .attr("y", -geometry.toothWidth / 2)
                .attr("width", geometry.statorOuterRadius - geometry.statorInnerRadius) // Approximation
                .attr("height", geometry.toothWidth)
                .attr("fill", "#334155")
                .attr("stroke", "#475569")
                .attr("stroke-width", 0.5 / scale)
                .attr("transform", `rotate(${(theta * 180 / Math.PI)}, 0, 0) translate(${geometry.statorInnerRadius}, 0)`);

            // Copper Windings (in slots)
            // Draw circles in the gaps between teeth
            const slotCenterTheta = theta + (toothAngle / 2);
            const windingRadius = (geometry.statorInnerRadius + geometry.slotDepth / 2);
            const { x, y } = pol2car(windingRadius, slotCenterTheta);

            g.append("circle")
                .attr("cx", x)
                .attr("cy", y)
                .attr("r", Math.min(geometry.slotDepth, (geometry.statorInnerRadius * toothAngle)) / 3)
                .attr("fill", "#b87333") // Copper
                .attr("opacity", 0.8);
        }

        // --- Draw Rotor ---

        // Rotor Iron
        g.append("circle")
            .attr("r", geometry.rotorOuterRadius)
            .attr("fill", "#cbd5e1") // Slate 300
            .attr("stroke", "#475569")
            .attr("stroke-width", 1 / scale);

        // Magnets
        const poleAngle = (2 * Math.PI) / geometry.poles;
        for (let i = 0; i < geometry.poles; i++) {
            const theta = i * poleAngle;
            // Arc segment for magnet
            const arcGen = d3.arc()
                .innerRadius(geometry.rotorOuterRadius - geometry.magnetThickness)
                .outerRadius(geometry.rotorOuterRadius)
                .startAngle(theta)
                .endAngle(theta + poleAngle)
                .padAngle(0.02);

            g.append("path")
                .attr("d", arcGen as any)
                .attr("fill", i % 2 === 0 ? "#ef4444" : "#3b82f6") // N/S poles
                .attr("opacity", 0.9);
        }

        // Shaft
        g.append("circle")
            .attr("r", geometry.rotorInnerRadius)
            .attr("fill", "#94a3b8") // Shaft color
            .attr("stroke", "#0f172a");

        // Center Crosshair
        g.append("line").attr("x1", -5).attr("y1", 0).attr("x2", 5).attr("y2", 0).attr("stroke", "cyan").attr("stroke-width", 2 / scale);
        g.append("line").attr("x1", 0).attr("y1", -5).attr("x2", 0).attr("y2", 5).attr("stroke", "cyan").attr("stroke-width", 2 / scale);

    }, [geometry]);

    return (
        <div className="w-full aspect-square bg-slate-900 rounded-lg border border-slate-700 overflow-hidden relative shadow-inner shadow-black/50">
            <svg ref={svgRef} width="100%" height="100%" viewBox="0 0 400 400" preserveAspectRatio="xMidYMid meet" />
            <div className="absolute bottom-2 right-2 text-xs text-slate-500 font-mono">
                D_out: {(geometry.statorOuterRadius * 2).toFixed(1)}mm
            </div>
        </div>
    );
};

export default MotorVisualizer;
