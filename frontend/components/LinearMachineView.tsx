"use client";

import React, { useEffect, useRef, useMemo } from 'react';
import * as d3 from 'd3';
import { MachineGeometry } from '../types';
import { useTheme } from '@/context/ThemeContext';

interface LinearMachineViewProps {
    Qs: number;  // Stator slot number
    p: number;   // Pole pair number (for torque)
    ps: number;  // Suspension pole pair number
    coilPitchY?: number | null;
    layer_X_phases?: string[] | null;  // Phase assignment for Layer X (upper conductors)
    layer_Y_phases?: string[] | null;  // Phase assignment for Layer Y (lower conductors)
    layer_X_signs?: string[] | null;  // Sign assignment for Layer X
    layer_Y_signs?: string[] | null;  // Sign assignment for Layer Y
    grouping_AC?: (number | boolean)[] | null;  // Grouping AC information
}

const LinearMachineView: React.FC<LinearMachineViewProps> = ({ 
    Qs, 
    p, 
    ps, 
    coilPitchY,
    layer_X_phases,
    layer_Y_phases,
    layer_X_signs,
    layer_Y_signs,
    grouping_AC
}) => {
    const svgRef = useRef<SVGSVGElement>(null);
    const { theme } = useTheme();

    // Calculate coil distribution (distributed winding)
    const coilDistribution = useMemo(() => {
        const slots = Qs;
        const poles = p * 2; // p is pole pairs, so poles = p * 2
        const slotsPerPole = slots / poles;
        const phases = 3; // 3-phase winding
        const slotsPerPolePerPhase = slotsPerPole / phases;
        
        // Generate coil sides for each slot (2 conductors per slot)
        // Each coil has two sides: one in top layer, one in bottom layer
        const coils: Array<{ slot: number; topConductor: number | null; bottomConductor: number | null }> = [];
        
        // Calculate total number of coils
        const totalCoils = slots; // Each slot contributes to coils
        
        // Assign conductors to slots
        // For a typical distributed winding, conductors are numbered sequentially
        // Each slot has 2 conductors (top and bottom layers)
        for (let slot = 0; slot < slots; slot++) {
            // Top conductor: slot number + 1 (1-indexed)
            // Bottom conductor: connects to another slot's top conductor
            // For simplicity, we'll number them sequentially
            const topConductor = slot * 2 + 1;
            const bottomConductor = slot * 2 + 2;
            
            coils.push({
                slot,
                topConductor,
                bottomConductor: bottomConductor <= totalCoils * 2 ? bottomConductor : null
            });
        }
        
        return coils;
    }, [Qs, p]);

    // Phase colors - different colors for U, V, W phases
    const phaseColors = useMemo(() => theme === 'dark' ? {
        U: { fill: "#ef4444", stroke: "#dc2626" },      // Red for U phase
        V: { fill: "#3b82f6", stroke: "#2563eb" },     // Blue for V phase
        W: { fill: "#10b981", stroke: "#059669" },     // Green for W phase
        default: { fill: "#b87333", stroke: "#d97706" } // Copper for unknown
    } : {
        U: { fill: "#dc2626", stroke: "#b91c1c" },     // Red for U phase
        V: { fill: "#2563eb", stroke: "#1d4ed8" },     // Blue for V phase
        W: { fill: "#059669", stroke: "#047857" },     // Green for W phase
        default: { fill: "#d97706", stroke: "#f59e0b" } // Amber for unknown
    }, [theme]);

    // Color schemes based on theme
    const colors = useMemo(() => theme === 'dark' ? {
        statorYoke: "#334155",      // Slate 700
        statorTeeth: "#475569",     // Slate 600
        slotBackground: "#0f172a",  // Slate 950
        conductor: "#b87333",       // Copper
        conductorStroke: "#d97706",  // Amber 600
        text: "#e2e8f0",             // Slate 200
        textSecondary: "#94a3b8",    // Slate 400
        border: "#475569",            // Slate 600
        background: "#0f172a"        // Slate 950
    } : {
        statorYoke: "#64748b",       // Slate 500
        statorTeeth: "#94a3b8",      // Slate 400
        slotBackground: "#f8fafc",   // Slate 50
        conductor: "#d97706",        // Amber 600 (brighter copper)
        conductorStroke: "#f59e0b",  // Amber 500
        text: "#1e293b",             // Slate 800
        textSecondary: "#475569",    // Slate 600
        border: "#cbd5e1",           // Slate 300
        background: "#f8fafc"        // Slate 50
    }, [theme]);

    useEffect(() => {
        if (!svgRef.current) return;

        const drawSVG = () => {
            if (!svgRef.current) return;

            const svg = d3.select(svgRef.current);
            svg.selectAll("*").remove();

            // Calculate required dimensions first
            const statorYokeHeight = 40;
            const slotDepth = 80;
            const toothHeight = 20;
            const airGapHeight = 10;
            const rotorHeight = 60;
            const secondRowHeight = 40; // Height for second row of magnets (ps)
            const netForceRowHeight = 35; // Height for net force visualization row
            const legendHeight = 30; // Space for legend
            const titleHeight = 30; // Space for title
            
            // Calculate margins
            const margin = { top: 40, right: 20, bottom: 50, left: 60 };
            
            // Calculate actual container dimensions
            const container = svgRef.current.parentElement;
            const containerWidth = container?.clientWidth || 800;
            const containerHeight = container?.clientHeight || 500;
            
            // Use container width, but ensure minimum width
            const minInnerWidth = 600;
            const width = Math.max(containerWidth, minInnerWidth + margin.left + margin.right);
            
            // Update SVG dimensions - set width first
            svg.attr("width", width);
            
            const innerWidth = width - margin.left - margin.right;

        const g = svg.append("g")
            .attr("transform", `translate(${margin.left}, ${margin.top})`);

        // Calculate dimensions
        const slotWidth = innerWidth / Qs;

        // Calculate the actual content height (without legend)
        const actualContentHeight = statorYokeHeight + slotDepth + toothHeight + airGapHeight + rotorHeight + secondRowHeight + netForceRowHeight;
        
        // Calculate total height needed including legend
        const legendSpace = 40; // Space needed for legend
        const totalNeededHeight = actualContentHeight + legendSpace;
        
        // Determine scaling strategy
        // If container is large enough, use it; otherwise, use full size and let container scroll
        const minContainerHeight = 400; // Minimum container height before we start scaling
        let scaleY: number;
        let requiredHeight: number;
        
        if (containerHeight >= minContainerHeight && containerHeight >= totalNeededHeight + margin.top + margin.bottom) {
            // Container is large enough - use container height, no scaling needed
            scaleY = 1.0;
            requiredHeight = containerHeight;
        } else if (containerHeight >= minContainerHeight) {
            // Container is reasonable size but content is too large - scale down to fit
            const availableInnerHeight = containerHeight - margin.top - margin.bottom;
            scaleY = Math.max((availableInnerHeight - legendSpace) / actualContentHeight, 0.5); // Minimum scale 0.5
            const scaledContentHeight = actualContentHeight * scaleY;
            requiredHeight = scaledContentHeight + legendSpace + margin.top + margin.bottom;
        } else {
            // Container is too small - use full size (1.0 scale), let container scroll
            scaleY = 1.0;
            requiredHeight = totalNeededHeight + margin.top + margin.bottom;
        }
        
        // Update SVG height to fit all content
        svg.attr("height", requiredHeight);
        
        // Calculate inner height for drawing
        const innerHeight = requiredHeight - margin.top - margin.bottom;
        const scaledYoke = statorYokeHeight * scaleY;
        const scaledSlot = slotDepth * scaleY;
        const scaledTooth = toothHeight * scaleY;
        const scaledAirGap = airGapHeight * scaleY;
        const scaledRotor = rotorHeight * scaleY;
        const scaledSecondRow = secondRowHeight * scaleY;
        const scaledNetForceRow = netForceRowHeight * scaleY;

        // Draw Stator Yoke (top)
        g.append("rect")
            .attr("x", 0)
            .attr("y", 0)
            .attr("width", innerWidth)
            .attr("height", scaledYoke)
            .attr("fill", colors.statorYoke)
            .attr("stroke", colors.border)
            .attr("stroke-width", 1);

        // Draw Slots and Teeth
        for (let i = 0; i < Qs; i++) {
            const x = i * slotWidth;
            const slotX = x;
            const toothX = x + slotWidth * 0.6; // Tooth takes 60% of slot pitch
            const toothW = slotWidth * 0.4;     // Tooth width

            // Slot (rectangular opening)
            g.append("rect")
                .attr("x", slotX)
                .attr("y", scaledYoke)
                .attr("width", slotWidth)
                .attr("height", scaledSlot)
                .attr("fill", colors.slotBackground)
                .attr("stroke", colors.border)
                .attr("stroke-width", 0.5);

            // Tooth
            g.append("rect")
                .attr("x", toothX)
                .attr("y", scaledYoke)
                .attr("width", toothW)
                .attr("height", scaledSlot)
                .attr("fill", colors.statorTeeth)
                .attr("stroke", colors.border)
                .attr("stroke-width", 0.5);

            // Conductors in slot (2 per slot - always show both)
            const coil = coilDistribution[i];
            const conductorRadius = Math.min(slotWidth * 0.15, scaledSlot * 0.15);
            const topY = scaledYoke + scaledSlot * 0.3;
            const bottomY = scaledYoke + scaledSlot * 0.7;

            // Get phase information for this slot
            const layerXPhase = layer_X_phases && i < layer_X_phases.length ? layer_X_phases[i] : null;
            const layerYPhase = layer_Y_phases && i < layer_Y_phases.length ? layer_Y_phases[i] : null;
            const layerXSign = layer_X_signs && i < layer_X_signs.length ? layer_X_signs[i] : null;
            const layerYSign = layer_Y_signs && i < layer_Y_signs.length ? layer_Y_signs[i] : null;
            const isAC = grouping_AC && i < grouping_AC.length ? (grouping_AC[i] === 1 || grouping_AC[i] === true) : false;

            // Determine colors for conductors based on phase
            const topPhaseColor = layerXPhase && phaseColors[layerXPhase as 'U' | 'V' | 'W'] 
                ? phaseColors[layerXPhase as 'U' | 'V' | 'W'] 
                : phaseColors.default;
            const bottomPhaseColor = layerYPhase && phaseColors[layerYPhase as 'U' | 'V' | 'W'] 
                ? phaseColors[layerYPhase as 'U' | 'V' | 'W'] 
                : phaseColors.default;

            // Top conductor (Layer X - upper conductor)
            g.append("circle")
                .attr("cx", slotX + slotWidth / 2)
                .attr("cy", topY)
                .attr("r", conductorRadius)
                .attr("fill", topPhaseColor.fill)
                .attr("stroke", topPhaseColor.stroke)
                .attr("stroke-width", 1.5);

            // AC label for top conductor if grouping_AC is true
            if (isAC) {
                g.append("text")
                    .attr("x", slotX + slotWidth / 2)
                    .attr("y", topY + 3)
                    .attr("text-anchor", "middle")
                    .attr("font-size", "7px")
                    .attr("font-weight", "bold")
                    .attr("fill", theme === 'dark' ? "#ffffff" : "#000000")
                    .text("AC");
            }

            // Slot number - marked above the upper conductor
            g.append("text")
                .attr("x", slotX + slotWidth / 2)
                .attr("y", topY - conductorRadius - 15)
                .attr("text-anchor", "middle")
                .attr("font-size", "9px")
                .attr("font-weight", "bold")
                .attr("fill", colors.textSecondary)
                .text(`S${i + 1}`);

            // Bottom conductor (Layer Y - lower conductor)
            g.append("circle")
                .attr("cx", slotX + slotWidth / 2)
                .attr("cy", bottomY)
                .attr("r", conductorRadius)
                .attr("fill", bottomPhaseColor.fill)
                .attr("stroke", bottomPhaseColor.stroke)
                .attr("stroke-width", 1.5);

            // AC label for bottom conductor if grouping_AC is true
            // Note: grouping_AC typically applies to Layer X, but we show it on both if needed
            if (isAC) {
                g.append("text")
                    .attr("x", slotX + slotWidth / 2)
                    .attr("y", bottomY + 3)
                    .attr("text-anchor", "middle")
                    .attr("font-size", "7px")
                    .attr("font-weight", "bold")
                    .attr("fill", theme === 'dark' ? "#ffffff" : "#000000")
                    .text("AC");
            }
        }

        // Air gap line
        const airGapY = scaledYoke + scaledSlot;
        g.append("line")
            .attr("x1", 0)
            .attr("y1", airGapY)
            .attr("x2", innerWidth)
            .attr("y2", airGapY)
            .attr("stroke", colors.border)
            .attr("stroke-width", 2)
            .attr("stroke-dasharray", "4,4");

        // Rotor (simplified as a rectangle)
        g.append("rect")
            .attr("x", 0)
            .attr("y", airGapY)
            .attr("width", innerWidth)
            .attr("height", scaledRotor)
            .attr("fill", theme === 'dark' ? "#cbd5e1" : "#e2e8f0")
            .attr("stroke", colors.border)
            .attr("stroke-width", 1);

        // First row of rotor poles (p - for torque)
        const poles = p * 2; // p is pole pairs, so total poles = p * 2
        const poleWidth = innerWidth / poles;
        for (let i = 0; i < poles; i++) {
            const poleX = i * poleWidth;
            g.append("rect")
                .attr("x", poleX)
                .attr("y", airGapY)
                .attr("width", poleWidth)
                .attr("height", scaledRotor)
                .attr("fill", i % 2 === 0 
                    ? (theme === 'dark' ? "#ef4444" : "#dc2626")
                    : (theme === 'dark' ? "#3b82f6" : "#2563eb"))
                .attr("opacity", 0.3)
                .attr("stroke", colors.border)
                .attr("stroke-width", 0.5);

            // Pole label
            g.append("text")
                .attr("x", poleX + poleWidth / 2)
                .attr("y", airGapY + scaledRotor / 2)
                .attr("text-anchor", "middle")
                .attr("font-size", "10px")
                .attr("font-weight", "bold")
                .attr("fill", colors.text)
                .text(i % 2 === 0 ? "N" : "S");
        }

        // Second row of magnets (ps - for suspension)
        const secondRowY = airGapY + scaledRotor;
        const suspensionPoles = ps * 2; // ps is pole pairs, so total poles = ps * 2
        const suspensionPoleWidth = innerWidth / suspensionPoles;
        for (let i = 0; i < suspensionPoles; i++) {
            const poleX = i * suspensionPoleWidth;
            g.append("rect")
                .attr("x", poleX)
                .attr("y", secondRowY)
                .attr("width", suspensionPoleWidth)
                .attr("height", scaledSecondRow)
                .attr("fill", i % 2 === 0 
                    ? (theme === 'dark' ? "#f59e0b" : "#d97706")
                    : (theme === 'dark' ? "#10b981" : "#059669"))
                .attr("opacity", 0.4)
                .attr("stroke", colors.border)
                .attr("stroke-width", 0.5);

            // Suspension pole label
            g.append("text")
                .attr("x", poleX + suspensionPoleWidth / 2)
                .attr("y", secondRowY + scaledSecondRow / 2)
                .attr("text-anchor", "middle")
                .attr("font-size", "9px")
                .attr("font-weight", "bold")
                .attr("fill", colors.text)
                .text(i % 2 === 0 ? "N" : "S");
        }

        // Net Force Row - visualize the interaction between torque and suspension poles
        const netForceRowY = secondRowY + scaledSecondRow;
        
        // Calculate net force segments across the width
        // We'll create segments based on the least common multiple of pole counts for accurate visualization
        const forceSegments = Math.max(poles, suspensionPoles) * 4; // Enough segments for smooth visualization
        const segmentWidth = innerWidth / forceSegments;
        
        // Draw net force visualization as segmented bars
        for (let i = 0; i < forceSegments; i++) {
            const x = i * segmentWidth;
            const position = x / innerWidth; // Normalized position (0 to 1)
            
            // Determine which torque pole and suspension pole are at this position
            const torquePoleIndex = Math.floor(position * poles) % poles;
            const suspensionPoleIndex = Math.floor(position * suspensionPoles) % suspensionPoles;
            
            // Determine pole polarities
            // Torque: even index = N (0), odd index = S (1)
            const torquePoleType = torquePoleIndex % 2 === 0 ? 'N' : 'S';
            // Suspension: even index = N (0), odd index = S (1)
            const suspensionPoleType = suspensionPoleIndex % 2 === 0 ? 'N' : 'S';
            
            // Calculate net force
            // N + N = repulsion (negative force, minimum) - Red
            // S + S = repulsion (negative force, minimum) - Red
            // N + S = attraction (positive force, maximum) - Green
            // S + N = attraction (positive force, maximum) - Green
            const isRepulsion = (torquePoleType === 'N' && suspensionPoleType === 'N') || 
                               (torquePoleType === 'S' && suspensionPoleType === 'S');
            
            // Draw force segment
            // Repulsion: extends upward from center (negative force)
            // Attraction: extends downward from center (positive force)
            const centerY = netForceRowY + scaledNetForceRow / 2;
            const forceHeight = scaledNetForceRow * 0.4; // 40% of row height for force visualization
            
            if (isRepulsion) {
                // Repulsion - red, extends upward
                g.append("rect")
                    .attr("x", x)
                    .attr("y", centerY - forceHeight)
                    .attr("width", segmentWidth)
                    .attr("height", forceHeight)
                    .attr("fill", theme === 'dark' ? "#ef4444" : "#dc2626")
                    .attr("opacity", 0.7)
                    .attr("stroke", colors.border)
                    .attr("stroke-width", 0.5);
                
                // Up arrow indicator for repulsion
                if (i % Math.floor(forceSegments / 8) === 0) {
                    g.append("path")
                        .attr("d", `M ${x + segmentWidth / 2} ${centerY - forceHeight * 0.8} L ${x + segmentWidth / 2 - 3} ${centerY - forceHeight * 0.6} L ${x + segmentWidth / 2 + 3} ${centerY - forceHeight * 0.6} Z`)
                        .attr("fill", theme === 'dark' ? "#ef4444" : "#dc2626")
                        .attr("opacity", 0.9);
                }
            } else {
                // Attraction - green, extends downward
                g.append("rect")
                    .attr("x", x)
                    .attr("y", centerY)
                    .attr("width", segmentWidth)
                    .attr("height", forceHeight)
                    .attr("fill", theme === 'dark' ? "#10b981" : "#059669")
                    .attr("opacity", 0.7)
                    .attr("stroke", colors.border)
                    .attr("stroke-width", 0.5);
                
                // Down arrow indicator for attraction
                if (i % Math.floor(forceSegments / 8) === 0) {
                    g.append("path")
                        .attr("d", `M ${x + segmentWidth / 2} ${centerY + forceHeight * 0.8} L ${x + segmentWidth / 2 - 3} ${centerY + forceHeight * 0.6} L ${x + segmentWidth / 2 + 3} ${centerY + forceHeight * 0.6} Z`)
                        .attr("fill", theme === 'dark' ? "#10b981" : "#059669")
                        .attr("opacity", 0.9);
                }
            }
        }
        
        // Draw center line (neutral force reference)
        g.append("line")
            .attr("x1", 0)
            .attr("y1", netForceRowY + scaledNetForceRow / 2)
            .attr("x2", innerWidth)
            .attr("y2", netForceRowY + scaledNetForceRow / 2)
            .attr("stroke", colors.border)
            .attr("stroke-width", 1.5)
            .attr("stroke-dasharray", "3,3")
            .attr("opacity", 0.6);
        
        // Add label for net force row
        g.append("text")
            .attr("x", 5)
            .attr("y", netForceRowY + scaledNetForceRow / 2 + 4)
            .attr("font-size", "10px")
            .attr("font-weight", "bold")
            .attr("fill", colors.textSecondary)
            .text("Net Force");
        
        // Add force magnitude labels
        g.append("text")
            .attr("x", innerWidth - 80)
            .attr("y", netForceRowY + scaledNetForceRow * 0.25)
            .attr("font-size", "8px")
            .attr("fill", theme === 'dark' ? "#ef4444" : "#dc2626")
            .text("Min (Repulsion)");
        
        g.append("text")
            .attr("x", innerWidth - 80)
            .attr("y", netForceRowY + scaledNetForceRow * 0.75)
            .attr("font-size", "8px")
            .attr("fill", theme === 'dark' ? "#10b981" : "#059669")
            .text("Max (Attraction)");

        // Title
        const titleText = coilPitchY !== null && coilPitchY !== undefined
            ? `Linear Machine View - Qs=${Qs}, p=${p}, ps=${ps}, coil_pitch_y=${coilPitchY}`
            : `Linear Machine View - Qs=${Qs}, p=${p}, ps=${ps}`;
        g.append("text")
            .attr("x", innerWidth / 2)
            .attr("y", -10)
            .attr("text-anchor", "middle")
            .attr("font-size", "14px")
            .attr("font-weight", "bold")
            .attr("fill", colors.text)
            .text(titleText);

        // Legend - position at the bottom of the scaled content area
        const scaledTotalContentHeight = actualContentHeight * scaleY;
        const legendY = scaledTotalContentHeight + 10;
        const legendItems = [
            { label: "Stator Yoke", color: colors.statorYoke },
            { label: "Stator Teeth", color: colors.statorTeeth },
            { label: `Torque Poles (p=${p})`, color: theme === 'dark' ? "#ef4444" : "#dc2626" },
            { label: `Suspension Poles (ps=${ps})`, color: theme === 'dark' ? "#f59e0b" : "#d97706" },
            { label: "Net Force (Repulsion)", color: theme === 'dark' ? "#ef4444" : "#dc2626" },
            { label: "Net Force (Attraction)", color: theme === 'dark' ? "#10b981" : "#059669" }
        ];

        // Add phase colors to legend if winding data is available
        if (layer_X_phases || layer_Y_phases) {
            legendItems.push(
                { label: "Phase U", color: phaseColors.U.fill },
                { label: "Phase V", color: phaseColors.V.fill },
                { label: "Phase W", color: phaseColors.W.fill }
            );
        }

        legendItems.forEach((item, i) => {
            const x = 20 + i * 150;
            g.append("rect")
                .attr("x", x)
                .attr("y", legendY)
                .attr("width", 12)
                .attr("height", 12)
                .attr("fill", item.color)
                .attr("stroke", colors.border);

            g.append("text")
                .attr("x", x + 18)
                .attr("y", legendY + 9)
                .attr("font-size", "10px")
                .attr("fill", colors.textSecondary)
                .text(item.label);
        });
        };

        // Initial draw
        drawSVG();

        // Set up ResizeObserver to redraw when container size changes
        const container = svgRef.current.parentElement;
        if (container && typeof ResizeObserver !== 'undefined') {
            const resizeObserver = new ResizeObserver(() => {
                drawSVG();
            });
            resizeObserver.observe(container);

            return () => {
                resizeObserver.disconnect();
            };
        }

        // Fallback: listen to window resize if ResizeObserver is not available
        const handleResize = () => {
            drawSVG();
        };
        window.addEventListener('resize', handleResize);

        return () => {
            window.removeEventListener('resize', handleResize);
        };
    }, [Qs, p, ps, colors, coilDistribution, theme, coilPitchY, layer_X_phases, layer_Y_phases, layer_X_signs, layer_Y_signs, grouping_AC, phaseColors]);

    return (
        <div className={`w-full rounded-lg border overflow-auto ${
            theme === 'dark' 
                ? 'bg-slate-900 border-slate-700' 
                : 'bg-white border-slate-200'
        }`} style={{ minHeight: '500px' }}>
            <svg 
                ref={svgRef} 
                width="100%" 
                className="w-full"
                style={{ display: 'block' }}
                preserveAspectRatio="none"
            />
        </div>
    );
};

export default LinearMachineView;

