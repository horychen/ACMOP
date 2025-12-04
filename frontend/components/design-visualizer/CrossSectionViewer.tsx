"use client";

import React, { useEffect, useRef, useMemo } from 'react';
import * as d3 from 'd3';
import { GeometricComponentsObjects, GeometricComponent, GP } from '@/lib/DesignData';
import { useTheme } from '@/context/ThemeContext';

interface CrossSectionViewerProps {
  geometry: GeometricComponentsObjects;
  showParameters?: boolean;
  gpData?: GP;
}

export default function CrossSectionViewer({ 
  geometry, 
  showParameters = false, 
  gpData 
}: CrossSectionViewerProps) {
  const svgRef = useRef<SVGSVGElement>(null);
  const containerRef = useRef<HTMLDivElement>(null);
  const { theme } = useTheme();

  // Color mapping for components
  const componentColors = useMemo(() => {
    const baseColors: Record<string, string> = {
      rotorCore: theme === 'dark' ? '#cbd5e1' : '#e2e8f0',
      shaft: theme === 'dark' ? '#94a3b8' : '#64748b',
      rotorMagnet: theme === 'dark' ? '#ef4444' : '#dc2626',
      sleeve: theme === 'dark' ? '#64748b' : '#475569',
      statorCore: theme === 'dark' ? '#334155' : '#64748b',
      coils: theme === 'dark' ? '#b87333' : '#d97706',
    };
    return baseColors;
  }, [theme]);

  // Calculate bounds and scale
  const { bounds, scale, center } = useMemo(() => {
    if (!geometry) return { bounds: null, scale: 1, center: { x: 0, y: 0 } };

    let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;

    // Extract all points from all components
    Object.values(geometry).forEach((component: GeometricComponent | null) => {
      if (!component) return;
      
      // Try to get points from list_region first
      if (component.list_region && Array.isArray(component.list_region)) {
        component.list_region.forEach((region: any[]) => {
          if (Array.isArray(region)) {
            region.forEach((item: any) => {
              if (item.move_to && Array.isArray(item.move_to)) {
                const [x, y] = item.move_to;
                minX = Math.min(minX, x);
                minY = Math.min(minY, y);
                maxX = Math.max(maxX, x);
                maxY = Math.max(maxY, y);
              }
              if (item.line_to && Array.isArray(item.line_to)) {
                const [x, y] = item.line_to;
                minX = Math.min(minX, x);
                minY = Math.min(minY, y);
                maxX = Math.max(maxX, x);
                maxY = Math.max(maxY, y);
              }
            });
          }
        });
      }

      // Also check for point properties (P1, P2, etc.)
      Object.keys(component).forEach(key => {
        if (/^P\d+/.test(key) || /^P[A-Za-z]+/.test(key)) {
          const val = component[key];
          if (Array.isArray(val) && val.length === 2) {
            const [x, y] = val;
            if (typeof x === 'number' && typeof y === 'number') {
              minX = Math.min(minX, x);
              minY = Math.min(minY, y);
              maxX = Math.max(maxX, x);
              maxY = Math.max(maxY, y);
            }
          }
        }
      });
    });

    if (minX === Infinity) {
      return { bounds: null, scale: 1, center: { x: 0, y: 0 } };
    }

    const width = maxX - minX;
    const height = maxY - minY;
    const centerX = (minX + maxX) / 2;
    const centerY = (minY + maxY) / 2;
    const maxDim = Math.max(width, height);

    return {
      bounds: { minX, minY, maxX, maxY, width, height },
      scale: maxDim > 0 ? 1 : 1,
      center: { x: centerX, y: centerY }
    };
  }, [geometry]);

  useEffect(() => {
    if (!svgRef.current || !geometry || !bounds) return;

    const svg = d3.select(svgRef.current);
    svg.selectAll("*").remove();

    const container = containerRef.current;
    if (!container) return;

    const width = container.clientWidth;
    const height = container.clientHeight;
    const padding = 20;

    // Calculate scale to fit
    const scaleX = (width - 2 * padding) / bounds.width;
    const scaleY = (height - 2 * padding) / bounds.height;
    const viewScale = Math.min(scaleX, scaleY, 1) * 0.9; // 90% to leave some padding

    // Center the view
    const viewBoxX = center.x - (width / viewScale) / 2;
    const viewBoxY = center.y - (height / viewScale) / 2;

    svg.attr("width", width)
       .attr("height", height)
       .attr("viewBox", `${viewBoxX} ${viewBoxY} ${width / viewScale} ${height / viewScale}`)
       .attr("preserveAspectRatio", "xMidYMid meet");

    const g = svg.append("g");

    // Draw components in order (back to front)
    const drawOrder: (keyof GeometricComponentsObjects)[] = [
      'statorCore',
      'rotorCore',
      'rotorMagnet',
      'sleeve',
      'shaft',
      'coils'
    ];

    drawOrder.forEach(componentKey => {
      const component = geometry[componentKey];
      if (!component) return;

      const color = component.color || componentColors[componentKey] || '#888888';
      const strokeColor = theme === 'dark' ? '#94a3b8' : '#475569';
      const strokeWidth = 0.5;

      // Draw from list_region if available
      if (component.list_region && Array.isArray(component.list_region)) {
        component.list_region.forEach((region: any[]) => {
          if (!Array.isArray(region)) return;

          const pathData: string[] = [];
          let currentX = 0, currentY = 0;

          region.forEach((item: any, index: number) => {
            if (item.move_to && Array.isArray(item.move_to)) {
              const [x, y] = item.move_to;
              currentX = x;
              currentY = y;
              if (index === 0) {
                pathData.push(`M ${x} ${y}`);
              } else {
                pathData.push(`M ${x} ${y}`);
              }
            } else if (item.line_to && Array.isArray(item.line_to)) {
              const [x, y] = item.line_to;
              pathData.push(`L ${x} ${y}`);
              currentX = x;
              currentY = y;
            } else if (item.arc && Array.isArray(item.arc)) {
              // Handle arc: [centerX, centerY, radius, startAngle, endAngle]
              const [cx, cy, r, startAngle, endAngle] = item.arc;
              const startX = cx + r * Math.cos(startAngle);
              const startY = cy + r * Math.sin(startAngle);
              const endX = cx + r * Math.cos(endAngle);
              const endY = cy + r * Math.sin(endAngle);
              const largeArc = Math.abs(endAngle - startAngle) > Math.PI ? 1 : 0;
              const sweep = endAngle > startAngle ? 1 : 0;
              pathData.push(`A ${r} ${r} 0 ${largeArc} ${sweep} ${endX} ${endY}`);
            }
          });

          if (pathData.length > 0) {
            const pathString = pathData.join(' ');
            g.append("path")
              .attr("d", pathString)
              .attr("fill", color)
              .attr("stroke", strokeColor)
              .attr("stroke-width", strokeWidth)
              .attr("opacity", 0.9);
          }
        });
      } else {
        // Fallback: try to draw from point properties
        const points: [number, number][] = [];
        Object.keys(component).forEach(key => {
          if (/^P\d+/.test(key) || /^P[A-Za-z]+/.test(key)) {
            const val = component[key];
            if (Array.isArray(val) && val.length === 2) {
              const [x, y] = val;
              if (typeof x === 'number' && typeof y === 'number') {
                points.push([x, y]);
              }
            }
          }
        });

        if (points.length > 0) {
          // Sort points to form a closed path
          const sortedPoints = [...points].sort((a, b) => {
            const angleA = Math.atan2(a[1] - center.y, a[0] - center.x);
            const angleB = Math.atan2(b[1] - center.y, b[0] - center.x);
            return angleA - angleB;
          });

          const pathString = sortedPoints.map((p, i) => 
            i === 0 ? `M ${p[0]} ${p[1]}` : `L ${p[0]} ${p[1]}`
          ).join(' ') + ' Z';

          g.append("path")
            .attr("d", pathString)
            .attr("fill", color)
            .attr("stroke", strokeColor)
            .attr("stroke-width", strokeWidth)
            .attr("opacity", 0.9);
        }
      }
    });

    // Draw center crosshair
    g.append("line")
      .attr("x1", center.x - 5)
      .attr("y1", center.y)
      .attr("x2", center.x + 5)
      .attr("y2", center.y)
      .attr("stroke", theme === 'dark' ? '#06b6d4' : '#0891b2')
      .attr("stroke-width", 0.5);

    g.append("line")
      .attr("x1", center.x)
      .attr("y1", center.y - 5)
      .attr("x2", center.x)
      .attr("y2", center.y + 5)
      .attr("stroke", theme === 'dark' ? '#06b6d4' : '#0891b2')
      .attr("stroke-width", 0.5);

  }, [geometry, bounds, center, componentColors, theme]);

  if (!geometry || !bounds) {
    return (
      <div className="w-full h-full flex items-center justify-center text-muted-foreground text-sm">
        无法渲染几何图形
      </div>
    );
  }

  return (
    <div ref={containerRef} className="w-full h-full relative">
      <svg ref={svgRef} className="w-full h-full" />
      {showParameters && gpData && (
        <div className="absolute top-2 right-2 bg-background/80 backdrop-blur-sm border border-border rounded-lg p-2 text-xs max-w-xs">
          <div className="font-semibold mb-1">关键参数</div>
          <div className="space-y-1">
            {gpData.mm_r_so && (
              <div className="flex justify-between">
                <span className="text-muted-foreground">r_so:</span>
                <span>{gpData.mm_r_so.value?.toFixed(2)}</span>
              </div>
            )}
            {gpData.mm_r_si && (
              <div className="flex justify-between">
                <span className="text-muted-foreground">r_si:</span>
                <span>{gpData.mm_r_si.value?.toFixed(2)}</span>
              </div>
            )}
            {gpData.mm_r_ro && (
              <div className="flex justify-between">
                <span className="text-muted-foreground">r_ro:</span>
                <span>{gpData.mm_r_ro.value?.toFixed(2)}</span>
              </div>
            )}
            {gpData.mm_d_pm && (
              <div className="flex justify-between">
                <span className="text-muted-foreground">d_pm:</span>
                <span>{gpData.mm_d_pm.value?.toFixed(2)}</span>
              </div>
            )}
          </div>
        </div>
      )}
    </div>
  );
}

