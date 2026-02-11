"use client";

import React, { useEffect, useState, useMemo } from "react";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/alert";
import { Loader2, AlertTriangle, RefreshCw } from "lucide-react";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";

interface MachineSpecs {
  geometry: {
    d_stator_outer: number;
    d_rotor_outer: number;
    d_shaft: number;
    air_gap: number;
    magnet_thickness: number;
    tooth_width: number;
    tooth_depth: number;
    num_slots?: number; // Sometimes in winding
  };
  winding: {
    num_slots: number;
    conductors_per_slot: number;
    wire_diameter_with_insulation: number;
    num_poles: number;
  };
}

export default function MachineVisualizer() {
  const [initialSpecs, setInitialSpecs] = useState<MachineSpecs | null>(null);
  const [specs, setSpecs] = useState<MachineSpecs | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [conductorError, setConductorError] = useState<string | null>(null);

  const fetchSpecs = () => {
    setLoading(true);
    fetch("http://localhost:8000/api/machine-specs")
      .then((res) => {
        if (!res.ok) throw new Error("Failed to fetch machine specs");
        return res.json();
      })
      .then((data) => {
        setInitialSpecs(data);
        setSpecs(data);
        setLoading(false);
      })
      .catch((err) => {
        setError(err.message);
        setLoading(false);
      });
  };

  useEffect(() => {
    fetchSpecs();
  }, []);

  const handleParamChange = (section: 'geometry' | 'winding', key: string, value: string) => {
    if (!specs) return;
    const numValue = parseFloat(value);
    if (isNaN(numValue)) return;

    setSpecs({
      ...specs,
      [section]: {
        ...specs[section],
        [key]: numValue
      }
    });
  };

  const visualizationData = useMemo(() => {
    if (!specs) return null;

    const { geometry, winding } = specs;
    
    const R_so = geometry.d_stator_outer / 2;
    const R_ro = geometry.d_rotor_outer / 2;
    const gap = geometry.air_gap;
    const R_si = R_ro + gap;
    const toothDepth = geometry.tooth_depth;
    const R_sb = R_si + toothDepth; // Slot bottom (outer radius of slot)
    const toothWidth = geometry.tooth_width;
    const numSlots = winding.num_slots;
    const numPoles = winding.num_poles;
    const conductorsPerSlot = winding.conductors_per_slot;
    const wireDiam = winding.wire_diameter_with_insulation;
    const wireRad = wireDiam / 2;
    const magnetThickness = geometry.magnet_thickness;
    const R_shaft = geometry.d_shaft / 2;

    // --- Magnet Geometry ---
    // Magnets are usually surface mounted (SPM) or interior (IPM). 
    // Assuming SPM based on d_rotor_outer.
    // Inner Radius of Magnet = R_ro - magnetThickness
    const R_mag_in = R_ro - magnetThickness;
    
    // Check if back iron exists
    // Back iron is between Shaft and Magnet Inner
    const hasBackIron = R_mag_in > R_shaft + 0.01;
    
    const magnetPaths: {d: string, color: string}[] = [];
    if (numPoles > 0) {
        const poleAngle = 360 / numPoles;
        for (let i = 0; i < numPoles; i++) {
            const startAngle = i * poleAngle;
            const endAngle = (i + 1) * poleAngle;
            
            // Draw arc segment
            // Start Outer
            const getPt = (r: number, theta: number) => {
                const rad = (theta * Math.PI) / 180;
                return { x: r * Math.cos(rad), y: r * Math.sin(rad) };
            };

            const p1 = getPt(R_ro, startAngle);
            const p2 = getPt(R_ro, endAngle);
            // End Inner
            const p3 = getPt(R_mag_in, endAngle);
            const p4 = getPt(R_mag_in, startAngle);
            
            // Path: M p1 A ... p2 L p3 A ... p4 Z
            // Large arc flag: poleAngle > 180? (Unlikely for poles)
            
            let d = `M ${p1.x} ${p1.y} `;
            d += `A ${R_ro} ${R_ro} 0 0 1 ${p2.x} ${p2.y} `;
            d += `L ${p3.x} ${p3.y} `;
            d += `A ${R_mag_in} ${R_mag_in} 0 0 0 ${p4.x} ${p4.y} `;
            d += `Z`;
            
            // Colors: Low saturation Red (#fca5a5) and Blue (#93c5fd)
            // i % 2 === 0 ? Red : Blue
            magnetPaths.push({
                d,
                color: i % 2 === 0 ? "#fca5a5" : "#93c5fd"
            });
        }
    }

    // Generate Stator Teeth and Slots
    // We'll generate a path for the stator lamination
    // Center a tooth at angle 0.
    // Tooth edges are parallel to radial line.
    
    const statorPathCommands: string[] = [];
    
    // Helper to get point at radius r, angle theta (degrees)
    const getPt = (r: number, theta: number) => {
      const rad = (theta * Math.PI) / 180;
      return { x: r * Math.cos(rad), y: r * Math.sin(rad) };
    };

    // Helper to intersect line (parallel to ray at angle, offset by d) with circle radius r
    // Line: x*cos(angle) + y*sin(angle) = d? No.
    // Ray direction u = (cos a, sin a). Normal n = (-sin a, cos a).
    // Line: p . n = offset. 
    // If offset > 0, it's to the "left" of the ray.
    // Tooth centered at alpha. Width w.
    // Edge 1: Right of center (looking out). Offset = -w/2 along normal?
    // Let's stick to: Tooth center at alpha.
    // Edge 1 (CCW side): Parallel to alpha ray, shifted by -w/2 along normal (-sin a, cos a).
    // Edge 2 (CW side): Parallel to alpha ray, shifted by +w/2 along normal.
    // Wait, normal (-sin a, cos a) points 90 deg CCW.
    // So +w/2 is "left" (CCW), -w/2 is "right" (CW).
    
    // We need intersection of line P + t*U with circle |P+tU|^2 = R^2
    // Or simply: Line distance to origin is w/2.
    // Intersection with circle R:
    // h = w/2.
    // distance along line from projection of origin = sqrt(R^2 - h^2).
    // So point is: (d_along * u) + (h * n)
    
    const getToothPoints = (angle: number, r: number, w: number) => {
        const rad = (angle * Math.PI) / 180;
        const u = { x: Math.cos(rad), y: Math.sin(rad) };
        const n = { x: -Math.sin(rad), y: Math.cos(rad) }; // Normal pointing CCW
        
        const h = w / 2;
        // Check if R >= h
        if (r < h) return null; // Should not happen for valid geometry
        
        const d = Math.sqrt(r*r - h*h);
        
        // Two edges
        // Edge 1 (CCW side of tooth): offset +h along n
        const p1 = {
            x: d * u.x + h * n.x,
            y: d * u.y + h * n.y
        };
        // Edge 2 (CW side of tooth): offset -h along n
        const p2 = {
            x: d * u.x - h * n.x,
            y: d * u.y - h * n.y
        };
        return { cw: p2, ccw: p1 };
    };

    let d_path = "";
    
    // Start at first tooth (angle 0)
    // We need to trace the inner contour: Tooth Tip -> Slot Bottom -> Next Tooth Tip
    // Actually, usually stator is:
    // Outer Circle (CCW)
    // Inner contour (CW) to make a hole
    
    // Let's draw the lamination shape.
    // Outer circle
    d_path += `M ${R_so} 0 A ${R_so} ${R_so} 0 1 1 ${-R_so} 0 A ${R_so} ${R_so} 0 1 1 ${R_so} 0 Z `;
    
    // Inner contour (cutout)
    // We'll go CW.
    // For each slot i:
    // Tooth i center at i * 360/N.
    // Tooth i CCW edge at R_si.
    // Tooth i CCW edge at R_sb? No, slot bottom is arc.
    
    // Let's define the sequence of points for the inner bore + slots.
    // We want to subtract this from the outer circle.
    // Or just draw the stator as a filled path with evenodd rule.
    
    const innerPoints: {x: number, y: number, type: 'L'|'A', r?: number, sweep?: number}[] = [];
    
    for (let i = 0; i < numSlots; i++) {
        const angle = i * 360 / numSlots;
        // Tooth i
        const tipPoints = getToothPoints(angle, R_si, toothWidth);
        const rootPoints = getToothPoints(angle, R_sb, toothWidth);
        
        if (!tipPoints || !rootPoints) continue;
        
        // We are moving CW.
        // Previous Slot was i-1.
        // We arrive at Tooth i from the CCW side (p1).
        // Go to Tooth i CW side (p2) along the tip arc?
        // Wait, "Tooth Width" is the metal width.
        // "Slot" is the space.
        // So we trace:
        // 1. Along Tooth i tip (R_si) from CCW edge to CW edge? 
        //    No, if we go CW, we hit CCW edge first, then CW edge.
        //    Wait, angle increases CCW.
        //    So if we go CW (decreasing angle), we hit CCW edge (higher angle) then CW edge (lower angle).
        
        // Let's list vertices in CW order for the hole.
        // Start at Tooth 0 CCW edge (approx +something deg).
        
        // Vertices for Tooth i:
        // 1. Tip CCW (R_si)
        // 2. Root CCW (R_sb) -> Line
        // 3. Root CW (R_sb) -> Arc (Slot Bottom) - WAIT.
        // The slot is BETWEEN teeth.
        // So the "Slot Bottom" is between Tooth i and Tooth i+1?
        // No, the slot is the space. The metal is the tooth.
        // The "Yoke" connects the teeth at the bottom (outer radius).
        // Wait, usually teeth point INWARDS.
        // Stator OD = 13. ID = 8.3.
        // Yoke is at OD? Or ID?
        // Usually Stator Yoke is at the Outer Diameter. Teeth point inwards.
        // So Slot Bottom is at R_so - yoke_thickness?
        // Or R_si + tooth_depth?
        // R_si = 4.15. Tooth depth = 2.0. R_sb = 6.15.
        // R_so = 6.5.
        // So Yoke is from 6.15 to 6.5.
        // So the "Slot Bottom" is at R_sb = 6.15.
        // The "Slot Opening" is at R_si = 4.15.
        
        // Path trace (CW):
        // Start at Tooth i, CCW edge, Tip (R_si).
        // Line to Tooth i, CW edge, Tip (R_si)? No, that's the tooth tip face.
        // Then Line to Tooth i, CW edge, Root (R_sb)? No, that would cut through the tooth.
        // We are tracing the AIR boundary (hole).
        // The hole is the bore + slots.
        
        // Let's trace the METAL boundary (Stator).
        // Outer Circle (CCW).
        // Inner Boundary (CW):
        // Start at Tooth 0 CCW edge, Tip (R_si).
        // Arc to Tooth 0 CW edge, Tip (R_si). (This is the tooth tip face).
        // Line to Tooth 0 CW edge, Root (R_sb). (Side of slot).
        // Arc to Tooth -1 (or 11) CCW edge, Root (R_sb). (Slot bottom).
        // Line to Tooth -1 CCW edge, Tip (R_si). (Side of slot).
        // ... repeat.
        
        // Let's do this loop for i = 0 to 11 (CCW), but generate path commands.
        // We can just use a single path for the stator.
        // M (Start Point)
        // For each tooth i:
        //   L (Tip CW) -> This is the side of the previous slot?
        // Let's go CCW.
        // Start Tooth 0 CW edge, Tip.
        // L -> Tooth 0 CCW edge, Tip. (Tooth Tip Arc)
        // L -> Tooth 0 CCW edge, Root. (Tooth Side)
        // L -> Tooth 1 CW edge, Root. (Slot Bottom Arc)
        // L -> Tooth 1 CW edge, Tip. (Tooth Side)
        // ...
        
        // Points for Tooth i:
        const tip = getToothPoints(angle, R_si, toothWidth);
        const root = getToothPoints(angle, R_sb, toothWidth);
        
        if (!tip || !root) continue;
        
        // We need next tooth for slot bottom
        const nextAngle = (i + 1) * 360 / numSlots;
        const nextRoot = getToothPoints(nextAngle, R_sb, toothWidth);
        const nextTip = getToothPoints(nextAngle, R_si, toothWidth); // For next iteration start
        
        if (!nextRoot) continue;

        // Command sequence for one pitch (Tooth + Slot):
        // 1. Move to Tip CW (if first)
        // 2. Arc to Tip CCW (Tooth Tip)
        // 3. Line to Root CCW (Tooth Side)
        // 4. Arc to Next Root CW (Slot Bottom)
        // 5. Line to Next Tip CW (Tooth Side - Next)
        
        if (i === 0) {
            d_path += `M ${tip.cw.x} ${tip.cw.y} `;
        }
        
        // Arc Tip CW -> Tip CCW
        // Radius R_si. Large arc? No. Sweep? Yes (CCW).
        d_path += `A ${R_si} ${R_si} 0 0 1 ${tip.ccw.x} ${tip.ccw.y} `;
        
        // Line Tip CCW -> Root CCW
        d_path += `L ${root.ccw.x} ${root.ccw.y} `;
        
        // Arc Root CCW -> Next Root CW (Slot Bottom)
        // Radius R_sb.
        d_path += `A ${R_sb} ${R_sb} 0 0 1 ${nextRoot.cw.x} ${nextRoot.cw.y} `;
        
        // Line Next Root CW -> Next Tip CW
        // This will be the start of next iteration's arc
        d_path += `L ${nextTip!.cw.x} ${nextTip!.cw.y} `;
    }
    
    d_path += "Z"; // Close inner loop
    
    // Combine with outer circle for evenodd fill
    const fullStatorPath = `M ${R_so} 0 A ${R_so} ${R_so} 0 1 1 ${-R_so} 0 A ${R_so} ${R_so} 0 1 1 ${R_so} 0 Z ` + d_path;

    // --- Conductors ---
    // Calculate for one slot (between Tooth 0 and Tooth 1)
    // Slot boundaries:
    // 1. Tooth 0 CCW Edge (Line from Tip.ccw to Root.ccw)
    // 2. Tooth 1 CW Edge (Line from NextTip.cw to NextRoot.cw)
    // 3. Slot Bottom Arc (Root.ccw to NextRoot.cw)
    // 4. Slot Opening Arc (Tip.ccw to NextTip.cw) - effectively open or airgap
    
    // We'll pack in the polygon defined by these 4 curves.
    // For simplicity, we can approximate the arcs as lines for collision, or use exact math.
    
    const conductors: {x: number, y: number}[] = [];
    
    // Slot 0 is between Tooth 0 and Tooth 1.
    // Tooth 0 center 0 deg. Tooth 1 center 30 deg.
    // Slot center 15 deg.
    
    // Boundary Lines:
    // Left Wall (Tooth 0 side): Tooth 0 at 0 deg.
    // Normal to ray 0 is (0, 1).
    // CCW edge is shifted by +w/2 along normal. So y = w/2.
    const halfTooth = toothWidth / 2;
    
    // Right Wall (Tooth 1 side): Tooth 1 at angleStep.
    const angleStep = 360 / numSlots;
    const radStep = (angleStep * Math.PI) / 180;
    const n2 = { x: -Math.sin(radStep), y: Math.cos(radStep) }; // Normal to Tooth 1 ray
    // CW edge is shifted by -w/2 along normal.
    // Line eq: P . n2 = -w/2.
    
    // Algorithm:
    // 1. Generate candidate points.
    
    const validCircles: {x: number, y: number}[] = [];
    const r = wireRad;
    const padding = 0.02; // Insulation/Air margin
    const effR = r + padding;
    
    // Brute force packing with heuristic
    // Generate a hexagonal grid in the bounding box of the slot
    // Filter points inside the region.
    // Sort by distance to walls (preference for "hugging").
    
    // Bounding box for Slot 0 (approx):
    // x from R_si*cos(angleStep/2) to R_sb*cos(angleStep/2).
    // But we scan a bit wider to be safe.
    
    const candidates: {x: number, y: number}[] = [];
    const step = effR * 1.732; // Hex spacing (sqrt(3)*r for tightest packing is 2r * sin(60)? No. Vertical spacing for hex is r*sqrt(3). Horizontal is 2r.)
    // Let's use a fine grid or standard hex grid.
    // Hex grid: rows separated by r*sqrt(3). Points in row separated by 2r. Odd rows shifted by r.
    const rowHeight = effR * Math.sqrt(3);
    const colWidth = effR * 2.0;

    // Scan range
    // x is roughly radial. y is tangential.
    // We scan x from R_si to R_sb. y from 0 to R_sb.
    
    for (let x = R_si - 1; x < R_sb + 1; x += rowHeight) {
        for (let y = -1; y < R_sb + 1; y += colWidth) {
             // Hex offset for x? Usually we stack layers radially?
             // Let's treat 'x' as the "layer" direction (radial).
             // If we want layers to hug the tooth walls, maybe we should align grid to walls?
             // But walls are not parallel.
             // Standard hex grid is fine if dense enough.
             
             // Let's shift y based on x row index to make it hexagonal
             const xIndex = Math.round((x - (R_si-1)) / rowHeight);
             const yOff = (xIndex % 2 === 0) ? 0 : colWidth/2;
             
             candidates.push({x, y: y + yOff});
        }
    }
    
    // Filter function
    const isInsideSlot = (p: {x: number, y: number}) => {
        const pr = Math.sqrt(p.x*p.x + p.y*p.y);
        if (pr < R_si + effR || pr > R_sb - effR) return false;
        
        // Check Wall 1 (Tooth 0 CCW edge, y = w/2)
        // We need y > w/2 + effR
        if (p.y < halfTooth + effR) return false;
        
        // Check Wall 2 (Tooth 1 CW edge)
        // Line P.n2 = -w/2.
        // We need P.n2 < -w/2 - effR (since normal points CCW, and slot is CW of Tooth 1 center, so "negative" side relative to normal?
        // Let's re-verify sign.
        // Tooth 1 center ray C. C.n2 = 0.
        // Slot is CW of Tooth 1.
        // Angle of slot < Angle of Tooth 1.
        // Normal n2 is 90 deg CCW from Tooth 1.
        // So Slot is "approaching" the normal direction? No.
        // Tooth 1 is at angle alpha. Normal is alpha + 90.
        // Slot is at alpha - delta.
        // Dot product: cos( (alpha+90) - (alpha-delta) ) = cos(90 + delta) = -sin(delta).
        // Since delta > 0, -sin(delta) is negative.
        // So Slot is on the negative side of the normal plane passing through origin.
        // So P.n2 is negative.
        // The wall is at -w/2.
        // We want to be "more negative" (further CW) than the wall?
        // Or "less negative" (closer to Tooth 0)?
        // Tooth 0 is at 0. Tooth 1 at 30. Slot at 15.
        // Wall is "left" side of Tooth 1 (looking from origin).
        // We want to be to the "right" of that wall (towards Tooth 0)?
        // Wait, "Right Wall" in my previous comment was Tooth 1 side.
        // Looking from origin out:
        // Left is Tooth 1 (30 deg). Right is Tooth 0 (0 deg).
        // My coordinate system: y is "up" (towards 90 deg).
        // So Tooth 1 (30 deg) has higher y than Tooth 0 (0 deg).
        // So Tooth 1 is the "Top/Left" wall. Tooth 0 is "Bottom/Right" wall.
        // Let's stick to "Tooth 0 Wall" and "Tooth 1 Wall".
        
        // Tooth 0 Wall: y = w/2. We want y > w/2. (Since slot angle > 0). Correct.
        
        // Tooth 1 Wall: P.n2 = -w/2.
        // Slot is "below" Tooth 1 (smaller angle).
        // Normal n2 points "above" Tooth 1 (larger angle).
        // So we want P.n2 < -w/2.
        // And with spacing: P.n2 < -w/2 - effR.
        
        const val = p.x * n2.x + p.y * n2.y;
        if (val > -halfTooth - effR) return false;
        
        return true;
    };
    
    // Sort candidates by "hugging" score.
    // Score = min(dist to Wall 1, dist to Wall 2).
    // We want to fill from walls inwards?
    // Or maybe fill from Slot Bottom (R_sb) inwards?
    // "Tooth internal point can have one more layer" -> Suggests packing density is key.
    // Let's prioritize: 1. Back of slot (R_sb), 2. Walls.
    
    const candidatesWithScore = candidates.filter(isInsideSlot).map(p => {
        const dist1 = Math.abs(p.y - halfTooth); // Dist to Tooth 0
        const val = p.x * n2.x + p.y * n2.y;
        const dist2 = Math.abs(val - (-halfTooth)); // Dist to Tooth 1
        
        const distToBack = R_sb - Math.sqrt(p.x*p.x + p.y*p.y);
        
        // Heuristic: Minimize distance to walls and back.
        // score small = better.
        // Let's try to fill layers from back to front.
        // Primary sort: distToBack. Secondary: min(dist1, dist2).
        
        // Actually, just minimizing distance to ANY boundary (walls or back) usually gives good packing.
        const minWall = Math.min(dist1, dist2);
        return { p, score: distToBack + minWall * 0.1 }; 
    });
    
    candidatesWithScore.sort((a, b) => a.score - b.score);
    
    // Greedy placement
    const placed: {x: number, y: number}[] = [];
    for (const item of candidatesWithScore) {
        let overlap = false;
        for (const existing of placed) {
            const dx = item.p.x - existing.x;
            const dy = item.p.y - existing.y;
            if (dx*dx + dy*dy < (2*effR)*(2*effR) - 0.0001) {
                overlap = true;
                break;
            }
        }
        if (!overlap) {
            placed.push(item.p);
        }
    }
    
    // Check count
    let maxFit = placed.length;
    let finalConductors = placed;
    let errorMsg = null;
    
    if (placed.length < conductorsPerSlot) {
        errorMsg = `Cannot fit ${conductorsPerSlot} conductors. Max capacity: ${maxFit}.`;
    } else {
        finalConductors = placed.slice(0, conductorsPerSlot);
    }
    
    // Replicate for all slots
    const allConductors: {x: number, y: number}[] = [];
    for (let i = 0; i < numSlots; i++) {
        const angle = i * 360 / numSlots;
        const rad = angle * Math.PI / 180;
        const cos = Math.cos(rad);
        const sin = Math.sin(rad);
        
        finalConductors.forEach(c => {
            allConductors.push({
                x: c.x * cos - c.y * sin,
                y: c.x * sin + c.y * cos
            });
        });
    }

    return {
        R_so, R_si, R_ro, R_sb, R_shaft, R_mag_in, hasBackIron,
        fullStatorPath,
        allConductors,
        magnetPaths,
        errorMsg,
        maxFit,
        conductorsPerSlot
    };
  }, [specs]);

  if (loading) return <div className="flex justify-center p-8"><Loader2 className="animate-spin" /></div>;
  if (error) return <Alert variant="destructive"><AlertTitle>Error</AlertTitle><AlertDescription>{error}</AlertDescription></Alert>;
  if (!specs || !visualizationData) return <div>No data</div>;

  const { R_so, fullStatorPath, allConductors, magnetPaths, errorMsg, maxFit, conductorsPerSlot, R_ro, R_shaft, R_mag_in, hasBackIron } = visualizationData;
  const viewBoxSize = R_so * 2.4;
  
  return (
    <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
      {/* Controls */}
      <Card className="md:col-span-1">
        <CardHeader>
            <CardTitle className="flex items-center justify-between">
                Parameters
                <Button variant="outline" size="icon" onClick={fetchSpecs} title="Reload from File">
                    <RefreshCw className="h-4 w-4" />
                </Button>
            </CardTitle>
        </CardHeader>
        <CardContent className="space-y-4">
            <div className="space-y-2">
                <Label>Stator OD (mm)</Label>
                <Input 
                    type="number" 
                    value={specs.geometry.d_stator_outer} 
                    onChange={(e) => handleParamChange('geometry', 'd_stator_outer', e.target.value)}
                />
            </div>
            <div className="space-y-2">
                <Label>Rotor OD (mm)</Label>
                <Input 
                    type="number" 
                    value={specs.geometry.d_rotor_outer} 
                    onChange={(e) => handleParamChange('geometry', 'd_rotor_outer', e.target.value)}
                />
            </div>
            <div className="space-y-2">
                <Label>Tooth Width (mm)</Label>
                <Input 
                    type="number" 
                    value={specs.geometry.tooth_width} 
                    onChange={(e) => handleParamChange('geometry', 'tooth_width', e.target.value)}
                />
            </div>
            <div className="space-y-2">
                <Label>Tooth Depth (mm)</Label>
                <Input 
                    type="number" 
                    value={specs.geometry.tooth_depth} 
                    onChange={(e) => handleParamChange('geometry', 'tooth_depth', e.target.value)}
                />
            </div>
            <div className="space-y-2">
                <Label>Magnet Thickness (mm)</Label>
                <Input 
                    type="number" 
                    value={specs.geometry.magnet_thickness} 
                    onChange={(e) => handleParamChange('geometry', 'magnet_thickness', e.target.value)}
                />
            </div>
             <div className="space-y-2">
                <Label>Shaft Diameter (mm)</Label>
                <Input 
                    type="number" 
                    value={specs.geometry.d_shaft} 
                    onChange={(e) => handleParamChange('geometry', 'd_shaft', e.target.value)}
                />
            </div>
        </CardContent>
      </Card>

      {/* Visualization */}
      <Card className="md:col-span-2">
        <CardHeader>
            <CardTitle>Cross-Section View</CardTitle>
        </CardHeader>
        <CardContent className="flex flex-col items-center">
            <div className="relative border rounded p-4 bg-white">
                <svg width="600" height="600" viewBox={`${-viewBoxSize/2} ${-viewBoxSize/2} ${viewBoxSize} ${viewBoxSize}`}>
                    {/* Scale Bar - 1mm */}
                    <line x1={-viewBoxSize/2 + 1} y1={viewBoxSize/2 - 2} x2={-viewBoxSize/2 + 1 + 1} y2={viewBoxSize/2 - 2} stroke="black" strokeWidth="0.1" />
                    <text x={-viewBoxSize/2 + 1} y={viewBoxSize/2 - 2.5} fontSize="0.5" fill="black">1mm</text>

                    {/* Stator Core - Grey */}
                    <path d={fullStatorPath} fill="#9ca3af" stroke="#4b5563" strokeWidth="0.05" fillRule="evenodd" />
                    
                    {/* Rotor Back Iron (if exists) - Grey */}
                    {hasBackIron && (
                        <circle cx="0" cy="0" r={R_mag_in} fill="#9ca3af" stroke="#4b5563" strokeWidth="0.05" />
                    )}

                    {/* Magnets */}
                    {magnetPaths.map((m, i) => (
                        <path key={i} d={m.d} fill={m.color} stroke="none" />
                    ))}
                    
                    {/* Shaft - White */}
                    <circle cx="0" cy="0" r={R_shaft} fill="white" stroke="#4b5563" strokeWidth="0.05" />
                    
                    {/* Conductors */}
                    {allConductors.map((c, i) => (
                        <circle key={i} cx={c.x} cy={c.y} r={specs.winding.wire_diameter_with_insulation/2} fill="#b91c1c" stroke="none" />
                    ))}
                </svg>
            </div>
            
            {errorMsg && (
                <Alert variant="destructive" className="mt-4">
                    <AlertTriangle className="h-4 w-4" />
                    <AlertTitle>Fitting Error</AlertTitle>
                    <AlertDescription>{errorMsg}</AlertDescription>
                </Alert>
            )}
            
            <div className="mt-4 text-sm text-muted-foreground">
                <p>Stator OD: {specs.geometry.d_stator_outer}mm</p>
                <p>Conductors/Slot: {conductorsPerSlot}</p>
                <p>Wire Diameter: {specs.winding.wire_diameter_with_insulation}mm</p>
                {errorMsg && <p className="font-bold text-red-500">Max Fit: {maxFit}</p>}
            </div>
        </CardContent>
      </Card>
    </div>
  );
}
