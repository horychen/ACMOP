"use client";

import React, { useState, useEffect, useRef, useMemo } from 'react';
import Editor from '@monaco-editor/react';
import * as d3 from 'd3';
import { useTheme } from '@/context/ThemeContext';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select';
import { Label } from '@/components/ui/label';

interface Point {
    name: string;
    x: number;
    y: number;
}

interface DrawCommand {
    type: 'line' | 'arc';
    points: number[][];
    center?: number[];
}

interface InteractiveDrawingEditorProps {
    initialCode?: string;
    parameters?: Record<string, number>;
    availableComponents?: string[];
    selectedComponent?: string;
    onComponentChange?: (component: string) => void;
    componentCodeMap?: Record<string, string>;
}

// JavaScript drawer implementation
class JSDrawer {
    private commands: DrawCommand[] = [];
    private points: Map<string, Point> = new Map();
    private pointReferences: Map<string, number[]> = new Map(); // Track point references by name

    drawLine(p1: number[] | Point, p2: number[] | Point): DrawCommand[] {
        const start = Array.isArray(p1) ? p1 : [p1.x, p1.y];
        const end = Array.isArray(p2) ? p2 : [p2.x, p2.y];
        
        // Track point references - try to identify point names from the call stack or context
        // For now, we'll track the actual point values and try to match them later
        this.trackPoint(start);
        this.trackPoint(end);
        
        const cmd: DrawCommand = {
            type: 'line',
            points: [start, end]
        };
        this.commands.push(cmd);
        return [cmd];
    }

    drawArc(center: number[], start: number[] | Point, end: number[] | Point): DrawCommand[] {
        const startPt = Array.isArray(start) ? start : [start.x, start.y];
        const endPt = Array.isArray(end) ? end : [end.x, end.y];
        
        // Track point references
        if (center) {
            this.trackPoint(center);
        }
        this.trackPoint(startPt);
        this.trackPoint(endPt);
        
        const cmd: DrawCommand = {
            type: 'arc',
            center: center,
            points: [startPt, endPt]
        };
        this.commands.push(cmd);
        return [cmd];
    }

    private trackPoint(point: number[]): void {
        // Create a key from the point coordinates for tracking
        const key = `${point[0].toFixed(6)},${point[1].toFixed(6)}`;
        this.pointReferences.set(key, point);
    }

    getSketch(name: string, color: string): void {
        // Store sketch info if needed
    }

    getCommands(): DrawCommand[] {
        return this.commands;
    }

    clear(): void {
        this.commands = [];
        this.points.clear();
        this.pointReferences.clear();
    }

    registerPoint(name: string, x: number, y: number): void {
        this.points.set(name, { name, x, y });
    }

    getPoints(): Map<string, Point> {
        return this.points;
    }
    
    // Get all unique points from commands
    getAllPointsFromCommands(): Map<string, Point> {
        const pointMap = new Map<string, Point>();
        const seenPoints = new Map<string, string>(); // Map coordinate string to point name
        
        // First, collect all points from registered points
        this.points.forEach((point, name) => {
            const coordKey = `${point.x.toFixed(6)},${point.y.toFixed(6)}`;
            seenPoints.set(coordKey, name);
            pointMap.set(name, point);
        });
        
        // Then, collect points from commands that aren't registered yet
        let pointCounter = 1;
        this.commands.forEach(cmd => {
            if (cmd.type === 'line') {
                cmd.points.forEach((pt, idx) => {
                    const coordKey = `${pt[0].toFixed(6)},${pt[1].toFixed(6)}`;
                    if (!seenPoints.has(coordKey)) {
                        const name = `P${pointCounter++}`;
                        seenPoints.set(coordKey, name);
                        pointMap.set(name, { name, x: pt[0], y: pt[1] });
                    }
                });
            } else if (cmd.type === 'arc') {
                if (cmd.center) {
                    const centerKey = `${cmd.center[0].toFixed(6)},${cmd.center[1].toFixed(6)}`;
                    if (!seenPoints.has(centerKey) && cmd.center[0] !== 0 && cmd.center[1] !== 0) {
                        const name = `C${pointCounter++}`;
                        seenPoints.set(centerKey, name);
                        pointMap.set(name, { name, x: cmd.center[0], y: cmd.center[1] });
                    }
                }
                cmd.points.forEach((pt, idx) => {
                    const coordKey = `${pt[0].toFixed(6)},${pt[1].toFixed(6)}`;
                    if (!seenPoints.has(coordKey)) {
                        const name = `P${pointCounter++}`;
                        seenPoints.set(coordKey, name);
                        pointMap.set(name, { name, x: pt[0], y: pt[1] });
                    }
                });
            }
        });
        
        return pointMap;
    }
}

// Convert Python-like code to JavaScript
function convertPythonToJS(code: string): string {
    let jsCode = code;
    
    // Replace Python math functions
    jsCode = jsCode.replace(/np\.pi/g, 'Math.PI');
    jsCode = jsCode.replace(/np\.cos/g, 'Math.cos');
    jsCode = jsCode.replace(/np\.sin/g, 'Math.sin');
    jsCode = jsCode.replace(/np\.sqrt/g, 'Math.sqrt');
    jsCode = jsCode.replace(/np\.abs/g, 'Math.abs');
    jsCode = jsCode.replace(/np\.arctan2/g, 'Math.atan2');
    jsCode = jsCode.replace(/np\.arccos/g, 'Math.acos');
    
    // Replace standalone cos/sin (avoid replacing Math.cos/sin)
    jsCode = jsCode.replace(/\bcos\(/g, 'Math.cos(');
    jsCode = jsCode.replace(/\bsin\(/g, 'Math.sin(');
    // Fix double Math.Math
    jsCode = jsCode.replace(/Math\.Math\./g, 'Math.');
    
    // Replace list syntax [x, y] stays the same in JS
    // Replace variable assignments (P1 = [x, y] becomes const P1 = [x, y])
    // But only if not already const/let/var
    jsCode = jsCode.replace(/^(\s*)(?!(?:const|let|var)\s)([A-Z]\w+)\s*=\s*\[/gm, '$1const $2 = [');
    
    // Replace intermediate variables (r_P2, alpha_P3, etc.)
    jsCode = jsCode.replace(/^(\s*)(?!(?:const|let|var)\s)([a-z_]\w+)\s*=\s*/gm, '$1const $2 = ');
    
    // Replace list_segments += drawer.drawLine(...) with drawer.drawLine(...)
    jsCode = jsCode.replace(/list_segments\s*\+=\s*drawer\./g, 'drawer.');
    
    // Remove list_segments initialization if present
    jsCode = jsCode.replace(/list_segments\s*=\s*\[\];?\s*/g, '');
    
    // Replace Python comments (but preserve // comments)
    jsCode = jsCode.replace(/(?<!\/\/\s*)#(.*)$/gm, (match, comment) => {
        return match.startsWith('//') ? match : `// ${comment}`;
    });
    
    // Replace abs() with Math.abs() (avoid Math.Math.abs)
    jsCode = jsCode.replace(/\babs\(/g, 'Math.abs(');
    jsCode = jsCode.replace(/Math\.Math\./g, 'Math.');
    
    return jsCode;
}

// Extract point definitions from code by executing a simplified version
function extractPoints(code: string, parameters: Record<string, any>): Map<string, Point> {
    const points = new Map<string, Point>();
    
    try {
        // Create a context for evaluation
        const context: any = {
            ...parameters,
            Math,
            drawer: {
                drawLine: () => {},
                drawArc: () => {},
                getSketch: () => {}
            }
        };
        
        // Convert and execute code to get point values
        const jsCode = convertPythonToJS(code);
        
        // Extract point assignments and evaluate them
        const pointRegex = /(?:const\s+)?([A-Z]\d+)\s*=\s*\[([^\]]+)\]/g;
        let match;
        
        while ((match = pointRegex.exec(code)) !== null) {
            const name = match[1];
            const expr = match[2];
            
            try {
                // Try to evaluate the expression
                const evalCode = `(${expr})`;
                const result = new Function(...Object.keys(context), `return ${evalCode}`)(...Object.values(context));
                
                if (Array.isArray(result) && result.length >= 2) {
                    const x = typeof result[0] === 'number' ? result[0] : parseFloat(result[0]);
                    const y = typeof result[1] === 'number' ? result[1] : parseFloat(result[1]);
                    
                    if (!isNaN(x) && !isNaN(y)) {
                        points.set(name, { name, x, y });
                    }
                }
            } catch (e) {
                // If evaluation fails, try parsing as literal
                const coords = expr.split(',').map(s => {
                    const trimmed = s.trim();
                    // Try to evaluate if it's an expression
                    try {
                        return new Function(...Object.keys(context), `return ${trimmed}`)(...Object.values(context));
                    } catch {
                        return parseFloat(trimmed);
                    }
                });
                
                if (coords.length >= 2 && !isNaN(coords[0]) && !isNaN(coords[1])) {
                    points.set(name, { name, x: coords[0], y: coords[1] });
                }
            }
        }
    } catch (e) {
        console.warn('Error extracting points:', e);
    }
    
    return points;
}

// Default code templates for different components
function getDefaultCodeForComponent(component: string): string {
    const templates: Record<string, string> = {
        rotorCore: `// Rotor Core - Point calculations based on CrossSectInnerNotchedRotor
// Parameters: r_ri, d_ri, d_rp, alpha_rp, alpha_rm, alpha_rs

const P1 = [r_ri, 0];

const r_P2 = r_ri + d_ri + d_rp;
const P2 = [r_P2, 0];

const alpha_P3 = alpha_rp - alpha_rm;
const P3 = [r_P2 * Math.cos(alpha_P3), r_P2 * -Math.sin(alpha_P3)];

const r_P4 = r_ri + d_ri;
const P4 = [r_P4 * Math.cos(alpha_P3), r_P4 * -Math.sin(alpha_P3)];

const alpha_P5 = alpha_P3 + alpha_rs;
const P5 = [r_P4 * Math.cos(alpha_P5), r_P4 * -Math.sin(alpha_P5)];

const P6 = [r_ri * Math.cos(alpha_P5), r_ri * -Math.sin(alpha_P5)];

// Drawing commands
drawer.drawLine(P1, P2);
drawer.drawArc([0, 0], P3, P2);
drawer.drawLine(P3, P4);
drawer.drawArc([0, 0], P5, P4);
drawer.drawLine(P5, P6);
drawer.drawArc([0, 0], P6, P1);`,
        shaft: `// Shaft - Simple circular cross-section
// Parameters: r_ri

const P1 = [r_ri, 0];
const NP1 = [-r_ri, 0];

// Drawing commands
drawer.drawArc([0, 0], NP1, P1);
drawer.drawArc([0, 0], P1, NP1);`,
        rotorMagnet: `// Rotor Magnet - Permanent magnet geometry
// Parameters: r_ri, d_ri, d_pm, alpha_rp, alpha_rm, alpha_rs

const P1 = [r_ri, 0];
const r_P2 = r_ri + d_ri + d_rp;
const P2 = [r_P2, 0];

const alpha_P3 = alpha_rp - alpha_rm;
const r_P4 = r_ri + d_ri;
const P4 = [r_P4 * Math.cos(alpha_P3), r_P4 * -Math.sin(alpha_P3)];

const P3_extra = [(r_P4 + d_pm) * Math.cos(alpha_P3), (r_P4 + d_pm) * -Math.sin(alpha_P3)];

const alpha_P5 = alpha_P3 + alpha_rs;
const P5 = [r_P4 * Math.cos(alpha_P5), r_P4 * -Math.sin(alpha_P5)];
const P6_extra = [(r_P4 + d_pm) * Math.cos(alpha_P5), (r_P4 + d_pm) * -Math.sin(alpha_P5)];

// Drawing commands
drawer.drawLine(P3_extra, P4);
drawer.drawArc([0, 0], P5, P4);
drawer.drawLine(P5, P6_extra);
drawer.drawArc([0, 0], P6_extra, P3_extra);`,
        sleeve: `// Sleeve - Outer protective layer
// Parameters: r_or, d_sleeve, p

const P1 = [r_or, 0];
const P2 = [r_or + d_sleeve, 0];

const P3 = [Math.cos(Math.PI / p) * P1[0] + Math.sin(Math.PI / p) * P1[1],
            -Math.sin(Math.PI / p) * P1[0] + Math.cos(Math.PI / p) * P1[1]];
const P4 = [Math.cos(Math.PI / p) * P2[0] + Math.sin(Math.PI / p) * P2[1],
            -Math.sin(Math.PI / p) * P2[0] + Math.cos(Math.PI / p) * P2[1]];

// Drawing commands
drawer.drawLine(P1, P2);
drawer.drawArc([0, 0], P4, P2);
drawer.drawLine(P4, P3);
drawer.drawArc([0, 0], P3, P1);`,
        statorCore: `// Stator Core - Stator iron geometry
// Parameters: r_si, r_so, Qs (number of slots)

const slotAngle = 2 * Math.PI / Qs;
const P1 = [r_si, 0];
const P2 = [r_so, 0];

// Drawing commands (simplified - one slot/tooth pair)
drawer.drawLine(P1, P2);
drawer.drawArc([0, 0], [r_so * Math.cos(slotAngle), r_so * Math.sin(slotAngle)], P2);
drawer.drawLine([r_so * Math.cos(slotAngle), r_so * Math.sin(slotAngle)], 
                [r_si * Math.cos(slotAngle), r_si * Math.sin(slotAngle)]);
drawer.drawArc([0, 0], [r_si * Math.cos(slotAngle), r_si * Math.sin(slotAngle)], P1);`,
        coils: `// Coils - Stator winding geometry
// Parameters: r_si, r_so, Qs, slotDepth

const slotAngle = 2 * Math.PI / Qs;
const r_slot_bottom = r_si + slotDepth;
const P1 = [r_slot_bottom, 0];
const P2 = [r_si, 0];

// Drawing commands (winding in one slot)
drawer.drawLine(P1, P2);
drawer.drawArc([0, 0], [r_si * Math.cos(slotAngle), r_si * Math.sin(slotAngle)], P2);
drawer.drawLine([r_si * Math.cos(slotAngle), r_si * Math.sin(slotAngle)],
                [r_slot_bottom * Math.cos(slotAngle), r_slot_bottom * Math.sin(slotAngle)]);
drawer.drawArc([0, 0], [r_slot_bottom * Math.cos(slotAngle), r_slot_bottom * Math.sin(slotAngle)], P1);`
    };
    
    return templates[component] || templates.rotorCore;
}

const InteractiveDrawingEditor: React.FC<InteractiveDrawingEditorProps> = ({ 
    initialCode,
    parameters = {
        r_ri: 40,
        d_ri: 8,
        d_rp: 5,
        alpha_rp: Math.PI / 2,
        alpha_rm: Math.PI / 3,
        alpha_rs: Math.PI / 18,
        r_or: 48,
        d_sleeve: 1,
        p: 2,
        r_si: 42,
        r_so: 50,
        Qs: 12,
        slotDepth: 5,
        d_pm: 2
    },
    availableComponents = ["rotorCore", "shaft", "rotorMagnet", "sleeve", "statorCore", "coils"],
    selectedComponent: externalSelectedComponent,
    onComponentChange,
    componentCodeMap
}) => {
    const { theme } = useTheme();
    const svgRef = useRef<SVGSVGElement>(null);
    const [internalSelectedComponent, setInternalSelectedComponent] = useState(availableComponents[0] || "rotorCore");
    
    // Use external selectedComponent if provided, otherwise use internal state
    const selectedComponent = externalSelectedComponent ?? internalSelectedComponent;
    
    // Get initial code for selected component
    const getCodeForComponent = (component: string): string => {
        if (componentCodeMap && componentCodeMap[component]) {
            return componentCodeMap[component];
        }
        if (initialCode && component === availableComponents[0]) {
            return initialCode;
        }
        return getDefaultCodeForComponent(component);
    };
    
    const [code, setCode] = useState(() => getCodeForComponent(selectedComponent));
    const [error, setError] = useState<string | null>(null);
    const [drawer] = useState(() => new JSDrawer());
    const [points, setPoints] = useState<Map<string, Point>>(new Map());
    const [zoomLevel, setZoomLevel] = useState(1);
    const [panX, setPanX] = useState(0);
    const [panY, setPanY] = useState(0);
    const containerRef = useRef<HTMLDivElement>(null);
    
    // Update code when component changes
    useEffect(() => {
        const newCode = getCodeForComponent(selectedComponent);
        setCode(newCode);
    }, [selectedComponent, componentCodeMap, initialCode, availableComponents]);
    
    // Handle wheel event to prevent page scrolling
    useEffect(() => {
        const container = containerRef.current;
        if (!container) return;
        
        const handleWheel = (e: WheelEvent) => {
            e.preventDefault();
            e.stopPropagation();
            const delta = e.deltaY > 0 ? 0.9 : 1.1;
            setZoomLevel(prev => Math.max(0.1, Math.min(10, prev * delta)));
        };
        
        // Use capture phase and non-passive listener to ensure preventDefault works
        container.addEventListener('wheel', handleWheel, { passive: false, capture: true });
        
        return () => {
            container.removeEventListener('wheel', handleWheel, { capture: true } as EventListenerOptions);
        };
    }, []);
    
    // Memoize parameters to avoid unnecessary re-renders
    const paramsString = useMemo(() => JSON.stringify(parameters || {}), [parameters]);
    
    const handleComponentChange = (component: string) => {
        if (onComponentChange) {
            onComponentChange(component);
        } else {
            setInternalSelectedComponent(component);
        }
    };

    // Execute code and render
    useEffect(() => {
        if (!svgRef.current) return;

        drawer.clear();
        setError(null);
        setPoints(new Map());

        try {
            // Convert Python-like code to JavaScript
            const jsCode = convertPythonToJS(code);
            
            // Create execution context - ensure all values are serializable
            const context: Record<string, any> = {
                drawer,
                Math,
                console
            };
            
            // Add parameters safely
            if (parameters) {
                Object.entries(parameters).forEach(([key, value]) => {
                    // Only add primitive values and arrays
                    if (typeof value === 'number' || typeof value === 'string' || typeof value === 'boolean' || Array.isArray(value)) {
                        context[key] = value;
                    }
                });
            }

            // Execute the code - build function with proper parameter names
            // Validate parameter names are valid JavaScript identifiers
            const paramNames = Object.keys(context).filter(name => {
                // Check if it's a valid JavaScript identifier
                return /^[a-zA-Z_$][a-zA-Z0-9_$]*$/.test(name);
            });
            
            if (paramNames.length !== Object.keys(context).length) {
                throw new Error('Some parameter names are invalid JavaScript identifiers');
            }
            
            const paramValues = paramNames.map(key => context[key]);
            
            // Create function with parameter names as strings
            // Wrap in try-catch for better error handling
            let func: Function;
            try {
                func = new Function(...paramNames, jsCode);
            } catch (funcError) {
                throw new Error(`Failed to create function: ${funcError instanceof Error ? funcError.message : 'Unknown error'}. Check your code syntax.`);
            }
            
            // First, extract all point definitions from code and calculate their values
            // This gives us a map of point names to coordinates
            // Improved regex to match various point naming patterns: P1, P2, P3_extra, NP1, C1, etc.
            const pointRegex = /(?:const\s+|let\s+|var\s+)?([A-Z][A-Za-z0-9_]*)\s*=\s*\[([^\]]+)\]/g;
            let match;
            const pointValueMap = new Map<string, number[]>();
            
            // Reset regex lastIndex
            pointRegex.lastIndex = 0;
            while ((match = pointRegex.exec(code)) !== null) {
                const name = match[1];
                const expr = match[2].trim();
                // Only process if it looks like a point (starts with P, NP, C, or similar uppercase letter followed by alphanumeric/underscore)
                // This matches: P1, P2, P3_extra, NP1, C1, etc.
                if (name.match(/^[A-Z][A-Za-z0-9_]*$/) && (name.startsWith('P') || name.startsWith('NP') || name.startsWith('C'))) {
                    try {
                        // Split expression by comma and evaluate each part
                        const parts = expr.split(',').map(p => p.trim());
                        if (parts.length >= 2) {
                            // Evaluate each coordinate
                            const xFunc = new Function(...paramNames, `return ${parts[0]}`);
                            const yFunc = new Function(...paramNames, `return ${parts[1]}`);
                            const x = xFunc(...paramValues);
                            const y = yFunc(...paramValues);
                            
                            if (typeof x === 'number' && typeof y === 'number' && !isNaN(x) && !isNaN(y)) {
                                pointValueMap.set(name, [x, y]);
                                console.debug(`Extracted point ${name} = [${x.toFixed(6)}, ${y.toFixed(6)}]`);
                            }
                        }
                    } catch (e) {
                        // Skip if evaluation fails
                        console.debug(`Could not evaluate point ${name}:`, e);
                    }
                }
            }
            
            // Execute the code - this will call drawer.drawLine and drawer.drawArc
            func(...paramValues);

            // Extract points from drawer commands - these are the points actually used in drawing
            const pointsFromCommands = drawer.getAllPointsFromCommands();
            
            // Match point names with coordinates from commands
            const finalPoints = new Map<string, Point>();
            const coordToName = new Map<string, string>();
            // Increased tolerance for coordinate matching (handles floating point precision issues)
            const tolerance = 1e-3; // Coordinate matching tolerance
            
            // First, map coordinates to names from point definitions
            pointValueMap.forEach((coords, name) => {
                const coordKey = `${coords[0].toFixed(6)},${coords[1].toFixed(6)}`;
                coordToName.set(coordKey, name);
            });
            
            // Then, match points from commands with names
            pointsFromCommands.forEach((point, defaultName) => {
                // Try to find matching point name by coordinates
                let matchedName = defaultName;
                let minDistance = Infinity;
                
                pointValueMap.forEach((coords, name) => {
                    const dx = Math.abs(coords[0] - point.x);
                    const dy = Math.abs(coords[1] - point.y);
                    const distance = Math.sqrt(dx * dx + dy * dy);
                    
                    if (distance < tolerance && distance < minDistance) {
                        minDistance = distance;
                        matchedName = name;
                    }
                });
                
                // Always add points from commands (they are used in drawLine/drawArc)
                finalPoints.set(matchedName, { name: matchedName, x: point.x, y: point.y });
                console.debug(`Matched point: ${matchedName} at [${point.x.toFixed(6)}, ${point.y.toFixed(6)}]`);
            });
            
            // Also add any points from pointValueMap that weren't matched (in case they're defined but not used in drawing)
            // But only if they're actually used in the code (we check if they appear in drawLine/drawArc calls)
            const codeLower = code.toLowerCase();
            pointValueMap.forEach((coords, name) => {
                // Check if this point name appears in drawing commands
                if (codeLower.includes(name.toLowerCase() + ',') || codeLower.includes(name.toLowerCase() + ')')) {
                    // Only add if not already in finalPoints
                    if (!finalPoints.has(name)) {
                        finalPoints.set(name, { name, x: coords[0], y: coords[1] });
                        console.debug(`Added unused point: ${name} at [${coords[0].toFixed(6)}, ${coords[1].toFixed(6)}]`);
                    }
                }
            });
            
            // Register points with drawer
            finalPoints.forEach((point, name) => {
                drawer.registerPoint(name, point.x, point.y);
            });
            
            // Update points state for table display
            setPoints(new Map(finalPoints));
            
            // Debug: log extracted points
            console.log(`Extracted ${finalPoints.size} points from draw commands:`, 
                Array.from(finalPoints.entries()).map(([name, p]) => `${name}(${p.x.toFixed(2)}, ${p.y.toFixed(2)})`).join(', '));

            // Render visualization
            renderVisualization(drawer, svgRef.current, theme === 'dark', zoomLevel, panX, panY);
        } catch (err) {
            setError(err instanceof Error ? err.message : 'Unknown error');
            console.error('Execution error:', err);
        }
    }, [code, theme, paramsString, drawer, zoomLevel, panX, panY]);

    return (
        <Card>
            <CardHeader>
                <div className="flex items-center justify-between">
                    <div>
                        <CardTitle>Interactive Drawing Editor</CardTitle>
                        <CardDescription>Edit Python-like code to draw geometry. Points, lines, and arcs are rendered in real-time.</CardDescription>
                    </div>
                    {availableComponents.length > 0 && (
                        <div className="flex items-center gap-2">
                            <Label htmlFor="component-select" className="text-sm font-medium">Component:</Label>
                            <Select value={selectedComponent} onValueChange={handleComponentChange}>
                                <SelectTrigger id="component-select" className="w-[180px]">
                                    <SelectValue placeholder="Select component" />
                                </SelectTrigger>
                                <SelectContent>
                                    {availableComponents.map((component) => (
                                        <SelectItem key={component} value={component}>
                                            {component}
                                        </SelectItem>
                                    ))}
                                </SelectContent>
                            </Select>
                        </div>
                    )}
                </div>
            </CardHeader>
            <CardContent>
                <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
                    {/* Code Editor */}
                    <div className="space-y-2">
                        <label className="text-sm font-medium">Code Editor - {selectedComponent}</label>
                        <div className="border rounded-lg overflow-hidden" style={{ height: '600px' }}>
                            <Editor
                                height="600px"
                                defaultLanguage="python"
                                value={code}
                                onChange={(value) => setCode(value || '')}
                                theme={theme === 'dark' ? 'vs-dark' : 'light'}
                                options={{
                                    minimap: { enabled: false },
                                    fontSize: 14,
                                    wordWrap: 'on',
                                    automaticLayout: true,
                                }}
                            />
                        </div>
                        {error && (
                            <div className="text-sm text-red-500 bg-red-50 dark:bg-red-950/30 p-2 rounded">
                                Error: {error}
                            </div>
                        )}
                    </div>

                    {/* Visualization */}
                    <div className="space-y-2">
                        <label className="text-sm font-medium">Visualization</label>
                        <div 
                            ref={containerRef}
                            className="border rounded-lg bg-background p-4 relative" 
                            style={{ height: '600px' }}
                        >
                            <svg
                                ref={svgRef}
                                width="100%"
                                height="100%"
                                className="w-full h-full"
                            />
                        </div>
                        
                        {/* Points Coordinates Table - Below Visualization */}
                        <div className="border rounded-lg bg-background p-4">
                            <h3 className="text-sm font-semibold mb-3 text-foreground">点坐标 (mm)</h3>
                            {points.size > 0 ? (
                                <div className="overflow-x-auto max-h-64 overflow-y-auto">
                                    <table className="w-full text-sm border-collapse">
                                        <thead className="sticky top-0 bg-background">
                                            <tr className="border-b">
                                                <th className="text-left p-2 font-semibold text-foreground">点名称</th>
                                                <th className="text-right p-2 font-semibold text-foreground">X (mm)</th>
                                                <th className="text-right p-2 font-semibold text-foreground">Y (mm)</th>
                                            </tr>
                                        </thead>
                                        <tbody>
                                            {Array.from(points.values())
                                                .sort((a, b) => {
                                                    // Sort by point name (P1, P2, P3, etc.)
                                                    const numA = parseInt(a.name.replace(/\D/g, '')) || 0;
                                                    const numB = parseInt(b.name.replace(/\D/g, '')) || 0;
                                                    if (numA !== numB) return numA - numB;
                                                    return a.name.localeCompare(b.name);
                                                })
                                                .map((point) => {
                                                    // Parameters are already in millimeters, no conversion needed
                                                    const x_mm = point.x.toFixed(3);
                                                    const y_mm = point.y.toFixed(3);
                                                    return (
                                                        <tr key={point.name} className="border-b hover:bg-muted/50">
                                                            <td className="p-2 font-mono font-medium text-foreground">{point.name}</td>
                                                            <td className="p-2 text-right font-mono text-foreground">{x_mm}</td>
                                                            <td className="p-2 text-right font-mono text-foreground">{y_mm}</td>
                                                        </tr>
                                                    );
                                                })}
                                        </tbody>
                                    </table>
                                </div>
                            ) : (
                                <div className="text-sm text-muted-foreground text-center py-8">
                                    暂无点坐标数据
                                </div>
                            )}
                        </div>
                    </div>
                </div>
            </CardContent>
        </Card>
    );
};

function renderVisualization(drawer: JSDrawer, svgElement: SVGSVGElement, isDark: boolean, zoomLevel: number = 1, panX: number = 0, panY: number = 0) {
    const svg = d3.select(svgElement);
    svg.selectAll("*").remove();

    const width = svgElement.clientWidth || 500;
    const height = svgElement.clientHeight || 500;
    const centerX = width / 2;
    const centerY = height / 2;

    // Calculate bounds
    const commands = drawer.getCommands();
    const points = drawer.getPoints();
    
    let minX = Infinity, maxX = -Infinity;
    let minY = Infinity, maxY = -Infinity;

    // Find bounds from commands
    commands.forEach(cmd => {
        if (cmd.type === 'line') {
            cmd.points.forEach(([x, y]) => {
                minX = Math.min(minX, x);
                maxX = Math.max(maxX, x);
                minY = Math.min(minY, y);
                maxY = Math.max(maxY, y);
            });
        } else if (cmd.type === 'arc') {
            if (cmd.center) {
                const [cx, cy] = cmd.center;
                cmd.points.forEach(([x, y]) => {
                    const dx = x - cx;
                    const dy = y - cy;
                    const r = Math.sqrt(dx * dx + dy * dy);
                    minX = Math.min(minX, cx - r, x);
                    maxX = Math.max(maxX, cx + r, x);
                    minY = Math.min(minY, cy - r, y);
                    maxY = Math.max(maxY, cy + r, y);
                });
            }
        }
    });

    // Find bounds from points
    points.forEach(point => {
        minX = Math.min(minX, point.x);
        maxX = Math.max(maxX, point.x);
        minY = Math.min(minY, point.y);
        maxY = Math.max(maxY, point.y);
    });

    if (!isFinite(minX) || !isFinite(maxX)) {
        minX = -50;
        maxX = 50;
        minY = -50;
        maxY = 50;
    }

    // Add padding
    const padding = 50;
    const rangeX = maxX - minX || 100;
    const rangeY = maxY - minY || 100;
    const baseScale = Math.min((width - padding * 2) / rangeX, (height - padding * 2) / rangeY);
    const scale = baseScale * zoomLevel;

    // Create transform with zoom and pan
    const g = svg.append("g")
        .attr("transform", `translate(${centerX + panX}, ${centerY + panY}) scale(${scale}) translate(${-(minX + maxX) / 2}, ${-(minY + maxY) / 2})`);

    // Colors
    const colors = isDark ? {
        line: '#94a3b8',
        arc: '#3b82f6',
        point: '#ef4444',
        pointLabel: '#fbbf24',
        background: '#0f172a'
    } : {
        line: '#475569',
        arc: '#2563eb',
        point: '#dc2626',
        pointLabel: '#d97706',
        background: '#f8fafc'
    };

    // Draw commands
    commands.forEach((cmd, idx) => {
        if (cmd.type === 'line') {
            const [start, end] = cmd.points;
            g.append("line")
                .attr("x1", start[0])
                .attr("y1", start[1])
                .attr("x2", end[0])
                .attr("y2", end[1])
                .attr("stroke", colors.line)
                .attr("stroke-width", 2 / scale)
                .attr("marker-end", "url(#arrowhead)");
        } else if (cmd.type === 'arc' && cmd.center) {
            const [cx, cy] = cmd.center;
            const [start, end] = cmd.points;
            
            // Calculate arc parameters
            const dx1 = start[0] - cx;
            const dy1 = start[1] - cy;
            const dx2 = end[0] - cx;
            const dy2 = end[1] - cy;
            
            const r1 = Math.sqrt(dx1 * dx1 + dy1 * dy1);
            const r2 = Math.sqrt(dx2 * dx2 + dy2 * dy2);
            const r = Math.max(r1, r2);
            
            const angle1 = Math.atan2(dy1, dx1);
            const angle2 = Math.atan2(dy2, dx2);
            
            let sweepFlag = 1;
            let largeArc = 0;
            let deltaAngle = angle2 - angle1;
            
            // Normalize angle
            if (deltaAngle > Math.PI) deltaAngle -= 2 * Math.PI;
            if (deltaAngle < -Math.PI) deltaAngle += 2 * Math.PI;
            
            if (Math.abs(deltaAngle) > Math.PI) {
                largeArc = 1;
            }
            if (deltaAngle < 0) {
                sweepFlag = 0;
            }

            g.append("path")
                .attr("d", `M ${start[0]} ${start[1]} A ${r} ${r} 0 ${largeArc} ${sweepFlag} ${end[0]} ${end[1]}`)
                .attr("fill", "none")
                .attr("stroke", colors.arc)
                .attr("stroke-width", 2 / scale);
        }
    });

    // Draw points with enhanced labels
    points.forEach(point => {
        const labelOffset = 10 / scale;
        const pointRadius = 4 / scale;
        
        // Point circle - larger and more visible
        g.append("circle")
            .attr("cx", point.x)
            .attr("cy", point.y)
            .attr("r", pointRadius)
            .attr("fill", colors.point)
            .attr("stroke", colors.pointLabel)
            .attr("stroke-width", 2 / scale)
            .attr("opacity", 0.9);

        // Outer ring for better visibility
        g.append("circle")
            .attr("cx", point.x)
            .attr("cy", point.y)
            .attr("r", pointRadius * 1.5)
            .attr("fill", "none")
            .attr("stroke", colors.point)
            .attr("stroke-width", 1 / scale)
            .attr("opacity", 0.5);

        // Calculate label position (offset to avoid overlap)
        const labelX = point.x + labelOffset;
        const labelY = point.y - labelOffset;
        
        // Create a group for the label to manage z-order
        const labelGroup = g.append("g");
        
        // Create temporary text to measure bounding box
        const tempText = labelGroup.append("text")
            .attr("x", labelX)
            .attr("y", labelY)
            .attr("font-size", `${14 / scale}px`)
            .attr("font-weight", "bold")
            .attr("font-family", "monospace")
            .attr("opacity", 0)
            .text(point.name);
        
        // Get text bounding box for background
        const textBBox = (tempText.node() as SVGTextElement)?.getBBox();
        tempText.remove();
        
        // Point name text (no background box, with stroke for visibility)
        labelGroup.append("text")
            .attr("x", labelX)
            .attr("y", labelY)
            .attr("fill", colors.pointLabel)
            .attr("font-size", `${14 / scale}px`)
            .attr("font-weight", "bold")
            .attr("font-family", "monospace")
            .attr("text-anchor", "middle")
            .attr("dominant-baseline", "middle")
            .attr("stroke", isDark ? "rgba(15, 23, 42, 0.8)" : "rgba(248, 250, 252, 0.8)")
            .attr("stroke-width", 3 / scale)
            .attr("paint-order", "stroke")
            .text(point.name);
        
        // Connecting line from point to label
        g.append("line")
            .attr("x1", point.x + pointRadius)
            .attr("y1", point.y - pointRadius)
            .attr("x2", labelX)
            .attr("y2", labelY)
            .attr("stroke", colors.pointLabel)
            .attr("stroke-width", 0.5 / scale)
            .attr("stroke-dasharray", `${2 / scale},${2 / scale}`)
            .attr("opacity", 0.6);
    });

    // Draw scale bar (in screen coordinates, not transformed)
    const scaleBarLengthPx = 100; // 100 pixels
    const scaleBarLengthModel = scaleBarLengthPx / scale; // Convert to model coordinates
    const scaleBarX = 20;
    const scaleBarY = height - 40;
    
    // Scale bar line
    svg.append("line")
        .attr("x1", scaleBarX)
        .attr("y1", scaleBarY)
        .attr("x2", scaleBarX + scaleBarLengthPx)
        .attr("y2", scaleBarY)
        .attr("stroke", colors.pointLabel)
        .attr("stroke-width", 2)
        .attr("marker-end", "url(#arrowhead-scale)");
    
    // Scale bar label
    svg.append("text")
        .attr("x", scaleBarX + scaleBarLengthPx / 2)
        .attr("y", scaleBarY - 8)
        .attr("fill", colors.pointLabel)
        .attr("font-size", "12px")
        .attr("font-weight", "bold")
        .attr("font-family", "monospace")
        .attr("text-anchor", "middle")
        .attr("stroke", isDark ? "rgba(15, 23, 42, 0.8)" : "rgba(248, 250, 252, 0.8)")
        .attr("stroke-width", 2)
        .attr("paint-order", "stroke")
        .text(`${scaleBarLengthModel.toFixed(1)} mm`);

    // Add arrow markers
    const defs = svg.append("defs");
    
    // Arrow marker for lines
    defs.append("marker")
        .attr("id", "arrowhead")
        .attr("markerWidth", 10)
        .attr("markerHeight", 10)
        .attr("refX", 9)
        .attr("refY", 3)
        .attr("orient", "auto")
        .append("polygon")
        .attr("points", "0 0, 10 3, 0 6")
        .attr("fill", colors.line);
    
    // Arrow marker for scale bar
    defs.append("marker")
        .attr("id", "arrowhead-scale")
        .attr("markerWidth", 10)
        .attr("markerHeight", 10)
        .attr("refX", 9)
        .attr("refY", 3)
        .attr("orient", "auto")
        .append("polygon")
        .attr("points", "0 0, 10 3, 0 6")
        .attr("fill", colors.pointLabel);
}

export default InteractiveDrawingEditor;

