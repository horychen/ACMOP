"use client";

import React, { useEffect, useState, useRef } from "react";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Loader2, AlertTriangle, RefreshCw } from "lucide-react";
import { Button } from "@/components/ui/button";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/alert";

export default function AnimatedSvgVisualizer() {
    const [svgContent, setSvgContent] = useState<string | null>(null);
    const [loading, setLoading] = useState(true);
    const [error, setError] = useState<string | null>(null);
    const [spinSpeed, setSpinSpeed] = useState(4); // seconds per revolution
    const svgContainerRef = useRef<HTMLDivElement>(null);

    const svgUrl = `http://localhost:8001/api/stator-svg?t=${Date.now()}`;

    const fetchSvg = () => {
        setLoading(true);
        setError(null);
        fetch(svgUrl)
            .then(res => {
                const contentType = res.headers.get("content-type");
                if (!res.ok) {
                    throw new Error(`Failed to fetch SVG (${res.status})`);
                }
                if (contentType && contentType.includes("application/json")) {
                    return res.json().then(data => {
                        throw new Error(data.error || "Failed to fetch SVG");
                    });
                }
                return res.text();
            })
            .then(text => {
                setSvgContent(text);
                setLoading(false);
            })
            .catch(err => {
                setError(err.message);
                setLoading(false);
            });
    };

    useEffect(() => {
        fetchSvg();
    }, []);

    useEffect(() => {
        if (!svgContent || !svgContainerRef.current) return;

        // Simple SVG parsing and transformation 
        // We will find all <path> tags, and put them in stator or rotor groups based on stroke color
        const parser = new DOMParser();
        const doc = parser.parseFromString(svgContent, "image/svg+xml");
        const rootSvg = doc.querySelector("svg");

        if (!rootSvg) return;

        // Grouping paths
        const statorGroup = document.createElementNS("http://www.w3.org/2000/svg", "g");
        statorGroup.setAttribute("class", "stator-group");

        const rotorGroup = document.createElementNS("http://www.w3.org/2000/svg", "g");
        rotorGroup.setAttribute("class", "rotor-group");
        rotorGroup.style.transformOrigin = "250px 250px";
        rotorGroup.style.animation = `spin-cw ${spinSpeed}s linear infinite`;

        const defs = document.createElementNS("http://www.w3.org/2000/svg", "defs");
        defs.innerHTML = `
            <filter id="glow" x="-20%" y="-20%" width="140%" height="140%">
                <feGaussianBlur stdDeviation="2" result="blur"/>
                <feMerge>
                    <feMergeNode in="blur"/>
                    <feMergeNode in="SourceGraphic"/>
                </feMerge>
            </filter>
        `;
        rootSvg.insertBefore(defs, rootSvg.firstChild);

        // Add a style tag for the animation heartbeat and other overrides
        const style = document.createElementNS("http://www.w3.org/2000/svg", "style");
        style.innerHTML = `
            @keyframes spin-cw {
                from { transform: rotate(0deg); }
                to { transform: rotate(360deg); }
            }
            .pulse {
                animation: pulse-op 2s ease-in-out infinite alternate;
            }
            @keyframes pulse-op {
                0% { stroke-width: 0.05; stroke-opacity: 0.8; }
                100% { stroke-width: 0.15; stroke-opacity: 1; }
            }
        `;
        rootSvg.insertBefore(style, rootSvg.firstChild);

        // Classify paths into Stator or Rotor
        const elementsToMove = Array.from(rootSvg.childNodes);
        elementsToMove.forEach((node) => {
            if (node.nodeType !== Node.ELEMENT_NODE) return;
            const el = node as Element;

            // Background rect or other non-path shapes we can leave at root or put in stator
            if (el.tagName === "rect" && el.getAttribute("width") === "600") {
                // Background
                el.setAttribute("fill", "#0f172a"); // Dark mode background
                statorGroup.appendChild(node);
                return;
            }

            if (el.tagName === "path") {
                const stroke = el.getAttribute("stroke") || "";

                // Adjust stroke width for better visibility
                el.setAttribute("stroke-width", "0.05");
                el.setAttribute("filter", "url(#glow)");

                if (stroke.includes("rgb(13.3") || stroke.includes("rgb(33.3")) {
                    // Rotor magnet / shaft
                    if (stroke.includes("rgb(13.3")) {
                        el.setAttribute("stroke", "#06b6d4"); // Cyan
                    } else {
                        el.setAttribute("stroke", "#94a3b8"); // Gray shaft
                    }
                    rotorGroup.appendChild(node);
                } else if (stroke.includes("rgb(40") || stroke.includes("rgb(72.1")) {
                    // Stator core / winding
                    if (stroke.includes("rgb(40")) {
                        el.setAttribute("stroke", "#64748b"); // Slate
                    } else {
                        el.setAttribute("stroke", "#f59e0b"); // Orange/Amber
                        el.setAttribute("class", "pulse");
                    }
                    statorGroup.appendChild(node);
                } else {
                    // Default to stator
                    statorGroup.appendChild(node);
                }
            } else if (el.tagName !== "defs" && el.tagName !== "style") {
                statorGroup.appendChild(node);
            }
        });

        // Clear existing children that were moved, and append our groups
        // We already moved them via appendChild, so they are detached from rootSvg
        rootSvg.appendChild(statorGroup);
        rootSvg.appendChild(rotorGroup);

        svgContainerRef.current.innerHTML = "";
        svgContainerRef.current.appendChild(rootSvg);

    }, [svgContent, spinSpeed]);

    return (
        <Card className="w-full">
            <CardHeader>
                <CardTitle className="flex items-center justify-between text-base">
                    Animated Machine Geometry Visualization
                    <div className="flex items-center gap-4">
                        <div className="flex items-center gap-2 text-sm font-normal">
                            <label htmlFor="speed" className="text-muted-foreground">Speed (s/rev):</label>
                            <input
                                id="speed"
                                type="range"
                                min="0.5"
                                max="10"
                                step="0.5"
                                value={spinSpeed}
                                onChange={(e) => setSpinSpeed(Number(e.target.value))}
                                className="w-24"
                            />
                            <span className="w-8 text-right font-mono">{spinSpeed}s</span>
                        </div>
                        <Button variant="outline" size="icon" className="h-8 w-8" onClick={fetchSvg} title="Reload SVG">
                            <RefreshCw className="h-4 w-4" />
                        </Button>
                    </div>
                </CardTitle>
            </CardHeader>
            <CardContent>
                {loading ? (
                    <div className="flex justify-center items-center h-[500px]">
                        <Loader2 className="h-8 w-8 animate-spin text-muted-foreground" />
                    </div>
                ) : error ? (
                    <Alert variant="destructive">
                        <AlertTriangle className="h-4 w-4" />
                        <AlertTitle>Error Loading SVG</AlertTitle>
                        <AlertDescription>{error}</AlertDescription>
                    </Alert>
                ) : (
                    <div
                        className="w-full max-w-[600px] mx-auto aspect-square relative rounded-xl overflow-hidden shadow-inner border border-slate-800"
                        style={{ backgroundColor: "#0f172a" }}
                        ref={svgContainerRef}
                    />
                )}
            </CardContent>
        </Card>
    );
}
