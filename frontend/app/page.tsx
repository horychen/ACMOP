"use client"

import { useState } from "react";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { ProjectSelector } from "@/components/ProjectSelector";
import { useProject } from "@/context/ProjectContext";
import MotorVisualizer from "@/components/MotorVisualizer";
import LinearMachineView from "@/components/LinearMachineView";
import { calculateMachineDesign } from "@/services/physicsEngine";
import { DesignSpecs } from "@/types";

const DEFAULT_SPECS: DesignSpecs = {
    ratedPower: 5, // 5 kW
    ratedSpeed: 3000, // 3000 RPM
    ratedVoltage: 400,
    outerDiameterLimit: 120,
    axialLengthLimit: 100,
    airGap: 1.0,
    slotCount: 12,
    poleCount: 4,
    currentDensity: 5
};

export default function Home() {
  const { currentProject } = useProject();
  const [viewMode, setViewMode] = useState<'circular' | 'linear'>('circular');

  // Calculate default geometry for visualization
  const defaultResult = calculateMachineDesign(DEFAULT_SPECS);

  return (
    <div className="space-y-6">
      <div className="flex items-center justify-between">
        <h1 className="text-3xl font-bold tracking-tight">Dashboard</h1>
        <ProjectSelector />
      </div>

      <div className="grid gap-4 md:grid-cols-2 lg:grid-cols-3">
        <Card>
          <CardHeader>
            <CardTitle>Current Project</CardTitle>
            <CardDescription>Active project selection</CardDescription>
          </CardHeader>
          <CardContent>
            <p className="text-lg font-semibold">{currentProject || "No project selected"}</p>
          </CardContent>
        </Card>
        <Card>
          <CardHeader>
            <CardTitle>Design Viewer</CardTitle>
            <CardDescription>View and modify machine geometry</CardDescription>
          </CardHeader>
          <CardContent>
            <p>Visualize the cross-section and winding layout.</p>
          </CardContent>
        </Card>
        <Card>
          <CardHeader>
            <CardTitle>Optimization</CardTitle>
            <CardDescription>Run and analyze optimizations</CardDescription>
          </CardHeader>
          <CardContent>
            <p>View Pareto fronts and select optimal designs.</p>
          </CardContent>
        </Card>
      </div>

      {/* Machine Visualizer Section */}
      <Card className="col-span-full">
        <CardHeader>
          <div className="flex items-center justify-between">
            <div>
              <CardTitle>Machine Geometry Visualization</CardTitle>
              <CardDescription>Interactive view of motor cross-section and linear layout</CardDescription>
            </div>
            <div className="flex items-center space-x-2">
              <button
                onClick={() => setViewMode('circular')}
                className={`text-xs px-3 py-1 rounded transition-colors ${
                  viewMode === 'circular'
                    ? 'bg-primary text-primary-foreground'
                    : 'bg-muted text-muted-foreground hover:bg-muted/80'
                }`}
              >
                Circular
              </button>
              <button
                onClick={() => setViewMode('linear')}
                className={`text-xs px-3 py-1 rounded transition-colors ${
                  viewMode === 'linear'
                    ? 'bg-primary text-primary-foreground'
                    : 'bg-muted text-muted-foreground hover:bg-muted/80'
                }`}
              >
                Linear
              </button>
            </div>
          </div>
        </CardHeader>
        <CardContent>
          <div className="w-full">
            {viewMode === 'circular' 
              ? <MotorVisualizer geometry={defaultResult.geometry} />
              : <LinearMachineView geometry={defaultResult.geometry} />
            }
          </div>
        </CardContent>
      </Card>
    </div>
  );
}
export const dynamic = 'force-dynamic'
