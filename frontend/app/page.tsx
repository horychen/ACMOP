"use client"

import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { ProjectSelector } from "@/components/ProjectSelector";
import { useProject } from "@/context/ProjectContext";

export default function Home() {
  const { currentProject } = useProject();

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
    </div>
  );
}
