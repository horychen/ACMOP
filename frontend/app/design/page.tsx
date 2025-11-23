"use client"

import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card"
import { ProjectSelector } from "@/components/ProjectSelector"
import { useProject } from "@/context/ProjectContext"

export default function DesignViewerPage() {
    const { currentProject } = useProject()

    return (
        <div className="space-y-6">
            <div className="flex items-center justify-between">
                <h1 className="text-3xl font-bold tracking-tight">Design Viewer</h1>
                <ProjectSelector />
            </div>

            <div className="grid gap-4 md:grid-cols-2">
                <Card className="col-span-2 md:col-span-1">
                    <CardHeader>
                        <CardTitle>Geometry: {currentProject}</CardTitle>
                    </CardHeader>
                    <CardContent>
                        <div className="aspect-square bg-muted rounded-md flex items-center justify-center">
                            <span className="text-muted-foreground">Geometry Visualization Placeholder</span>
                        </div>
                    </CardContent>
                </Card>
                <Card className="col-span-2 md:col-span-1">
                    <CardHeader>
                        <CardTitle>Parameters</CardTitle>
                    </CardHeader>
                    <CardContent>
                        <div className="space-y-2">
                            <div className="flex justify-between">
                                <span className="font-medium">Stator Outer Radius:</span>
                                <span>100 mm</span>
                            </div>
                            <div className="flex justify-between">
                                <span className="font-medium">Rotor Outer Radius:</span>
                                <span>60 mm</span>
                            </div>
                            {/* Add more parameters here */}
                        </div>
                    </CardContent>
                </Card>
            </div>
        </div>
    )
}
