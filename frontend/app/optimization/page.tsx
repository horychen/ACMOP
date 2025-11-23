"use client"

import { useState, useEffect } from "react"
import { Button } from "@/components/ui/button"
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card"
import { ProjectSelector } from "@/components/ProjectSelector"
import { ParetoChart } from "@/components/ParetoChart"
import { useProject } from "@/context/ProjectContext"
import axios from "axios"

export default function OptimizationPage() {
    const { currentProject } = useProject()
    const [data, setData] = useState<any[]>([])
    const [loading, setLoading] = useState(false)

    const fetchData = async (project: string) => {
        if (!project) return
        setLoading(true)
        try {
            const res = await axios.get(`http://localhost:8000/api/results/swarm/${project}`)
            const rawData = res.data
            console.log("Fetched data:", rawData)

            // Process data: Extract Performance metrics
            // Structure: { "Initial": { "DesignName": { "Performance": { ... } } } }
            // Or sometimes just { "DesignName": { "Performance": { ... } } } if "Initial" key is missing/different

            let designs = []
            // Handle the "Initial" wrapper if present
            const root = rawData["Initial"] ? rawData["Initial"] : rawData

            for (const key in root) {
                const design = root[key]
                if (design["Performance"]) {
                    const perf = design["Performance"]
                    // Try to find torque ripple and efficiency
                    // Adjust keys based on actual JSON content
                    designs.push({
                        name: key,
                        x: perf["normalized_torque_ripple"] || perf["TorqueRipple"] || perf["f1"] || 0,
                        y: perf["efficiency"] || perf["Efficiency"] || perf["RatedEfficiency"] || perf["f2"] || 0,
                        ...perf
                    })
                }
            }
            setData(designs)
        } catch (err) {
            console.error("Failed to fetch swarm data", err)
            setData([])
        } finally {
            setLoading(false)
        }
    }

    useEffect(() => {
        if (currentProject) {
            fetchData(currentProject)
        }
    }, [currentProject])

    return (
        <div className="space-y-6">
            <div className="flex items-center justify-between">
                <h1 className="text-3xl font-bold tracking-tight">Optimization</h1>
                <div className="flex items-center gap-4">
                    <ProjectSelector />
                    <Button onClick={() => fetchData(currentProject)} disabled={loading}>
                        {loading ? "Loading..." : "Refresh Data"}
                    </Button>
                </div>
            </div>

            <div className="grid gap-4 md:grid-cols-2">
                <Card className="col-span-2">
                    <CardHeader>
                        <CardTitle>Pareto Front: {currentProject}</CardTitle>
                    </CardHeader>
                    <CardContent>
                        {data.length > 0 ? (
                            <ParetoChart
                                data={data}
                                xKey="x"
                                yKey="y"
                                xLabel="Torque Ripple / f1"
                                yLabel="Efficiency / f2"
                            />
                        ) : (
                            <div className="flex h-[400px] items-center justify-center border rounded-md">
                                <p className="text-muted-foreground">No data available or failed to load.</p>
                            </div>
                        )}
                    </CardContent>
                </Card>
            </div>
        </div>
    )
}
