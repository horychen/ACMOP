"use client"

import { useState, useEffect } from "react"
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card"
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs"
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "@/components/ui/table"
import { Badge } from "@/components/ui/badge"
import { ScrollArea } from "@/components/ui/scroll-area"
import { Loader2, AlertCircle } from "lucide-react"
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/alert"

interface InspectionData {
    m_spec: {
        fixed_parameters: Record<string, any>
        materials: Record<string, any>
    }
    m_para: {
        winding_parameters: Record<string, any>
        other_derived: Record<string, any>
    }
    design_parameters: Record<string, any>
    jmag_study: Record<string, any>
}

export default function InspectionPage() {
    const [data, setData] = useState<InspectionData | null>(null)
    const [loading, setLoading] = useState(true)
    const [error, setError] = useState<string | null>(null)

    useEffect(() => {
        const fetchData = async () => {
            try {
                const response = await fetch("http://localhost:8000/api/debug/inspection")
                if (!response.ok) {
                    throw new Error(`Failed to fetch inspection data: ${response.statusText}`)
                }
                const jsonData = await response.json()
                setData(jsonData)
            } catch (err: any) {
                setError(err.message || "An unknown error occurred")
            } finally {
                setLoading(false)
            }
        }

        fetchData()
    }, [])

    if (loading) {
        return (
            <div className="flex h-[80vh] items-center justify-center">
                <Loader2 className="h-8 w-8 animate-spin text-primary" />
                <span className="ml-2 text-lg">Loading inspection data...</span>
            </div>
        )
    }

    if (error) {
        return (
            <div className="p-8">
                <Alert variant="destructive">
                    <AlertCircle className="h-4 w-4" />
                    <AlertTitle>Error</AlertTitle>
                    <AlertDescription>{error}</AlertDescription>
                </Alert>
            </div>
        )
    }

    if (!data) return null

    const renderTable = (obj: Record<string, any>) => (
        <Table>
            <TableHeader>
                <TableRow>
                    <TableHead className="w-1/2">Variable</TableHead>
                    <TableHead>Value</TableHead>
                </TableRow>
            </TableHeader>
            <TableBody>
                {Object.entries(obj).map(([key, value]) => (
                    <TableRow key={key}>
                        <TableCell className="font-mono text-sm">{key}</TableCell>
                        <TableCell>
                            {typeof value === "object" && value !== null ? (
                                <pre className="text-xs bg-muted p-2 rounded max-h-40 overflow-auto">
                                    {JSON.stringify(value, null, 2)}
                                </pre>
                            ) : (
                                <span className="font-medium">
                                    {typeof value === "number" ? value.toFixed(4).replace(/\.?0+$/, "") : String(value)}
                                </span>
                            )}
                        </TableCell>
                    </TableRow>
                ))}
            </TableBody>
        </Table>
    )

    return (
        <div className="flex-1 space-y-4 p-8 pt-6">
            <div className="flex items-center justify-between space-y-2">
                <h2 className="text-3xl font-bold tracking-tight">Machine Inspection</h2>
            </div>
            <Tabs defaultValue="m_spec" className="space-y-4">
                <TabsList>
                    <TabsTrigger value="m_spec">Machine Spec</TabsTrigger>
                    <TabsTrigger value="m_para">Machine Para</TabsTrigger>
                    <TabsTrigger value="design">Design Parameters</TabsTrigger>
                    <TabsTrigger value="jmag">JMAG Study</TabsTrigger>
                </TabsList>

                <TabsContent value="m_spec" className="space-y-4">
                    <div className="grid gap-4 md:grid-cols-2 lg:grid-cols-7">
                        <Card className="col-span-4">
                            <CardHeader>
                                <CardTitle>Fixed Parameters</CardTitle>
                                <CardDescription>Geometric constraints and sizing</CardDescription>
                            </CardHeader>
                            <CardContent>
                                <ScrollArea className="h-[500px] pr-4">
                                    {renderTable(data.m_spec.fixed_parameters)}
                                </ScrollArea>
                            </CardContent>
                        </Card>
                        <Card className="col-span-3">
                            <CardHeader>
                                <CardTitle>Materials</CardTitle>
                                <CardDescription>Steel and magnet properties</CardDescription>
                            </CardHeader>
                            <CardContent>
                                <ScrollArea className="h-[500px] pr-4">
                                    {renderTable(data.m_spec.materials)}
                                </ScrollArea>
                            </CardContent>
                        </Card>
                    </div>
                </TabsContent>

                <TabsContent value="m_para" className="space-y-4">
                    <div className="grid gap-4 md:grid-cols-2 lg:grid-cols-7">
                        <Card className="col-span-4">
                            <CardHeader>
                                <CardTitle>Winding Parameters</CardTitle>
                                <CardDescription>Electromagnetic specifications</CardDescription>
                            </CardHeader>
                            <CardContent>
                                <ScrollArea className="h-[500px] pr-4">
                                    {renderTable(data.m_para.winding_parameters)}
                                </ScrollArea>
                            </CardContent>
                        </Card>
                        <Card className="col-span-3">
                            <CardHeader>
                                <CardTitle>Other Derived</CardTitle>
                                <CardDescription>Calculated physical properties</CardDescription>
                            </CardHeader>
                            <CardContent>
                                <ScrollArea className="h-[500px] pr-4">
                                    {renderTable(data.m_para.other_derived)}
                                </ScrollArea>
                            </CardContent>
                        </Card>
                    </div>
                </TabsContent>

                <TabsContent value="design" className="space-y-4">
                    <Card>
                        <CardHeader>
                            <CardTitle>Design Parameters</CardTitle>
                            <CardDescription>Free variables in the search space</CardDescription>
                        </CardHeader>
                        <CardContent>
                            <ScrollArea className="h-[600px] pr-4">
                                {renderTable(data.design_parameters)}
                            </ScrollArea>
                        </CardContent>
                    </Card>
                </TabsContent>

                <TabsContent value="jmag" className="space-y-4">
                    <Card>
                        <CardHeader>
                            <CardTitle>JMAG Study Configuration</CardTitle>
                            <CardDescription>FEA simulation settings and mesh details</CardDescription>
                        </CardHeader>
                        <CardContent>
                            <ScrollArea className="h-[600px] pr-4">
                                {renderTable(data.jmag_study)}
                            </ScrollArea>
                        </CardContent>
                    </Card>
                </TabsContent>
            </Tabs>
        </div>
    )
}
