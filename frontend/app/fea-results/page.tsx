"use client"

import { useState, useEffect } from "react";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "@/components/ui/table";
import { Loader2, AlertCircle } from "lucide-react";
import axios from "axios";

// Standard port for V2 backend according to previous implementation
const BACKEND_URL = process.env.NEXT_PUBLIC_BACKEND_URL || "http://localhost:8000";

export default function FEAResultsPage() {
    const [data, setData] = useState<any>(null);
    const [selectedInd, setSelectedInd] = useState<string>("");
    const [loading, setLoading] = useState(true);
    const [error, setError] = useState<string | null>(null);

    useEffect(() => {
        const fetchData = async () => {
            try {
                const response = await axios.get(`${BACKEND_URL}/api/fea-results`);
                if (response.data.error) {
                    setError(response.data.error);
                } else {
                    setData(response.data);
                    const keys = Object.keys(response.data);
                    if (keys.length > 0) {
                        setSelectedInd(keys[0]);
                    }
                }
            } catch (err: any) {
                console.error("Failed to fetch FEA results:", err);
                setError(err.message || "Failed to connect to backend");
            } finally {
                setLoading(false);
            }
        };
        fetchData();
    }, []);

    if (loading) {
        return (
            <div className="flex h-[calc(100vh-4rem)] items-center justify-center">
                <Loader2 className="w-8 h-8 animate-spin text-primary" />
                <span className="ml-2">Loading FEA Results...</span>
            </div>
        );
    }

    if (error) {
        return (
            <div className="flex h-[calc(100vh-4rem)] items-center justify-center">
                <Card className="max-w-md w-full border-destructive">
                    <CardHeader>
                        <CardTitle className="flex items-center text-destructive">
                            <AlertCircle className="w-5 h-5 mr-2" />
                            Error Loading Data
                        </CardTitle>
                    </CardHeader>
                    <CardContent>
                        <p className="text-muted-foreground">{error}</p>
                    </CardContent>
                </Card>
            </div>
        );
    }

    if (!data || Object.keys(data).length === 0) {
        return (
            <div className="flex h-[calc(100vh-4rem)] items-center justify-center">
                <Card className="max-w-md w-full">
                    <CardContent className="pt-6 text-center text-muted-foreground">
                        No FEA results found. Run a JMAG evaluation first.
                    </CardContent>
                </Card>
            </div>
        );
    }

    const individualKeys = Object.keys(data);
    const currentData = data[selectedInd];

    return (
        <div className="container mx-auto py-8 flex flex-col gap-6 max-h-[calc(100vh-4rem)] overflow-y-auto">
            <div className="flex items-center justify-between">
                <div>
                    <h1 className="text-3xl font-bold tracking-tight">FEA Results Visualization</h1>
                    <p className="text-muted-foreground mt-1">Review validation results exported from JMAG JSON.</p>
                </div>
                {individualKeys.length > 1 && (
                    <div className="flex items-center gap-2">
                        <span className="text-sm font-medium">Select Individual:</span>
                        <Select value={selectedInd} onValueChange={setSelectedInd}>
                            <SelectTrigger className="w-[200px]">
                                <SelectValue placeholder="Select..." />
                            </SelectTrigger>
                            <SelectContent>
                                {individualKeys.map(key => (
                                    <SelectItem key={key} value={key}>{key}</SelectItem>
                                ))}
                            </SelectContent>
                        </Select>
                    </div>
                )}
            </div>

            {currentData && (
                <>
                    <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-4">
                        <MetricCard title="Average Torque" value={currentData.torque_average} unit="Nm" />
                        <MetricCard title="Torque Ripple" value={currentData.normalized_torque_ripple ? currentData.normalized_torque_ripple * 100 : null} unit="%" />
                        <MetricCard title="Total Loss" value={currentData.rated_total_loss} unit="W" />
                        <MetricCard title="Material Cost" value={currentData.Cost} unit="$" />
                    </div>

                    <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
                        <Card className="col-span-1 md:col-span-2 shadow-sm">
                            <CardHeader>
                                <CardTitle>Detailed Specifications</CardTitle>
                                <CardDescription>Comprehensive metrics output from the JMAG runner</CardDescription>
                            </CardHeader>
                            <CardContent className="max-h-[500px] overflow-y-auto w-full pr-4">
                                <Table>
                                    <TableHeader>
                                        <TableRow>
                                            <TableHead>Metric</TableHead>
                                            <TableHead>Value</TableHead>
                                        </TableRow>
                                    </TableHeader>
                                    <TableBody>
                                        {Object.entries(currentData).map(([key, val]) => {
                                            if (key === 'x_denorm_dict' || typeof val === 'object') return null;
                                            return (
                                                <TableRow key={key}>
                                                    <TableCell className="font-mono text-xs text-muted-foreground">{key}</TableCell>
                                                    <TableCell className="font-medium">
                                                        {typeof val === 'number' ? parseFloat(val.toFixed(4)) : String(val)}
                                                    </TableCell>
                                                </TableRow>
                                            );
                                        })}
                                    </TableBody>
                                </Table>
                            </CardContent>
                        </Card>

                        <Card className="col-span-1 shadow-sm">
                            <CardHeader>
                                <CardTitle>Design Parameters</CardTitle>
                                <CardDescription>Target Free Parameters</CardDescription>
                            </CardHeader>
                            <CardContent className="max-h-[500px] overflow-y-auto">
                                <div className="space-y-4">
                                    {currentData.x_denorm_dict && Object.keys(currentData.x_denorm_dict).length > 0 ? (
                                        Object.entries(currentData.x_denorm_dict).map(([k, v]: [string, any]) => (
                                            <div key={k} className="flex justify-between items-center border-b pb-2">
                                                <span className="text-sm text-muted-foreground font-mono truncate mr-2" title={k}>{k}</span>
                                                <span className="text-sm font-semibold">{typeof v === 'number' ? v.toFixed(3) : String(v)}</span>
                                            </div>
                                        ))
                                    ) : (
                                        <div className="text-sm text-muted-foreground">No parametric config available.</div>
                                    )}
                                </div>
                            </CardContent>
                        </Card>
                    </div>
                </>
            )}
        </div>
    );
}

function MetricCard({ title, value, unit }: { title: string, value: any, unit: string }) {
    const formattedValue = typeof value === 'number' ? value.toFixed(2) : value || 'N/A';
    return (
        <Card>
            <CardHeader className="pb-2">
                <CardTitle className="text-sm font-medium text-muted-foreground uppercase tracking-wider">{title}</CardTitle>
            </CardHeader>
            <CardContent>
                <div className="text-2xl font-bold">
                    {formattedValue} <span className="text-sm font-normal text-muted-foreground">{unit}</span>
                </div>
            </CardContent>
        </Card>
    );
}
