"use client";

import React, { useState, useEffect } from 'react';
import axios from 'axios';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer } from 'recharts';
import { Loader2, Filter, Check } from 'lucide-react';
import * as d3 from 'd3';
import { Popover, PopoverContent, PopoverTrigger } from "@/components/ui/popover";
import { Button } from "@/components/ui/button";
import { cn } from "@/lib/utils";

interface MultiCsvVisualizerProps {
    projectName: string;
}

interface ChartData {
    filename: string;
    data: any[];
    keys: string[];
}

export default function MultiCsvVisualizer({ projectName }: MultiCsvVisualizerProps) {
    const [csvFiles, setCsvFiles] = useState<string[]>([]);
    const [selectedFiles, setSelectedFiles] = useState<string[]>([]);
    const [chartsData, setChartsData] = useState<ChartData[]>([]);
    const [loading, setLoading] = useState(false);
    const [open, setOpen] = useState(false);

    // Fetch file list
    useEffect(() => {
        const fetchCsvList = async () => {
            try {
                const response = await axios.get(`http://localhost:8000/api/results/csv/list/${projectName}`);
                const files = response.data;
                setCsvFiles(files);
                setSelectedFiles(files); // Default select all
            } catch (err) {
                console.error("Failed to fetch CSV list", err);
            }
        };

        if (projectName) {
            fetchCsvList();
        }
    }, [projectName]);

    // Fetch content for selected files
    useEffect(() => {
        const fetchAllContent = async () => {
            if (selectedFiles.length === 0) {
                setChartsData([]);
                return;
            }

            setLoading(true);

            try {
                const promises = selectedFiles.map(async (file) => {
                    try {
                        const response = await axios.get(`http://localhost:8000/api/results/csv/content/${projectName}/${file}`);
                        const csvContent = response.data.content;

                        const lines = csvContent.split('\n');
                        const headerIndex = lines.findIndex((line: string) => line.trim().startsWith('Time(s)'));

                        if (headerIndex === -1) return null;

                        const cleanCsvContent = lines.slice(headerIndex).join('\n');
                        const parsedData = d3.csvParse(cleanCsvContent);

                        if (parsedData.length > 0) {
                            const keys = Object.keys(parsedData[0]).filter(key => key !== 'Time(s)' && !isNaN(parseFloat(parsedData[0][key] as string)));

                            const formattedData = parsedData.map(row => {
                                const newRow: any = { ...row };
                                Object.keys(row).forEach(key => {
                                    const val = parseFloat(row[key] as string);
                                    if (!isNaN(val)) {
                                        newRow[key] = val;
                                    }
                                });
                                return newRow;
                            });

                            return {
                                filename: file,
                                data: formattedData,
                                keys: keys
                            };
                        }
                    } catch (e) {
                        console.error(`Error fetching ${file}`, e);
                    }
                    return null;
                });

                const results = await Promise.all(promises);
                const validResults = results.filter(r => r !== null) as ChartData[];
                setChartsData(validResults);

            } catch (err) {
                console.error("Failed to fetch CSV contents", err);
            } finally {
                setLoading(false);
            }
        };

        // Debounce slightly to avoid thrashing if user selects quickly
        const timeoutId = setTimeout(() => {
            fetchAllContent();
        }, 500);

        return () => clearTimeout(timeoutId);
    }, [selectedFiles, projectName]);

    const toggleFile = (file: string) => {
        setSelectedFiles(prev =>
            prev.includes(file)
                ? prev.filter(f => f !== file)
                : [...prev, file]
        );
    };

    return (
        <div className="flex flex-col h-full space-y-4">
            <div className="flex items-center justify-between bg-card p-4 rounded-lg border border-border">
                <h2 className="text-lg font-semibold">CSV Results Explorer</h2>

                <Popover open={open} onOpenChange={setOpen}>
                    <PopoverTrigger asChild>
                        <Button variant="outline" role="combobox" aria-expanded={open} className="w-[300px] justify-between">
                            <Filter className="mr-2 h-4 w-4" />
                            {selectedFiles.length} files selected
                        </Button>
                    </PopoverTrigger>
                    <PopoverContent className="w-[300px] p-4">
                        <div className="space-y-2">
                            <div className="flex items-center space-x-2 pb-2 border-b border-border">
                                <div
                                    className={cn(
                                        "flex h-4 w-4 items-center justify-center rounded-sm border border-primary cursor-pointer",
                                        selectedFiles.length === csvFiles.length ? "bg-primary text-primary-foreground" : "opacity-50"
                                    )}
                                    onClick={() => {
                                        if (selectedFiles.length === csvFiles.length) setSelectedFiles([]);
                                        else setSelectedFiles([...csvFiles]);
                                    }}
                                >
                                    <Check className="h-3 w-3" />
                                </div>
                                <span className="text-sm font-medium">Select All</span>
                            </div>
                            <div className="max-h-[300px] overflow-y-auto space-y-1">
                                {csvFiles.map((file) => (
                                    <div key={file} className="flex items-center space-x-2 py-1">
                                        <div
                                            className={cn(
                                                "flex h-4 w-4 items-center justify-center rounded-sm border border-primary cursor-pointer",
                                                selectedFiles.includes(file)
                                                    ? "bg-primary text-primary-foreground"
                                                    : "opacity-50"
                                            )}
                                            onClick={() => toggleFile(file)}
                                        >
                                            <Check className="h-3 w-3" />
                                        </div>
                                        <span className="text-sm truncate cursor-pointer" onClick={() => toggleFile(file)} title={file}>
                                            {file}
                                        </span>
                                    </div>
                                ))}
                            </div>
                        </div>
                    </PopoverContent>
                </Popover>
            </div>

            {loading && (
                <div className="flex items-center justify-center py-8">
                    <Loader2 className="h-8 w-8 animate-spin text-primary" />
                    <span className="ml-2 text-muted-foreground">Loading data...</span>
                </div>
            )}

            <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-4 pb-8">
                {chartsData.map((chart) => (
                    <div key={chart.filename} className="bg-card border border-border rounded-lg p-3 flex flex-col h-[250px]">
                        <div className="text-xs font-medium text-muted-foreground mb-2 truncate" title={chart.filename}>
                            {chart.filename}
                        </div>
                        <div className="flex-1 min-h-0">
                            <ResponsiveContainer width="100%" height="100%">
                                <LineChart data={chart.data}>
                                    <CartesianGrid strokeDasharray="3 3" stroke="#333" />
                                    <XAxis
                                        dataKey="Time(s)"
                                        type="number"
                                        domain={['auto', 'auto']}
                                        tick={{ fontSize: 10 }}
                                        tickFormatter={(tick) => tick.toFixed(3)}
                                    />
                                    <YAxis
                                        width={50}
                                        tick={{ fontSize: 10 }}
                                        domain={['auto', 'auto']}
                                        label={{
                                            value: chart.keys[0],
                                            angle: -90,
                                            position: 'insideLeft',
                                            style: { textAnchor: 'middle', fontSize: '10px' }
                                        }}
                                    />
                                    <Tooltip
                                        contentStyle={{ backgroundColor: '#1f2937', border: '1px solid #374151', color: '#f3f4f6', fontSize: '12px' }}
                                        labelFormatter={(label) => `Time: ${Number(label).toFixed(4)}s`}
                                    />
                                    {chart.keys.map((key, index) => (
                                        <Line
                                            key={key}
                                            type="monotone"
                                            dataKey={key}
                                            stroke={`hsl(${index * 137.5}, 70%, 50%)`} // Golden angle for distinct colors
                                            dot={false}
                                            strokeWidth={1.5}
                                        />
                                    ))}
                                </LineChart>
                            </ResponsiveContainer>
                        </div>
                    </div>
                ))}
            </div>
        </div>
    );
}
