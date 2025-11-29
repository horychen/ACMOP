"use client";

import React, { useState, useEffect } from 'react';
import axios from 'axios';
import { EfficiencyChart } from './Charts';
import { Loader2, FileText, AlertCircle } from 'lucide-react';
import * as d3 from 'd3';

interface CsvChartVisualizerProps {
    path2FEACsv?: string;
    projectName?: string;
}

const BACKEND_URL = process.env.NEXT_PUBLIC_BACKEND_URL || "http://localhost:8000";

interface ChartDataPoint {
    time: number;
    value: number;
}

export default function CsvChartVisualizer({ path2FEACsv, projectName }: CsvChartVisualizerProps) {
    const [csvFiles, setCsvFiles] = useState<string[]>([]);
    const [selectedFile, setSelectedFile] = useState<string>('');
    const [chartData, setChartData] = useState<ChartDataPoint[]>([]);
    const [loading, setLoading] = useState(false);
    const [error, setError] = useState<string | null>(null);

    useEffect(() => {
        const fetchCsvList = async () => {
            try {
                let response;
                if (path2FEACsv) {
                    response = await axios.get(`${BACKEND_URL}/api/results/csv/list-from-path`, {
                        params: { path: path2FEACsv }
                    });
                } else if (projectName) {
                    response = await axios.get(`${BACKEND_URL}/api/results/csv/list/${projectName}`);
                } else {
                    setError("Either projectName or path2FEACsv must be provided");
                    return;
                }
                
                setCsvFiles(response.data);
                if (response.data.length > 0) {
                    setSelectedFile(response.data[0]);
                }
            } catch (err) {
                console.error("Failed to fetch CSV list", err);
                setError("Failed to load CSV files.");
            }
        };

        if (projectName || path2FEACsv) {
            fetchCsvList();
        }
    }, [projectName, path2FEACsv]);

    useEffect(() => {
        const fetchCsvContent = async () => {
            if (!selectedFile) return;

            setLoading(true);
            setError(null);
            try {
                let response;
                if (path2FEACsv) {
                    response = await axios.get(`${BACKEND_URL}/api/results/csv/content-from-path`, {
                        params: { path: path2FEACsv, filename: selectedFile }
                    });
                } else if (projectName) {
                    response = await axios.get(`${BACKEND_URL}/api/results/csv/content/${projectName}/${selectedFile}`);
                } else {
                    setError("Either projectName or path2FEACsv must be provided");
                    return;
                }
                
                const csvContent = response.data.content;

                // Find the start of the actual data (header line starts with "Time(s)")
                const lines = csvContent.split('\n');
                const headerIndex = lines.findIndex((line: string) => line.trim().startsWith('Time(s)'));

                if (headerIndex === -1) {
                    throw new Error("Could not find data header 'Time(s)' in CSV file.");
                }

                const cleanCsvContent = lines.slice(headerIndex).join('\n');

                // Parse CSV using d3
                const parsedData = d3.csvParse(cleanCsvContent);

                if (parsedData.length > 0) {
                    // Try to find speed and efficiency columns
                    // Common column names: Speed, RPM, speed, Efficiency, efficiency, etc.
                    const speedKey = Object.keys(parsedData[0]).find(
                        key => key.toLowerCase().includes('speed') || 
                               key.toLowerCase().includes('rpm') ||
                               key === 'Speed' || key === 'speed'
                    );
                    
                    const efficiencyKey = Object.keys(parsedData[0]).find(
                        key => key.toLowerCase().includes('efficiency') ||
                               key.toLowerCase().includes('eta') ||
                               key === 'Efficiency' || key === 'efficiency'
                    );

                    // Use Time(s) as x-axis and first numeric column (excluding Time(s)) as y-axis
                    const timeKey = 'Time(s)';
                    const numericKeys = Object.keys(parsedData[0]).filter(
                        key => key !== timeKey && !isNaN(parseFloat(parsedData[0][key] as string))
                    );
                    
                    if (numericKeys.length > 0) {
                        // Use first numeric column as y-axis value
                        const formattedData: ChartDataPoint[] = parsedData
                            .map(row => {
                                const time = parseFloat(row[timeKey] as string);
                                const value = parseFloat(row[numericKeys[0]] as string);
                                if (!isNaN(time) && !isNaN(value)) {
                                    return { time, value };
                                }
                                return null;
                            })
                            .filter((item): item is ChartDataPoint => item !== null);
                        
                        setChartData(formattedData);
                    } else {
                        setChartData([]);
                    }
                } else {
                    setChartData([]);
                }

            } catch (err: any) {
                console.error("Failed to fetch CSV content", err);
                setError(err.message || "Failed to load CSV content.");
            } finally {
                setLoading(false);
            }
        };

        fetchCsvContent();
    }, [selectedFile, projectName, path2FEACsv]);

    return (
        <div className="flex flex-col h-full">
            <div className="mb-4 flex items-center space-x-4">
                <div className="flex items-center space-x-2">
                    <FileText className="w-4 h-4 text-muted-foreground" />
                    <span className="text-sm font-medium text-foreground">选择结果文件:</span>
                </div>
                <select
                    value={selectedFile}
                    onChange={(e) => setSelectedFile(e.target.value)}
                    className="bg-card border border-border text-sm rounded px-3 py-1.5 text-foreground focus:ring-1 focus:ring-primary outline-none min-w-[250px]"
                >
                    {csvFiles.map(file => (
                        <option key={file} value={file}>{file}</option>
                    ))}
                </select>
                {loading && <Loader2 className="w-4 h-4 animate-spin text-primary" />}
            </div>

            {error && (
                <div className="bg-destructive/10 text-destructive text-sm p-3 rounded-md flex items-center mb-4">
                    <AlertCircle className="w-4 h-4 mr-2" />
                    {error}
                </div>
            )}

            <div className="flex-1 min-h-[300px]">
                {chartData.length > 0 ? (
                    <EfficiencyChart data={chartData} title={selectedFile || undefined} />
                ) : (
                    <div className="flex items-center justify-center h-full text-muted-foreground">
                        {!loading && "无数据可显示"}
                    </div>
                )}
            </div>
        </div>
    );
}

