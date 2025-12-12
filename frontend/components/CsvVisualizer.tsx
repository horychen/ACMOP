"use client";

import React, { useState, useEffect, useMemo } from 'react';
import axios from 'axios';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer } from 'recharts';
import { Loader2, FileText, AlertCircle } from 'lucide-react';
import * as d3 from 'd3';
import { useTheme } from '@/context/ThemeContext';

interface CsvVisualizerProps {
    projectName?: string;
    path2FEACsv?: string;
    onCurrentFileChange?: (filePath: string) => void; // 回调函数，通知父组件当前选中的 CSV 文件路径
}

const BACKEND_URL = process.env.NEXT_PUBLIC_BACKEND_URL || "http://localhost:8000";

export default function CsvVisualizer({ projectName, path2FEACsv, onCurrentFileChange }: CsvVisualizerProps) {
    const { theme } = useTheme();
    const [csvFiles, setCsvFiles] = useState<string[]>([]);
    const [selectedFile, setSelectedFile] = useState<string>('');
    const [chartData, setChartData] = useState<any[]>([]);
    const [dataKeys, setDataKeys] = useState<string[]>([]);
    const [loading, setLoading] = useState(false);
    const [error, setError] = useState<string | null>(null);

    // Chart colors based on theme
    const chartColors = useMemo(() => {
        const isDark = theme === 'dark';
        return {
            // Light mode: softer, more elegant colors with better contrast
            // Dark mode: maintains current professional look
            gridStroke: isDark ? '#4b5563' : '#e5e7eb', // Slate 600 for dark, Slate 200 for light (softer grid)
            tooltipBg: isDark ? '#1f2937' : '#ffffff', // Slate 800 for dark, pure white for light
            tooltipBorder: isDark ? '#374151' : '#d1d5db', // Slate 700 for dark, Slate 300 for light (subtle border)
            tooltipText: isDark ? '#f3f4f6' : '#111827', // Slate 100 for dark, Gray 900 for light (better contrast)
            axisStroke: isDark ? '#9ca3af' : '#6b7280', // Slate 400 for dark, Gray 500 for light (balanced visibility)
            axisTick: isDark ? '#9ca3af' : '#6b7280', // Same as axis stroke
            lineLightness: isDark ? 50 : 35, // Lighter lines for dark mode, darker for light mode (more visible)
            lineSaturation: isDark ? 70 : 85, // More vibrant, saturated colors in light mode
        };
    }, [theme]);

    useEffect(() => {
        const fetchCsvList = async () => {
            // Only fetch if we have a valid path2FEACsv
            if (!path2FEACsv) {
                setCsvFiles([]);
                setSelectedFile('');
                setError(null);
                if (onCurrentFileChange) {
                    onCurrentFileChange('');
                }
                return;
            }

            console.log('CsvVisualizer: Fetching CSV list from path:', path2FEACsv);
            try {
                setError(null);
                const response = await axios.get(`${BACKEND_URL}/api/results/csv/list-from-path`, {
                    params: { path: path2FEACsv }
                });
                
                console.log('CsvVisualizer: CSV files received:', response.data);
                
                setCsvFiles(response.data || []);
                if (response.data && response.data.length > 0) {
                    const firstFile = response.data[0];
                    setSelectedFile(firstFile);
                    // 通知父组件当前文件路径
                    if (onCurrentFileChange && path2FEACsv) {
                        onCurrentFileChange(`${path2FEACsv}/${firstFile}`);
                    }
                } else {
                    setSelectedFile('');
                    if (onCurrentFileChange) {
                        onCurrentFileChange('');
                    }
                }
            } catch (err: any) {
                console.error("Failed to fetch CSV list", err);
                console.error("Path used:", path2FEACsv);
                console.error("Error details:", err.response?.data);
                const errorMessage = err.response?.data?.detail || err.message || "无法加载 CSV 文件列表";
                setError(errorMessage);
                setCsvFiles([]);
                if (onCurrentFileChange) {
                    onCurrentFileChange('');
                }
            }
        };

        fetchCsvList();
    }, [path2FEACsv]);

    useEffect(() => {
        const fetchCsvContent = async () => {
            if (!selectedFile || !path2FEACsv) return;

            setLoading(true);
            setError(null);
            try {
                const response = await axios.get(`${BACKEND_URL}/api/results/csv/content-from-path`, {
                    params: { path: path2FEACsv, filename: selectedFile }
                });
                
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
                    // Filter out columns that are not suitable for plotting (e.g., non-numeric)
                    // We assume 'Time(s)' is the X-axis.
                    const keys = Object.keys(parsedData[0]).filter(key => key !== 'Time(s)' && !isNaN(parseFloat(parsedData[0][key] as string)));
                    setDataKeys(keys);

                    // Convert values to numbers
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
                    setChartData(formattedData);
                    // 确保通知父组件当前文件路径
                    if (onCurrentFileChange && path2FEACsv && selectedFile) {
                        onCurrentFileChange(`${path2FEACsv}/${selectedFile}`);
                    }
                } else {
                    setChartData([]);
                    setDataKeys([]);
                }

            } catch (err: any) {
                console.error("Failed to fetch CSV content", err);
                console.error("Path used:", path2FEACsv);
                console.error("File selected:", selectedFile);
                console.error("Error details:", err.response?.data);
                const errorMessage = err.response?.data?.detail || err.message || "无法加载 CSV 内容";
                setError(errorMessage);
                setChartData([]);
                setDataKeys([]);
            } finally {
                setLoading(false);
            }
        };

        fetchCsvContent();
    }, [selectedFile, path2FEACsv]);

    return (
        <div className="flex flex-col h-full">
            <div className="mb-4 flex items-center space-x-4">
                <div className="flex items-center space-x-2">
                    <FileText className="w-4 h-4 text-muted-foreground" />
                    <span className="text-sm font-medium text-foreground">Select Result:</span>
                </div>
                <select
                    value={selectedFile}
                    onChange={(e) => {
                        const newFile = e.target.value;
                        setSelectedFile(newFile);
                        // 通知父组件当前文件路径
                        if (onCurrentFileChange && path2FEACsv && newFile) {
                            onCurrentFileChange(`${path2FEACsv}/${newFile}`);
                        } else if (onCurrentFileChange) {
                            onCurrentFileChange('');
                        }
                    }}
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

            <div className="flex-1 min-h-[400px] bg-card rounded-lg border border-border p-4 relative">
                {chartData.length > 0 ? (
                    <ResponsiveContainer width="100%" height="100%">
                        <LineChart data={chartData}>
                            <CartesianGrid strokeDasharray="3 3" stroke={chartColors.gridStroke} />
                            <XAxis
                                dataKey="Time(s)"
                                type="number"
                                domain={['auto', 'auto']}
                                tickFormatter={(tick) => tick.toFixed(4)}
                                label={{ value: 'Time (s)', position: 'insideBottomRight', offset: -5 }}
                                stroke={chartColors.axisStroke}
                                tick={{ fill: chartColors.axisTick }}
                            />
                            <YAxis 
                                stroke={chartColors.axisStroke}
                                tick={{ fill: chartColors.axisTick }}
                            />
                            <Tooltip
                                contentStyle={{ 
                                    backgroundColor: chartColors.tooltipBg, 
                                    border: `1px solid ${chartColors.tooltipBorder}`, 
                                    color: chartColors.tooltipText,
                                    borderRadius: '6px'
                                }}
                            />
                            <Legend />
                            {dataKeys.map((key, index) => (
                                <Line
                                    key={key}
                                    type="monotone"
                                    dataKey={key}
                                    stroke={`hsl(${index * 60}, ${chartColors.lineSaturation}%, ${chartColors.lineLightness}%)`}
                                    dot={false}
                                    strokeWidth={2}
                                />
                            ))}
                        </LineChart>
                    </ResponsiveContainer>
                ) : (
                    <div className="absolute inset-0 flex items-center justify-center text-muted-foreground">
                        {!loading && "No data to display"}
                    </div>
                )}
            </div>
        </div>
    );
}
