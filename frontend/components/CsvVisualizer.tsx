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
    selectedFile?: string; // 外部控制的选中文件（用于记忆功能）
    onFileSelect?: (fileName: string) => void; // 回调函数，通知父组件用户选择了哪个文件
    onCurrentFileChange?: (filePath: string) => void; // 回调函数，通知父组件当前选中的 CSV 文件路径
    csvDataCache?: Map<string, any>; // CSV数据缓存：Map<fileName, parsedData>
    loadingCache?: boolean; // 是否正在加载缓存
}

const BACKEND_URL = process.env.NEXT_PUBLIC_BACKEND_URL || "http://localhost:8000";

export default function CsvVisualizer({ projectName, path2FEACsv, selectedFile: externalSelectedFile, onFileSelect, onCurrentFileChange, csvDataCache, loadingCache }: CsvVisualizerProps) {
    const { theme } = useTheme();
    const [csvFiles, setCsvFiles] = useState<string[]>([]);
    const [internalSelectedFile, setInternalSelectedFile] = useState<string>('');
    const [chartData, setChartData] = useState<any[]>([]);
    const [dataKeys, setDataKeys] = useState<string[]>([]);
    const [loading, setLoading] = useState(false);
    const [error, setError] = useState<string | null>(null);
    
    // 使用外部传入的 selectedFile，如果没有则使用内部状态
    const selectedFile = externalSelectedFile !== undefined ? externalSelectedFile : internalSelectedFile;

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
                if (externalSelectedFile === undefined) {
                    setInternalSelectedFile('');
                }
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
                
                const files = response.data || [];
                console.log('CsvVisualizer: CSV files received:', files);
                console.log('CsvVisualizer: Number of files:', files.length);
                
                setCsvFiles(files);
                if (response.data && response.data.length > 0) {
                    // 如果有外部传入的 selectedFile，且新文件列表中包含它，则使用它（记忆功能）
                    // 否则使用第一个文件
                    let fileToSelect: string;
                    if (externalSelectedFile && files.includes(externalSelectedFile)) {
                        fileToSelect = externalSelectedFile;
                    } else {
                        fileToSelect = response.data[0];
                    }
                    
                    // 更新内部状态（如果外部没有控制）
                    if (externalSelectedFile === undefined) {
                        setInternalSelectedFile(fileToSelect);
                    }
                    
                    // 通知父组件文件选择
                    if (onFileSelect) {
                        onFileSelect(fileToSelect);
                    }
                    
                    // 通知父组件当前文件路径
                    if (onCurrentFileChange && path2FEACsv) {
                        onCurrentFileChange(`${path2FEACsv}/${fileToSelect}`);
                    }
                } else {
                    if (externalSelectedFile === undefined) {
                        setInternalSelectedFile('');
                    }
                    if (onFileSelect) {
                        onFileSelect('');
                    }
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
    }, [path2FEACsv, externalSelectedFile]);

    useEffect(() => {
        const fetchCsvContent = async () => {
            if (!selectedFile || !path2FEACsv) return;

            // 优先从缓存中读取数据
            if (csvDataCache && csvDataCache.has(selectedFile)) {
                console.log('Using cached CSV data for:', selectedFile);
                const cachedData = csvDataCache.get(selectedFile);
                
                if (cachedData && cachedData.length > 0) {
                    // 使用缓存的数据
                    const keys = Object.keys(cachedData[0]).filter(key => key !== 'Time(s)' && !isNaN(parseFloat(cachedData[0][key] as string)));
                    setDataKeys(keys);

                    // Convert values to numbers
                    const formattedData = cachedData.map((row: any) => {
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
                    setError(null);
                    if (onCurrentFileChange && path2FEACsv) {
                        onCurrentFileChange(`${path2FEACsv}/${selectedFile}`);
                    }
                    return;
                }
            }

            // 如果缓存中没有，才从后端请求（这种情况应该很少发生）
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
    }, [selectedFile, path2FEACsv, csvDataCache, onCurrentFileChange]);

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
                        // 更新内部状态（如果外部没有控制）
                        if (externalSelectedFile === undefined) {
                            setInternalSelectedFile(newFile);
                        }
                        // 通知父组件文件选择
                        if (onFileSelect) {
                            onFileSelect(newFile);
                        }
                        // 通知父组件当前文件路径
                        if (onCurrentFileChange && path2FEACsv && newFile) {
                            onCurrentFileChange(`${path2FEACsv}/${newFile}`);
                        } else if (onCurrentFileChange) {
                            onCurrentFileChange('');
                        }
                    }}
                    className="bg-card border border-border text-sm rounded px-3 py-1.5 text-foreground focus:ring-1 focus:ring-primary outline-none min-w-[250px]"
                >
                    {csvFiles.length === 0 ? (
                        <option value="">暂无 CSV 文件</option>
                    ) : (
                        csvFiles.map(file => (
                            <option key={file} value={file}>{file}</option>
                        ))
                    )}
                </select>
                {loading && <Loader2 className="w-4 h-4 animate-spin text-primary" />}
            </div>

            {/* 显示当前状态信息 */}
            <div className="mb-4 text-xs text-muted-foreground space-y-1 bg-muted/30 p-3 rounded-md">
                {path2FEACsv && (
                    <div className="flex items-start space-x-2">
                        <span className="font-medium whitespace-nowrap">CSV 路径:</span>
                        <span className="font-mono break-all">{path2FEACsv}</span>
                    </div>
                )}
                {loading && (
                    <div className="flex items-center space-x-2 text-primary">
                        <Loader2 className="w-3 h-3 animate-spin" />
                        <span>正在加载 CSV 文件列表...</span>
                    </div>
                )}
                {!loading && csvFiles.length === 0 && path2FEACsv && (
                    <div className="flex items-center space-x-2 text-amber-600 dark:text-amber-400">
                        <AlertCircle className="w-3 h-3" />
                        <span>未找到 CSV 文件（路径可能不存在或为空）</span>
                    </div>
                )}
                {csvFiles.length > 0 && (
                    <div className="flex items-center space-x-2 text-green-600 dark:text-green-400">
                        <span>✓ 找到 {csvFiles.length} 个 CSV 文件</span>
                    </div>
                )}
                {selectedFile && (
                    <div className="flex items-start space-x-2">
                        <span className="font-medium whitespace-nowrap">当前文件:</span>
                        <span className="font-mono break-all">{selectedFile}</span>
                    </div>
                )}
                {chartData.length > 0 && (
                    <div className="flex items-center space-x-2 text-green-600 dark:text-green-400">
                        <span>✓ 已加载 {chartData.length} 行数据，{dataKeys.length} 个数据列: {dataKeys.join(', ')}</span>
                    </div>
                )}
            </div>

            {error && (
                <div className="bg-destructive/10 text-destructive text-sm p-3 rounded-md flex items-center mb-4">
                    <AlertCircle className="w-4 h-4 mr-2" />
                    <div className="flex flex-col">
                        <span className="font-medium">错误:</span>
                        <span>{error}</span>
                    </div>
                </div>
            )}

            <div className="flex-1 min-h-[400px] bg-card rounded-lg border border-border p-4 relative" style={{ height: '500px' }}>
                {(loading || loadingCache) ? (
                    <div className="absolute inset-0 flex items-center justify-center bg-background/80 backdrop-blur-sm z-10">
                        <div className="flex flex-col items-center space-y-2">
                            <Loader2 className="w-6 h-6 animate-spin text-primary" />
                            <span className="text-sm text-muted-foreground">
                                {loadingCache ? '正在加载 CSV 文件缓存...' : '加载 CSV 数据中...'}
                            </span>
                        </div>
                    </div>
                ) : chartData.length > 0 ? (
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
                    <div className="absolute inset-0 flex items-center justify-center">
                        <div className="flex flex-col items-center space-y-2 text-muted-foreground">
                            {!path2FEACsv ? (
                                <>
                                    <AlertCircle className="w-8 h-8" />
                                    <span className="text-sm">未提供 CSV 路径</span>
                                </>
                            ) : csvFiles.length === 0 ? (
                                <>
                                    <AlertCircle className="w-8 h-8" />
                                    <span className="text-sm">未找到 CSV 文件</span>
                                    <span className="text-xs">路径: {path2FEACsv}</span>
                                </>
                            ) : !selectedFile ? (
                                <>
                                    <AlertCircle className="w-8 h-8" />
                                    <span className="text-sm">请选择一个 CSV 文件</span>
                                </>
                            ) : (
                                <>
                                    <AlertCircle className="w-8 h-8" />
                                    <span className="text-sm">无法加载 CSV 数据</span>
                                    <span className="text-xs">文件: {selectedFile}</span>
                                </>
                            )}
                        </div>
                    </div>
                )}
            </div>
        </div>
    );
}
