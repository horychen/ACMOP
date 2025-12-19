"use client"

import { useState, useEffect, useMemo } from "react";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Button } from "@/components/ui/button";
import { Alert, AlertDescription } from "@/components/ui/alert";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { Loader2, AlertCircle, FolderOpen } from "lucide-react";
import { ScatterChart, Scatter, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, ReferenceLine, Cell, Legend } from "recharts";
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "@/components/ui/table";
import { ScrollArea } from "@/components/ui/scroll-area";
import { ParetoFront2p5D } from "@/components/ParetoFront2p5D";
import { useTheme } from "@/context/ThemeContext";
import axios from "axios";

interface SensitivityData {
  [key: string]: {
    project_name: string;
    f1?: number;
    f2?: number;
    f3?: number;
    [key: string]: any; // 其他性能参数
  };
}

interface ParsedIndividual {
  key: string;
  parameter: string;
  percentage: number;
  data: any;
  generation?: number;
  individual_index?: number;
}

interface SensitivityTableRow {
  parameter: string;
  percentage: number;
  f1?: number;
  f2?: number;
  f3?: number;
  [key: string]: any;
}

export function SensitivityAnalysisViewer() {
  const { theme } = useTheme();
  const isDark = theme === 'dark';
  const [folderPath, setFolderPath] = useState<string>("");
  const [swarmData, setSwarmData] = useState<SensitivityData | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [baselineValues, setBaselineValues] = useState<Record<string, number>>({});
  
  // 获取所有可用的性能指标（提前定义，用于设置默认值）
  const availableMetrics = useMemo(() => {
    if (!swarmData || typeof swarmData !== "object") return [];
    
    const metrics = new Set<string>();
    try {
      Object.values(swarmData).forEach(value => {
        if (value && typeof value === "object" && !Array.isArray(value)) {
          Object.keys(value).forEach(key => {
            if (key !== "project_name" && key !== "x_denorm_dict" && 
                key !== "individual_name" && key !== "number_current_generation" &&
                key !== "individual_index" && typeof value[key] === "number") {
              metrics.add(key);
            }
          });
        }
      });
    } catch (error) {
      console.error("Error processing metrics:", error);
    }
    
    return Array.from(metrics).sort();
  }, [swarmData]);


  // 解析 project_name，提取参数名和百分比
  const parseProjectName = (projectName: string): { parameter: string; percentage: number } | null => {
    // 格式: param-{parameter_name}-pct-{percentage_value}-ind{index}
    // 或者: {study_name}-param-{parameter_name}-pct-{percentage_value}
    const match = projectName.match(/param-([^-]+)-pct-([\d.-]+)/);
    if (match) {
      return {
        parameter: match[1],
        percentage: parseFloat(match[2])
      };
    }
    return null;
  };

  // 解析数据
  const parsedData = useMemo(() => {
    if (!swarmData) {
      console.log("parsedData: swarmData 为空");
      return [];
    }

    const parsed: ParsedIndividual[] = [];
    Object.entries(swarmData).forEach(([key, value]: [string, any]) => {
      if (!value || typeof value !== "object") {
        console.warn("无效的数据项:", key, value);
        return;
      }
      
      const projectName = value.project_name;
      if (!projectName) {
        console.warn("缺少 project_name:", key);
        return;
      }
      
      const parsedInfo = parseProjectName(projectName);
      if (parsedInfo) {
        parsed.push({
          key,
          parameter: parsedInfo.parameter,
          percentage: parsedInfo.percentage,
          data: value,
          generation: value.number_current_generation,
          individual_index: value.individual_index
        });
      } else {
        console.warn("无法解析 project_name:", projectName);
      }
    });

    console.log("解析后的数据:", parsed.length, "个个体");
    console.log("解析后的参数:", [...new Set(parsed.map(p => p.parameter))]);
    return parsed;
  }, [swarmData]);

  // 按参数分组
  const groupedByParameter = useMemo(() => {
    const grouped: Record<string, ParsedIndividual[]> = {};
    parsedData.forEach(item => {
      if (!grouped[item.parameter]) {
        grouped[item.parameter] = [];
      }
      grouped[item.parameter].push(item);
    });

    // 对每个参数的数据按百分比排序
    Object.keys(grouped).forEach(param => {
      grouped[param].sort((a, b) => a.percentage - b.percentage);
    });

    return grouped;
  }, [parsedData]);

  // 生成表格数据
  const tableData = useMemo(() => {
    const rows: SensitivityTableRow[] = [];
    parsedData.forEach(item => {
      rows.push({
        parameter: item.parameter,
        percentage: item.percentage,
        f1: item.data.f1,
        f2: item.data.f2,
        f3: item.data.f3,
        ...item.data // 包含所有其他性能参数
      });
    });

    // 按参数名和百分比排序
    rows.sort((a, b) => {
      if (a.parameter !== b.parameter) {
        return a.parameter.localeCompare(b.parameter);
      }
      return a.percentage - b.percentage;
    });

    return rows;
  }, [parsedData]);


  // 加载 SwarmData.json
  const loadSwarmData = async () => {
    if (!folderPath.trim()) {
      setError("请输入文件夹路径");
      return;
    }

    setLoading(true);
    setError(null);

    try {
      const response = await axios.get("/api/acmopv2/sensitivity-analysis-data", {
        params: {
          folderPath: folderPath.trim()
        }
      });

      const data = response.data;
      console.log("加载的 SwarmData:", data);
      console.log("数据项数量:", Object.keys(data).length);
      setSwarmData(data);
    } catch (err: any) {
      console.error("加载数据失败:", err);
      setError(err.response?.data?.error || err.message || "加载数据失败");
      setSwarmData(null);
    } finally {
      setLoading(false);
    }
  };

  // 为每个性能指标准备图表数据
  const chartsDataByMetric = useMemo(() => {
    if (!swarmData || !groupedByParameter || availableMetrics.length === 0) {
      return {};
    }

    const chartsData: Record<string, {
      data: Array<{
        x: number;
        y: number;
        parameter: string;
        percentage: number;
        individual_index: number;
        color: string;
      }>;
      paramRanges: Array<{ param: string; minIndex: number; maxIndex: number; indices: number[] }>;
      minIndex: number;
      maxIndex: number;
      yMin: number;
      yMax: number;
      yDomainMin: number;
      yDomainMax: number;
      yTicks: number[];
      ticks: number[];
    }> = {};

    const parameters = Object.keys(groupedByParameter).sort();
    if (parameters.length === 0) return {};

    // 为每个性能指标生成数据
    availableMetrics.forEach(metric => {
      const allDataPoints: Array<{
        x: number;
        y: number;
        parameter: string;
        percentage: number;
        individual_index: number;
        color: string;
      }> = [];

      const paramRanges: Array<{ param: string; minIndex: number; maxIndex: number; indices: number[] }> = [];

      // 为每个参数生成不同颜色
      const paramColors: Record<string, string> = {};
      const colorPalette = [
        "#0ea5e9", // sky blue
        "#10b981", // green
        "#f59e0b", // amber
        "#8b5cf6", // purple
        "#ef4444", // red
        "#06b6d4", // cyan
        "#f97316", // orange
        "#ec4899", // pink
      ];
      parameters.forEach((param, idx) => {
        paramColors[param] = colorPalette[idx % colorPalette.length];
      });

      parameters.forEach((param) => {
        const paramData = groupedByParameter[param];
        if (!paramData || paramData.length === 0) return;

        const indices = paramData
          .map(item => {
            const idx = item.individual_index ?? item.data?.individual_index ?? -1;
            if (typeof idx === 'string') {
              const num = parseFloat(idx);
              return isNaN(num) ? -1 : num;
            }
            return typeof idx === 'number' ? idx : -1;
          })
          .filter(idx => idx >= 0)
          .sort((a, b) => a - b);

        if (indices.length === 0) return;

        const minIndex = Math.min(...indices);
        const maxIndex = Math.max(...indices);
        paramRanges.push({ param, minIndex, maxIndex, indices });

        paramData.forEach(item => {
          const value = item.data?.[metric];
          let individualIndex = item.individual_index ?? item.data?.individual_index ?? -1;
          
          if (typeof individualIndex === 'string') {
            const cleaned = individualIndex.replace(/[%\s]/g, '');
            const num = parseFloat(cleaned);
            individualIndex = isNaN(num) ? -1 : num;
          } else if (typeof individualIndex !== 'number') {
            individualIndex = -1;
          }
          
          if (value !== undefined && value !== null && !isNaN(value) && individualIndex >= 0 && typeof individualIndex === 'number') {
            allDataPoints.push({
              x: individualIndex,
              y: typeof value === 'number' ? value : parseFloat(String(value)) || 0,
              parameter: param,
              percentage: item.percentage,
              individual_index: individualIndex,
              color: paramColors[param] || "#0ea5e9"
            });
          }
        });
      });

      if (allDataPoints.length === 0) return;

      allDataPoints.sort((a, b) => a.x - b.x);

      const minIndex = Math.min(...allDataPoints.map(d => d.x));
      const maxIndex = Math.max(...allDataPoints.map(d => d.x));
      const tickStep = Math.max(1, Math.floor((maxIndex - minIndex) / 10));
      const ticks: number[] = [];
      for (let i = minIndex; i <= maxIndex; i += tickStep) {
        ticks.push(i);
      }
      if (ticks[ticks.length - 1] !== maxIndex) {
        ticks.push(maxIndex);
      }

      const yValues = allDataPoints.map(d => d.y);
      const yMin = Math.min(...yValues);
      const yMax = Math.max(...yValues);
      const yRange = yMax - yMin;
      const yPadding = yRange * 0.1 || Math.abs(yMin) * 0.1 || 1;
      const yDomainMin = yMin - yPadding;
      const yDomainMax = yMax + yPadding;

      const yTickCount = 8;
      const yTicks: number[] = [];
      for (let i = 0; i < yTickCount; i++) {
        const tickValue = yDomainMin + (yDomainMax - yDomainMin) * (i / (yTickCount - 1));
        yTicks.push(tickValue);
      }

      chartsData[metric] = {
        data: allDataPoints,
        paramRanges,
        minIndex,
        maxIndex,
        yMin,
        yMax,
        yDomainMin,
        yDomainMax,
        yTicks,
        ticks
      };
    });

    return chartsData;
  }, [swarmData, groupedByParameter, availableMetrics]);

  // 计算每个性能指标的基准值
  useEffect(() => {
    if (parsedData.length > 0 && availableMetrics.length > 0) {
      const baselines: Record<string, number> = {};
      availableMetrics.forEach(metric => {
        const baseline = parsedData.reduce((closest, current) => {
          const currentDist = Math.abs(current.percentage);
          const closestDist = Math.abs(closest.percentage);
          return currentDist < closestDist ? current : closest;
        });
        const value = baseline.data[metric];
        if (value !== undefined) {
          baselines[metric] = value;
        }
      });
      setBaselineValues(baselines);
    }
  }, [parsedData, availableMetrics]);


  // 准备 Pareto Front 数据
  const paretoData = useMemo(() => {
    if (!parsedData.length) return null;

    const individuals = parsedData.map(item => ({
      key: item.key,
      generation: item.generation ?? 0,
      individual_index: item.individual_index ?? 0,
      objectives: {
        f1: item.data.f1,
        f2: item.data.f2,
        f3: item.data.f3
      },
      f1: item.data.f1,
      f2: item.data.f2,
      f3: item.data.f3,
      parameter: item.parameter,
      percentage: item.percentage
    }));

    return {
      allIndividuals: individuals,
      paretoFront: individuals, // 敏感性分析的所有个体都显示
      objectives: ["f1", "f2", "f3"]
    };
  }, [parsedData]);

  return (
    <div className="space-y-6">
      <Card>
        <CardHeader>
          <CardTitle>敏感性分析结果</CardTitle>
          <CardDescription>
            从指定目录读取 SwarmData.json 文件并分析敏感性分析结果
          </CardDescription>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="flex gap-4">
            <div className="flex-1">
              <Label htmlFor="folderPath">文件夹路径（相对于 backend 目录）</Label>
              <Input
                id="folderPath"
                value={folderPath}
                onChange={(e) => setFolderPath(e.target.value)}
                placeholder="例如: _default/TIA_prototype_sensitivity_analysis"
                onKeyDown={(e) => {
                  if (e.key === "Enter") {
                    loadSwarmData();
                  }
                }}
              />
            </div>
            <div className="flex items-end">
              <Button onClick={loadSwarmData} disabled={loading || !folderPath.trim()}>
                {loading ? (
                  <>
                    <Loader2 className="h-4 w-4 mr-2 animate-spin" />
                    加载中...
                  </>
                ) : (
                  <>
                    <FolderOpen className="h-4 w-4 mr-2" />
                    加载数据
                  </>
                )}
              </Button>
            </div>
          </div>

          {error && (
            <Alert variant="destructive">
              <AlertCircle className="h-4 w-4" />
              <AlertDescription>{error}</AlertDescription>
            </Alert>
          )}

          {swarmData && (
            <div className="text-sm text-muted-foreground">
              已加载 {Object.keys(swarmData).length} 个个体
            </div>
          )}
        </CardContent>
      </Card>

      {swarmData && parsedData.length > 0 && (
        <>
          {/* Pareto Front 图 */}
          <Card>
            <CardHeader>
              <CardTitle>Pareto Front</CardTitle>
              <CardDescription>敏感性分析结果的 Pareto 前沿可视化</CardDescription>
            </CardHeader>
            <CardContent>
              {paretoData && (
                <ParetoFront2p5D
                  individuals={paretoData.allIndividuals}
                  objectives={paretoData.objectives}
                  comp={[0, 1]}
                  upToRankNo={1}
                />
              )}
            </CardContent>
          </Card>

          {/* 数据表格 */}
          <Card>
            <CardHeader>
              <CardTitle>参数-百分比-性能参数表格</CardTitle>
              <CardDescription>所有个体的详细数据</CardDescription>
            </CardHeader>
            <CardContent>
              <ScrollArea className="h-[400px]">
                <Table>
                  <TableHeader>
                    <TableRow>
                      <TableHead>参数</TableHead>
                      <TableHead>百分比变化</TableHead>
                      <TableHead>f1</TableHead>
                      <TableHead>f2</TableHead>
                      <TableHead>f3</TableHead>
                      {availableMetrics.filter(m => !["f1", "f2", "f3"].includes(m)).slice(0, 10).map(metric => (
                        <TableHead key={metric}>{metric}</TableHead>
                      ))}
                    </TableRow>
                  </TableHeader>
                  <TableBody>
                    {tableData.map((row, idx) => (
                      <TableRow key={idx}>
                        <TableCell className="font-medium">{row.parameter}</TableCell>
                        <TableCell>{(row.percentage * 100).toFixed(1)}%</TableCell>
                        <TableCell>{row.f1?.toFixed(4) ?? "N/A"}</TableCell>
                        <TableCell>{row.f2?.toFixed(4) ?? "N/A"}</TableCell>
                        <TableCell>{row.f3?.toFixed(4) ?? "N/A"}</TableCell>
                        {availableMetrics.filter(m => !["f1", "f2", "f3"].includes(m)).slice(0, 10).map(metric => (
                          <TableCell key={metric}>
                            {typeof row[metric] === "number" ? row[metric].toFixed(4) : "N/A"}
                          </TableCell>
                        ))}
                      </TableRow>
                    ))}
                  </TableBody>
                </Table>
              </ScrollArea>
            </CardContent>
          </Card>

          {/* 敏感性分析图表 - 每个性能指标一张图，显示所有参数 */}
          {Object.keys(chartsDataByMetric).length > 0 && (
            <div className="space-y-8">
              {availableMetrics.map((metric) => {
                const chartData = chartsDataByMetric[metric];
                if (!chartData || chartData.data.length === 0) return null;

                // 按参数分组数据，用于图例和颜色
                const dataByParam: Record<string, typeof chartData.data> = {};
                chartData.data.forEach(d => {
                  if (!dataByParam[d.parameter]) {
                    dataByParam[d.parameter] = [];
                  }
                  dataByParam[d.parameter].push(d);
                });

                // 准备表格数据（所有参数的数据）
                const allTableData = chartData.data
                  .sort((a, b) => {
                    if (a.parameter !== b.parameter) {
                      return a.parameter.localeCompare(b.parameter);
                    }
                    return a.x - b.x;
                  })
                  .map(d => ({
                    parameter: d.parameter,
                    individual_index: d.x,
                    percentage: d.percentage,
                    value: d.y
                  }));

                return (
                  <Card key={metric}>
                    <CardHeader>
                      <CardTitle>{metric}</CardTitle>
                      <CardDescription>所有参数变化对 {metric} 的影响对比</CardDescription>
                    </CardHeader>
                    <CardContent>
                      <div className="grid grid-cols-[1fr_auto] gap-4">
                        {/* 图表 */}
                        <div className={`rounded-lg p-4 border ${
                          isDark 
                            ? 'bg-slate-800 border-slate-700' 
                            : 'bg-white border-slate-200'
                        }`} style={{ width: '100%', minWidth: 400, position: 'relative' }}>
                          <div style={{ width: '100%', height: 400, position: 'relative', minWidth: 400 }}>
                            <ResponsiveContainer width="100%" height={400}>
                              <ScatterChart 
                                margin={{ top: 20, right: 20, bottom: 60, left: 80 }}
                                data={chartData.data}
                              >
                                <CartesianGrid 
                                  strokeDasharray="3 3" 
                                  stroke={isDark ? "#334155" : "#e5e7eb"} 
                                />
                                <XAxis 
                                  type="number" 
                                  dataKey="x" 
                                  domain={[chartData.minIndex - 0.5, chartData.maxIndex + 0.5]}
                                  ticks={chartData.ticks}
                                  tickFormatter={(value) => {
                                    const num = typeof value === 'number' ? value : parseFloat(String(value)) || 0;
                                    return num.toString();
                                  }}
                                  stroke={isDark ? "#94a3b8" : "#6b7280"}
                                  tick={{ fontSize: 12, fill: isDark ? "#94a3b8" : "#6b7280" }}
                                  label={{ 
                                    value: "个体编号", 
                                    position: 'insideBottomRight', 
                                    offset: -5, 
                                    fill: isDark ? "#94a3b8" : "#6b7280" 
                                  }}
                                  allowDataOverflow={false}
                                />
                                <YAxis 
                                  type="number"
                                  domain={[chartData.yDomainMin, chartData.yDomainMax]}
                                  ticks={chartData.yTicks}
                                  tickFormatter={(value) => {
                                    const num = typeof value === 'number' ? value : parseFloat(String(value)) || 0;
                                    if (Math.abs(num) >= 1000) {
                                      return num.toFixed(0);
                                    } else if (Math.abs(num) >= 1) {
                                      return num.toFixed(2);
                                    } else {
                                      return num.toFixed(4);
                                    }
                                  }}
                                  stroke={isDark ? "#94a3b8" : "#6b7280"}
                                  tick={{ fontSize: 12, fill: isDark ? "#94a3b8" : "#6b7280" }}
                                  label={{ 
                                    value: metric, 
                                    angle: -90, 
                                    position: "insideLeft",
                                    fill: isDark ? "#94a3b8" : "#6b7280"
                                  }}
                                  width={80}
                                  allowDataOverflow={false}
                                />
                                <Tooltip
                                  contentStyle={{ 
                                    backgroundColor: isDark ? '#1e293b' : '#ffffff', 
                                    borderColor: isDark ? '#475569' : '#d1d5db', 
                                    color: isDark ? '#f1f5f9' : '#111827',
                                    borderRadius: '6px'
                                  }}
                                  content={({ active, payload }) => {
                                    if (active && payload && payload.length) {
                                      const data = payload[0].payload;
                                      return (
                                        <div className={`rounded-lg border p-3 ${
                                          isDark 
                                            ? 'bg-slate-800 border-slate-700 text-slate-200' 
                                            : 'bg-white border-slate-200 text-slate-800'
                                        }`}>
                                          <p className="font-semibold">参数: {data.parameter}</p>
                                          <p className="text-sm">个体编号: {data.x}</p>
                                          <p className="text-sm">百分比: {(data.percentage * 100).toFixed(1)}%</p>
                                          <p className="text-sm" style={{ color: isDark ? '#38bdf8' : '#0284c7' }}>
                                            {metric}: {data.y.toFixed(4)}
                                          </p>
                                        </div>
                                      );
                                    }
                                    return null;
                                  }}
                                />
                                <Legend 
                                  wrapperStyle={{ paddingTop: '20px' }}
                                  formatter={(value: string, entry: any) => (
                                    <span style={{ 
                                      color: isDark ? '#cbd5e1' : '#475569',
                                      fontSize: '12px'
                                    }}>
                                      {value}
                                    </span>
                                  )}
                                />
                                {/* 为每个参数创建一个Scatter，用不同颜色区分 */}
                                {Object.entries(dataByParam).map(([param, paramData]) => (
                                  <Scatter 
                                    key={param}
                                    data={paramData}
                                    dataKey="y"
                                    name={param}
                                    fill={paramData[0]?.color || "#0ea5e9"}
                                    shape="circle"
                                    isAnimationActive={false}
                                  >
                                    {paramData.map((entry, index) => (
                                      <Cell 
                                        key={`cell-${param}-${entry.individual_index}-${index}`} 
                                        fill={entry.color} 
                                        opacity={0.8}
                                        stroke={entry.color}
                                        strokeWidth={2}
                                      />
                                    ))}
                                  </Scatter>
                                ))}
                                {baselineValues[metric] !== undefined && (
                                  <ReferenceLine 
                                    y={baselineValues[metric]} 
                                    stroke="#000" 
                                    strokeWidth={2} 
                                    strokeDasharray="5 5"
                                    label={{ value: "基准值", position: "right" }}
                                  />
                                )}
                              </ScatterChart>
                            </ResponsiveContainer>
                          </div>
                        </div>

                        {/* 数据表格 */}
                        <div className={`w-[300px] rounded-lg p-4 border ${
                          isDark 
                            ? 'bg-slate-800 border-slate-700' 
                            : 'bg-white border-slate-200'
                        }`}>
                          <ScrollArea className="h-[400px]">
                            <Table>
                              <TableHeader>
                                <TableRow className={isDark ? 'hover:bg-slate-700' : 'hover:bg-slate-50'}>
                                  <TableHead className={`w-[100px] ${isDark ? 'text-slate-300' : 'text-slate-700'}`}>参数</TableHead>
                                  <TableHead className={`w-[60px] ${isDark ? 'text-slate-300' : 'text-slate-700'}`}>编号</TableHead>
                                  <TableHead className={`w-[80px] ${isDark ? 'text-slate-300' : 'text-slate-700'}`}>百分比</TableHead>
                                  <TableHead className={isDark ? 'text-slate-300' : 'text-slate-700'}>值</TableHead>
                                </TableRow>
                              </TableHeader>
                              <TableBody>
                                {allTableData.map((row, idx) => (
                                  <TableRow 
                                    key={`table-${metric}-${row.parameter}-${row.individual_index}-${idx}`}
                                    className={isDark ? 'hover:bg-slate-700' : 'hover:bg-slate-50'}
                                  >
                                    <TableCell className={`font-medium ${isDark ? 'text-slate-200' : 'text-slate-800'}`}>
                                      {row.parameter}
                                    </TableCell>
                                    <TableCell className={`font-medium ${isDark ? 'text-slate-200' : 'text-slate-800'}`}>
                                      {row.individual_index}
                                    </TableCell>
                                    <TableCell className={isDark ? 'text-slate-300' : 'text-slate-600'}>
                                      {(row.percentage * 100).toFixed(1)}%
                                    </TableCell>
                                    <TableCell className={isDark ? 'text-slate-200' : 'text-slate-800'}>
                                      {row.value.toFixed(4)}
                                    </TableCell>
                                  </TableRow>
                                ))}
                              </TableBody>
                            </Table>
                          </ScrollArea>
                        </div>
                      </div>
                    </CardContent>
                  </Card>
                );
              })}
            </div>
          )}
        </>
      )}
    </div>
  );
}
