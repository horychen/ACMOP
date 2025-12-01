"use client"

import { useState, useEffect } from "react";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { ScatterChart, Scatter, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer } from "recharts";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "@/components/ui/table";
import CsvVisualizer from "@/components/CsvVisualizer";
import axios from "axios";

interface OptimizationResultsViewerProps {
  data: any;
}

export function OptimizationResultsViewer({ data }: OptimizationResultsViewerProps) {
  const [paretoData, setParetoData] = useState<any>(null);
  const [selectedIndividual, setSelectedIndividual] = useState<any>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // 从 data 中获取 machine_designer_full.json 的路径信息
  useEffect(() => {
    const fetchParetoData = async () => {
      if (!data) return;

      setLoading(true);
      setError(null);

      try {
        // 调用前端 API 路由获取 Pareto 前沿数据
        const response = await axios.get("/api/acmopv2/pareto-front", {
          params: {
            path2MachineDesignerFull: "machine_designer_full.json"
          }
        });

        setParetoData(response.data);

        // 默认选择最后一个评估的个体（索引最大的）
        if (response.data.allIndividuals && response.data.allIndividuals.length > 0) {
          const lastIndividual = response.data.allIndividuals[response.data.allIndividuals.length - 1];
          setSelectedIndividual(lastIndividual);
        }
      } catch (err: any) {
        console.error("Failed to fetch Pareto front data", err);
        setError(err.response?.data?.detail || err.message || "获取Pareto前沿数据失败");
      } finally {
        setLoading(false);
      }
    };

    fetchParetoData();
  }, [data]);

  // 处理个体选择变化
  const handleIndividualChange = (individualKey: string) => {
    if (!paretoData) return;
    
    const individual = paretoData.allIndividuals?.find(
      (ind: any) => ind.key === individualKey
    );
    
    if (individual) {
      setSelectedIndividual(individual);
    }
  };

  // 获取选中个体的完整数据
  const getSelectedIndividualFullData = () => {
    if (!selectedIndividual || !paretoData) return null;

    // 从 paretoFront 或 allIndividuals 中查找完整数据
    const fullData = paretoData.paretoFront?.find(
      (ind: any) => ind.key === selectedIndividual.key
    ) || paretoData.allIndividuals?.find(
      (ind: any) => ind.key === selectedIndividual.key
    );

    return fullData;
  };

  const selectedFullData = getSelectedIndividualFullData();
  const path2FEACsv = paretoData?.path2FEACsv || "";
  const individualIndex = selectedIndividual?.individual_index ?? selectedIndividual?.index;

  // 构建 CSV 路径（根据 individual_index）
  // path2FEACsv 格式类似: "C:\_Codes\ACMOP\backend\_default\SPMSM\csv/0/"
  // 需要替换最后的数字目录为 individual_index
  const csvPath = individualIndex !== undefined && individualIndex !== null && path2FEACsv
    ? path2FEACsv.replace(/\/\d+\/?$/, `/${individualIndex}/`)
    : path2FEACsv;

  // 准备 Pareto 前沿图表数据
  const paretoChartData = paretoData?.paretoFront?.map((ind: any) => ({
    objective1: ind.objectives?.f1 || 0,
    objective2: ind.objectives?.f2 || 0,
    objective3: ind.objectives?.f3 || 0,
    name: `Gen${ind.generation}-Ind${ind.individual_index}`,
    key: ind.key
  })) || [];

  // 准备 Pareto 前沿表格数据
  const paretoTableData = paretoData?.paretoFront || [];

  return (
    <div className="grid grid-cols-2 gap-6 h-full">
      {/* 左列：个体选择和结果可视化 */}
      <div className="space-y-4 overflow-y-auto">
        <Card>
          <CardHeader>
            <CardTitle>个体选择</CardTitle>
            <CardDescription>从已评估的个体中选择查看详细结果</CardDescription>
          </CardHeader>
          <CardContent>
            {loading ? (
              <div className="text-center py-4 text-muted-foreground">加载中...</div>
            ) : error ? (
              <div className="text-center py-4 text-destructive">{error}</div>
            ) : paretoData?.allIndividuals ? (
              <Select
                value={selectedIndividual?.key || ""}
                onValueChange={handleIndividualChange}
              >
                <SelectTrigger>
                  <SelectValue placeholder="选择个体" />
                </SelectTrigger>
                <SelectContent>
                  {paretoData.allIndividuals.map((ind: any) => (
                    <SelectItem key={ind.key} value={ind.key}>
                      {ind.is_pareto && "⭐ "}
                      Gen{ind.generation}-Ind{ind.individual_index}
                      {ind.project_name && ` (${ind.project_name})`}
                    </SelectItem>
                  ))}
                </SelectContent>
              </Select>
            ) : (
              <div className="text-center py-4 text-muted-foreground">暂无个体数据</div>
            )}
          </CardContent>
        </Card>

        {/* 选中个体的 JSON 结果 */}
        {selectedFullData && (
          <Card>
            <CardHeader>
              <CardTitle>性能指标</CardTitle>
              <CardDescription>
                个体 {selectedFullData.generation !== undefined 
                  ? `Gen${selectedFullData.generation}-Ind${selectedFullData.individual_index}`
                  : selectedFullData.key}
              </CardDescription>
            </CardHeader>
            <CardContent>
              <div className="space-y-4">
                {/* 目标函数值 */}
                {selectedFullData.objectives && (
                  <div>
                    <h4 className="font-semibold mb-2">目标函数值</h4>
                    <div className="space-y-1 text-sm">
                      <div className="flex justify-between">
                        <span className="text-muted-foreground">f1 ({paretoData?.objectives?.[0] || "目标1"}):</span>
                        <span className="font-medium">{selectedFullData.objectives.f1?.toFixed(4) || "N/A"}</span>
                      </div>
                      <div className="flex justify-between">
                        <span className="text-muted-foreground">f2 ({paretoData?.objectives?.[1] || "目标2"}):</span>
                        <span className="font-medium">{selectedFullData.objectives.f2?.toFixed(4) || "N/A"}</span>
                      </div>
                      {selectedFullData.objectives.f3 !== undefined && (
                        <div className="flex justify-between">
                          <span className="text-muted-foreground">f3 ({paretoData?.objectives?.[2] || "目标3"}):</span>
                          <span className="font-medium">{selectedFullData.objectives.f3?.toFixed(4) || "N/A"}</span>
                        </div>
                      )}
                    </div>
                  </div>
                )}

                {/* 设计参数 */}
                {selectedFullData.parameters && Object.keys(selectedFullData.parameters).length > 0 && (
                  <div>
                    <h4 className="font-semibold mb-2">设计参数</h4>
                    <div className="space-y-1 text-sm max-h-48 overflow-y-auto">
                      {Object.entries(selectedFullData.parameters).map(([key, value]: [string, any]) => (
                        <div key={key} className="flex justify-between">
                          <span className="text-muted-foreground">{key}:</span>
                          <span className="font-medium">{typeof value === 'number' ? value.toFixed(4) : String(value)}</span>
                        </div>
                      ))}
                    </div>
                  </div>
                )}

                {/* 其他性能指标 */}
                {selectedFullData.performance && Object.keys(selectedFullData.performance).length > 0 && (
                  <div>
                    <h4 className="font-semibold mb-2">其他性能指标</h4>
                    <div className="space-y-1 text-sm max-h-48 overflow-y-auto">
                      {Object.entries(selectedFullData.performance).slice(0, 10).map(([key, value]: [string, any]) => (
                        <div key={key} className="flex justify-between">
                          <span className="text-muted-foreground">{key}:</span>
                          <span className="font-medium">{typeof value === 'number' ? value.toFixed(4) : String(value)}</span>
                        </div>
                      ))}
                    </div>
                  </div>
                )}
              </div>
            </CardContent>
          </Card>
        )}

        {/* CSV 结果可视化 */}
        {csvPath && (
          <Card>
            <CardHeader>
              <CardTitle>FEA 仿真结果</CardTitle>
              <CardDescription>CSV 数据可视化</CardDescription>
            </CardHeader>
            <CardContent>
              <CsvVisualizer path2FEACsv={csvPath} />
            </CardContent>
          </Card>
        )}
      </div>

      {/* 右列：Pareto前沿和优化配置 */}
      <div className="space-y-4 overflow-y-auto">
        {/* Pareto前沿图表 */}
        {paretoChartData.length > 0 && (
          <Card>
            <CardHeader>
              <CardTitle>Pareto前沿</CardTitle>
              <CardDescription>多目标优化的Pareto最优解集</CardDescription>
            </CardHeader>
            <CardContent>
              <ResponsiveContainer width="100%" height={400}>
                <ScatterChart data={paretoChartData}>
                  <CartesianGrid strokeDasharray="3 3" />
                  <XAxis 
                    dataKey="objective1" 
                    name={paretoData?.objectives?.[0] || "目标1"}
                    type="number"
                    label={{ value: paretoData?.objectives?.[0] || "目标1", position: "insideBottom", offset: -5 }}
                  />
                  <YAxis 
                    dataKey="objective2" 
                    name={paretoData?.objectives?.[1] || "目标2"}
                    type="number"
                    label={{ value: paretoData?.objectives?.[1] || "目标2", angle: -90, position: "insideLeft" }}
                  />
                  <Tooltip 
                    cursor={{ strokeDasharray: '3 3' }}
                    content={({ active, payload }) => {
                      if (active && payload && payload.length) {
                        const data = payload[0].payload;
                        return (
                          <div className="bg-background border border-border rounded-lg p-3 shadow-lg">
                            <p className="font-semibold">{data.name}</p>
                            <p className="text-sm">
                              {paretoData?.objectives?.[0] || "目标1"}: {data.objective1?.toFixed(4)}
                            </p>
                            <p className="text-sm">
                              {paretoData?.objectives?.[1] || "目标2"}: {data.objective2?.toFixed(4)}
                            </p>
                            {data.objective3 !== undefined && (
                              <p className="text-sm">
                                {paretoData?.objectives?.[2] || "目标3"}: {data.objective3?.toFixed(4)}
                              </p>
                            )}
                          </div>
                        );
                      }
                      return null;
                    }}
                  />
                  <Legend />
                  <Scatter 
                    name="Pareto解" 
                    data={paretoChartData} 
                    fill="#8884d8"
                  />
                </ScatterChart>
              </ResponsiveContainer>
            </CardContent>
          </Card>
        )}

        {/* Pareto前沿表格 */}
        {paretoTableData.length > 0 && (
          <Card>
            <CardHeader>
              <CardTitle>Pareto前沿个体列表</CardTitle>
              <CardDescription>所有Pareto最优解的详细信息</CardDescription>
            </CardHeader>
            <CardContent>
              <div className="max-h-96 overflow-y-auto">
                <Table>
                  <TableHeader>
                    <TableRow>
                      <TableHead>个体</TableHead>
                      <TableHead>{paretoData?.objectives?.[0] || "目标1"}</TableHead>
                      <TableHead>{paretoData?.objectives?.[1] || "目标2"}</TableHead>
                      {paretoData?.objectives?.[2] && (
                        <TableHead>{paretoData.objectives[2]}</TableHead>
                      )}
                    </TableRow>
                  </TableHeader>
                  <TableBody>
                    {paretoTableData.map((ind: any) => (
                      <TableRow 
                        key={ind.key}
                        className={selectedIndividual?.key === ind.key ? "bg-muted" : ""}
                        onClick={() => handleIndividualChange(ind.key)}
                        style={{ cursor: "pointer" }}
                      >
                        <TableCell>
                          Gen{ind.generation}-Ind{ind.individual_index}
                        </TableCell>
                        <TableCell>
                          {ind.objectives?.f1?.toFixed(4) || "N/A"}
                        </TableCell>
                        <TableCell>
                          {ind.objectives?.f2?.toFixed(4) || "N/A"}
                        </TableCell>
                        {paretoData?.objectives?.[2] && (
                          <TableCell>
                            {ind.objectives?.f3?.toFixed(4) || "N/A"}
                          </TableCell>
                        )}
                      </TableRow>
                    ))}
                  </TableBody>
                </Table>
              </div>
            </CardContent>
          </Card>
        )}

        {/* 优化配置 */}
        {paretoData?.mooConfig && (
          <Card>
            <CardHeader>
              <CardTitle>优化配置</CardTitle>
              <CardDescription>多目标优化配置参数</CardDescription>
            </CardHeader>
            <CardContent>
              <div className="space-y-2">
                <div className="text-sm">
                  <span className="text-muted-foreground">配置名称:</span>{" "}
                  <span className="font-medium">{paretoData.select_fea_config_dict || "N/A"}</span>
                </div>
                <div className="space-y-1 text-sm">
                  {Object.entries(paretoData.mooConfig).map(([key, value]: [string, any]) => (
                    <div key={key} className="flex justify-between">
                      <span className="text-muted-foreground">{key}:</span>
                      <span className="font-medium">{String(value)}</span>
                    </div>
                  ))}
                </div>
              </div>
            </CardContent>
          </Card>
        )}

        {!paretoData && !loading && !error && (
          <Card>
            <CardContent className="py-8 text-center text-muted-foreground">
              暂无优化结果数据
            </CardContent>
          </Card>
        )}
      </div>
    </div>
  );
}
