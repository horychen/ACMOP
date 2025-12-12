"use client"

import { useState, useEffect, useCallback, useMemo } from "react";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { ScatterChart, Scatter, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer } from "recharts";
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "@/components/ui/table";
import { RefreshCw, Loader2, AlertCircle, FileJson, FileText } from "lucide-react";
import CsvVisualizer from "@/components/CsvVisualizer";
import { ParetoFront2p5D } from "@/components/ParetoFront2p5D";
import axios from "axios";

export default function OptimizationPage() {
  const [folders, setFolders] = useState<string[]>([]);
  const [selectedFolder, setSelectedFolder] = useState<string>("");
  const [paretoData, setParetoData] = useState<any>(null);
  const [selectedIndividual, setSelectedIndividual] = useState<any>(null);
  const [loading, setLoading] = useState(false);
  const [loadingFolders, setLoadingFolders] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [zFilter, setZFilter] = useState<number | undefined>(undefined);
  const [swarmDataPath, setSwarmDataPath] = useState<string>("");
  const [currentCsvFilePath, setCurrentCsvFilePath] = useState<string>("");

  // 加载文件夹列表
  useEffect(() => {
    const fetchFolders = async () => {
      setLoadingFolders(true);
      try {
        const response = await axios.get("/api/acmopv2/list-optimization-folders");
        setFolders(response.data.folders || []);
        if (response.data.folders && response.data.folders.length > 0 && !selectedFolder) {
          setSelectedFolder(response.data.folders[0]);
        }
      } catch (err: any) {
        console.error("Failed to fetch folders", err);
        setError(err.response?.data?.error || "获取文件夹列表失败");
      } finally {
        setLoadingFolders(false);
      }
    };

    fetchFolders();
  }, []);

  // 加载 Pareto 前沿数据
  const fetchParetoData = useCallback(async () => {
    if (!selectedFolder) return;

    setLoading(true);
    setError(null);

    try {
      const response = await axios.get("/api/acmopv2/pareto-front", {
        params: {
          folderName: selectedFolder
        }
      });

      setParetoData(response.data);
      
      // 设置 SwarmData.json 路径
      if (response.data.path2SwarmData) {
        const fullPath = `${response.data.path2SwarmData}/SwarmData.json`;
        setSwarmDataPath(fullPath);
      } else if (selectedFolder) {
        setSwarmDataPath(`backend/_default/${selectedFolder}/SwarmData.json`);
      }

      // 默认选择最后一个评估的个体（索引最大的）
      if (response.data.allIndividuals && response.data.allIndividuals.length > 0) {
        const lastIndividual = response.data.allIndividuals[response.data.allIndividuals.length - 1];
        setSelectedIndividual(lastIndividual);
      }
    } catch (err: any) {
      console.error("Failed to fetch Pareto front data", err);
      const errorMessage = err.response?.data?.detail || err.response?.data?.error || err.message || "获取Pareto前沿数据失败";
      console.error("Error details:", JSON.stringify(err.response?.data, null, 2));
      setError(errorMessage);
      setParetoData(null);
      setSwarmDataPath("");
    } finally {
      setLoading(false);
    }
  }, [selectedFolder]);

  // 当 selectedFolder 变化时加载数据
  useEffect(() => {
    fetchParetoData();
  }, [fetchParetoData]);

  // 页面可见性变化时重新加载数据（刷新页面或切换标签页后回来时）
  useEffect(() => {
    const handleVisibilityChange = () => {
      if (document.visibilityState === 'visible' && selectedFolder) {
        // 页面变为可见时重新加载数据
        fetchParetoData();
      }
    };

    // 监听窗口焦点变化（页面刷新后）
    const handleFocus = () => {
      if (selectedFolder) {
        fetchParetoData();
      }
    };

    document.addEventListener('visibilitychange', handleVisibilityChange);
    window.addEventListener('focus', handleFocus);

    return () => {
      document.removeEventListener('visibilitychange', handleVisibilityChange);
      window.removeEventListener('focus', handleFocus);
    };
  }, [selectedFolder, fetchParetoData]);

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
  const csvPath = useMemo(() => {
    if (!selectedIndividual || !path2FEACsv) {
      return "";
    }
    if (individualIndex === undefined || individualIndex === null) {
      return path2FEACsv;
    }
    // 替换路径末尾的数字目录为 individual_index
    const newPath = path2FEACsv.replace(/\/\d+\/?$/, `/${individualIndex}/`);
    console.log('CSV Path calculation:', {
      original: path2FEACsv,
      individualIndex,
      newPath,
      selectedIndividual: selectedIndividual.key
    });
    return newPath;
  }, [path2FEACsv, individualIndex, selectedIndividual]);

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

  const handleRefresh = () => {
    if (selectedFolder) {
      fetchParetoData();
    }
  };

  return (
    <div className="h-[calc(100vh-4rem)] flex flex-col">
      {/* 状态栏（置顶） */}
      {(swarmDataPath || currentCsvFilePath) && (
        <div className="border-primary/50 bg-primary/5 flex-shrink-0 border-b px-4 py-1.5">
          <div className="flex items-center space-x-4 text-xs">
            {swarmDataPath && (
              <div className="flex items-center space-x-2">
                <FileJson className="h-3 w-3 text-primary flex-shrink-0" />
                <span className="text-muted-foreground whitespace-nowrap">正在读取:</span>
                <span className="font-mono text-primary font-medium truncate">{swarmDataPath}</span>
                {loading && (
                  <Loader2 className="h-3 w-3 animate-spin text-primary ml-1 flex-shrink-0" />
                )}
              </div>
            )}
            {currentCsvFilePath && (
              <div className="flex items-center space-x-2 border-l pl-4">
                <FileText className="h-3 w-3 text-primary flex-shrink-0" />
                <span className="text-muted-foreground whitespace-nowrap">当前CSV:</span>
                <span className="font-mono text-primary font-medium truncate">{currentCsvFilePath}</span>
              </div>
            )}
          </div>
        </div>
      )}

      <div className="flex-1 overflow-hidden flex flex-col space-y-6">
        {/* Header */}
        <div className="flex items-center justify-between flex-shrink-0 pt-6">
          <div>
            <h1 className="text-3xl font-bold tracking-tight">多目标优化监控</h1>
            <p className="text-muted-foreground mt-1">
              监控优化过程，查看 Pareto 前沿和个体性能
            </p>
          </div>
          <div className="flex items-center gap-4">
            <Select
              value={selectedFolder}
              onValueChange={setSelectedFolder}
              disabled={loadingFolders || folders.length === 0}
            >
              <SelectTrigger className="w-[300px]">
                <SelectValue placeholder="选择优化文件夹..." />
              </SelectTrigger>
              <SelectContent>
                {folders.map((folder) => (
                  <SelectItem key={folder} value={folder}>
                    {folder}
                  </SelectItem>
                ))}
              </SelectContent>
            </Select>
            <Button onClick={handleRefresh} disabled={loading || !selectedFolder}>
              <RefreshCw className={`h-4 w-4 mr-2 ${loading ? 'animate-spin' : ''}`} />
              刷新数据
            </Button>
          </div>
        </div>

        {/* Error Display */}
        {error && (
          <Card className="border-destructive flex-shrink-0">
            <CardContent className="pt-6">
              <div className="flex items-center space-x-2 text-destructive">
                <AlertCircle className="h-5 w-5" />
                <span>{error}</span>
              </div>
            </CardContent>
          </Card>
        )}

        {/* Main Content - Two Column Layout */}
        {selectedFolder && (
        <div className="flex-1 grid grid-cols-2 gap-6 overflow-hidden min-h-0">
          {/* Left Column: Individual Selection and Results */}
          <div className="space-y-4 overflow-y-auto">
            <Card>
              <CardHeader>
                <CardTitle>个体选择</CardTitle>
                <CardDescription>从已评估的个体中选择查看详细结果</CardDescription>
              </CardHeader>
              <CardContent>
                {loading ? (
                  <div className="text-center py-4 text-muted-foreground">
                    <Loader2 className="h-5 w-5 animate-spin mx-auto mb-2" />
                    加载中...
                  </div>
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

            {/* Selected Individual JSON Results */}
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

            {/* CSV Results Visualization */}
            {csvPath && selectedIndividual && (
              <Card>
                <CardHeader>
                  <CardTitle>FEA 仿真结果</CardTitle>
                  <CardDescription>
                    CSV 数据可视化 - {selectedIndividual.generation !== undefined 
                      ? `Gen${selectedIndividual.generation}-Ind${selectedIndividual.individual_index ?? selectedIndividual.index}`
                      : selectedIndividual.key}
                  </CardDescription>
                </CardHeader>
                <CardContent>
                  <CsvVisualizer 
                    key={`${csvPath}-${selectedIndividual.key}`} 
                    path2FEACsv={csvPath}
                    onCurrentFileChange={setCurrentCsvFilePath}
                  />
                </CardContent>
              </Card>
            )}
          </div>

          {/* Right Column: Pareto Front Information */}
          <div className="space-y-4 overflow-y-auto">
            {/* Pareto Front Chart - 2.5D Visualization */}
            {paretoData?.allIndividuals && paretoData.allIndividuals.length > 0 && (
              <Card>
                <CardHeader>
                  <CardTitle>Pareto前沿（三目标可视化）</CardTitle>
                  <CardDescription>
                    多目标优化的Pareto最优解集，颜色表示第三个目标函数值（f3）
                  </CardDescription>
                </CardHeader>
                <CardContent>
                  <ParetoFront2p5D
                    individuals={paretoData.allIndividuals.map((ind: any) => ({
                      ...ind,
                      f1: ind.objectives?.f1 ?? ind.f1 ?? 0,
                      f2: ind.objectives?.f2 ?? ind.f2 ?? 0,
                      f3: ind.objectives?.f3 ?? ind.f3 ?? 0,
                    }))}
                    objectives={paretoData.objectives || ["f1", "f2", "f3"]}
                    comp={[0, 1]} // f1 vs f2，f3 作为颜色
                    upToRankNo={1}
                    zFilter={zFilter}
                    onZFilterChange={setZFilter}
                    onIndividualSelect={handleIndividualChange}
                  />
                </CardContent>
              </Card>
            )}

            {/* Pareto Front Table */}
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

            {/* Optimization Configuration */}
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
        )}

        {/* Empty State */}
        {!selectedFolder && !loadingFolders && (
          <Card>
            <CardContent className="py-8 text-center text-muted-foreground">
              {folders.length === 0 
                ? "未找到优化文件夹，请确保 _default 目录下有包含 SwarmData.json 的文件夹"
                : "请从上方选择器中选择一个优化文件夹"}
            </CardContent>
          </Card>
        )}
      </div>
    </div>
  );
}

export const dynamic = 'force-dynamic';
