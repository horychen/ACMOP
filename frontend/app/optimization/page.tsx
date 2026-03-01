"use client"

import { useState, useEffect, useCallback, useMemo, useRef } from "react";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { ScatterChart, Scatter, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer } from "recharts";
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "@/components/ui/table";
import { RefreshCw, Loader2, AlertCircle, FileJson, FileText, Wrench } from "lucide-react";
import { useRouter } from "next/navigation";
import CsvVisualizer from "@/components/CsvVisualizer";
import { ParetoFront2p5D } from "@/components/ParetoFront2p5D";
import { ParameterHistogram } from "@/components/ParameterHistogram";
import axios from "axios";

const BACKEND_URL = process.env.NEXT_PUBLIC_BACKEND_URL || "http://localhost:8000";
const SELECTED_INDIVIDUAL_KEY = "fine-tune-selected-individual";
const SELECTED_FOLDER_KEY = "fine-tune-selected-folder";
const OPTIMIZATION_SELECTED_FOLDER_KEY = "optimization-selected-folder";
const OPTIMIZATION_SELECTED_INDIVIDUAL_KEY = "optimization-selected-individual";

export default function OptimizationPage() {
  const router = useRouter();
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
  const [selectedCsvFileName, setSelectedCsvFileName] = useState<string>(""); // 记忆选中的 CSV 文件名
  const [geometryPdfPath, setGeometryPdfPath] = useState<string>(""); // 几何 PDF 文件路径
  const [loadingGeometry, setLoadingGeometry] = useState<boolean>(false); // 加载几何 PDF 的状态
  const [highlightedKeys, setHighlightedKeys] = useState<string[]>([]); // 高亮的个体key列表（用于直方图选中）
  const [csvDataCache, setCsvDataCache] = useState<Map<string, Map<string, any>>>(new Map()); // CSV数据缓存：Map<individualPath, Map<fileName, parsedData>>
  const csvDataCacheRef = useRef<Map<string, Map<string, any>>>(new Map()); // CSV数据缓存引用，用于在useCallback中访问
  const [loadingCsvCache, setLoadingCsvCache] = useState<boolean>(false); // 加载CSV缓存的状态

  // 同步ref和state
  useEffect(() => {
    csvDataCacheRef.current = csvDataCache;
  }, [csvDataCache]);

  // 从 localStorage 加载上次选择
  useEffect(() => {
    if (typeof window !== "undefined") {
      const savedFolder = localStorage.getItem(OPTIMIZATION_SELECTED_FOLDER_KEY);
      const savedIndividual = localStorage.getItem(OPTIMIZATION_SELECTED_INDIVIDUAL_KEY);

      if (savedFolder) {
        setSelectedFolder(savedFolder);
      }

      if (savedIndividual) {
        try {
          const individual = JSON.parse(savedIndividual);
          setSelectedIndividual(individual);
        } catch (e) {
          console.error("Failed to parse saved individual:", e);
        }
      }
    }
  }, []);

  // 保存选中的文件夹到 localStorage
  useEffect(() => {
    if (selectedFolder && typeof window !== "undefined") {
      localStorage.setItem(OPTIMIZATION_SELECTED_FOLDER_KEY, selectedFolder);
    }
  }, [selectedFolder]);

  // 保存选中的个体到 localStorage
  useEffect(() => {
    if (selectedIndividual && typeof window !== "undefined") {
      localStorage.setItem(OPTIMIZATION_SELECTED_INDIVIDUAL_KEY, JSON.stringify(selectedIndividual));
    }
  }, [selectedIndividual]);

  // 加载文件夹列表
  useEffect(() => {
    const fetchFolders = async () => {
      setLoadingFolders(true);
      try {
        const response = await axios.get(`${BACKEND_URL}/api/acmopv2/list-optimization-folders`);
        setFolders(response.data.folders || []);
      } catch (err: any) {
        console.error("Failed to fetch folders", err);
        setError(err.response?.data?.error || "获取文件夹列表失败");
      } finally {
        setLoadingFolders(false);
      }
    };

    fetchFolders();
  }, []);

  // 加载个体对应的所有CSV文件到缓存
  const loadCsvDataForIndividual = useCallback(async (csvPath: string) => {
    if (!csvPath) return;

    // 检查缓存中是否已有该路径的数据
    if (csvDataCacheRef.current.has(csvPath)) {
      console.log('CSV data already cached for:', csvPath);
      return;
    }

    setLoadingCsvCache(true);
    try {

      // 先获取文件列表
      const listResponse = await axios.get(`${BACKEND_URL}/api/results/csv/list-from-path`, {
        params: { path: csvPath }
      });

      const files = listResponse.data || [];
      console.log(`Loading ${files.length} CSV files for individual:`, csvPath);

      // 并行加载所有CSV文件
      const csvDataMap = new Map<string, any>();
      const loadPromises = files.map(async (filename: string) => {
        try {
          const contentResponse = await axios.get(`${BACKEND_URL}/api/results/csv/content-from-path`, {
            params: { path: csvPath, filename }
          });

          const csvContent = contentResponse.data.content;
          const lines = csvContent.split('\n');
          const headerIndex = lines.findIndex((line: string) => line.trim().startsWith('Time(s)'));

          if (headerIndex === -1) {
            console.warn(`Could not find data header 'Time(s)' in CSV file: ${filename}`);
            return null;
          }

          const cleanCsvContent = lines.slice(headerIndex).join('\n');

          // 使用d3解析CSV（需要在组件中导入）
          const d3 = await import('d3');
          const parsedData = d3.csvParse(cleanCsvContent);

          if (parsedData.length > 0) {
            csvDataMap.set(filename, parsedData);
          }
        } catch (err: any) {
          console.error(`Failed to load CSV file ${filename}:`, err);
        }
      });

      await Promise.all(loadPromises);

      // 更新缓存
      setCsvDataCache(prev => {
        const newCache = new Map(prev);
        newCache.set(csvPath, csvDataMap);
        csvDataCacheRef.current = newCache; // 同步更新ref
        return newCache;
      });

      console.log(`Successfully cached ${csvDataMap.size} CSV files for:`, csvPath);
    } catch (err: any) {
      console.error("Failed to load CSV data cache", err);
    } finally {
      setLoadingCsvCache(false);
    }
  }, [BACKEND_URL]);

  // 加载 Pareto 前沿数据（只在用户更换文件夹时调用）
  const fetchParetoData = useCallback(async (clearIndividual: boolean = true) => {
    if (!selectedFolder) return;

    setLoading(true);
    setError(null);

    // 只有在明确要求时才清空个体选择（例如切换文件夹时）
    if (clearIndividual) {
      setSelectedIndividual(null);
      setGeometryPdfPath("");
      setCurrentCsvFilePath("");
      setSelectedCsvFileName("");
      setHighlightedKeys([]);
      setZFilter(undefined);

      // 清空CSV缓存（切换文件夹时）
      setCsvDataCache(new Map());
    }

    try {
      const response = await axios.get(`${BACKEND_URL}/api/acmopv2/pareto-front`, {
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
        setSwarmDataPath(`backend_v2/output/_default/${selectedFolder}/SwarmData.json`);
      }

      // 如果之前有保存的个体，尝试恢复选择（从 localStorage 读取）
      if (typeof window !== "undefined") {
        const savedIndividualStr = localStorage.getItem(OPTIMIZATION_SELECTED_INDIVIDUAL_KEY);
        if (savedIndividualStr) {
          try {
            const savedIndividual = JSON.parse(savedIndividualStr);
            const savedIndividualKey = savedIndividual.key;
            const fullData = response.data.paretoFront?.find(
              (ind: any) => ind.key === savedIndividualKey
            ) || response.data.allIndividuals?.find(
              (ind: any) => ind.key === savedIndividualKey
            );

            if (fullData) {
              setSelectedIndividual(fullData);
              // 自动生成几何 PDF
              setTimeout(() => {
                generateGeometryPdf(fullData);
              }, 100);
            }
          } catch (e) {
            console.error("Failed to restore individual:", e);
          }
        }
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
  }, [selectedFolder, loadCsvDataForIndividual]);

  // 当 selectedFolder 变化时加载数据（只在用户手动更换文件夹时读取）
  useEffect(() => {
    if (selectedFolder) {
      // 检查是否是首次加载（从 localStorage 恢复）
      const isInitialLoad = !paretoData;
      // 首次加载时不清空个体，后续切换文件夹时清空
      fetchParetoData(!isInitialLoad);
    }
  }, [selectedFolder]); // 移除 fetchParetoData 依赖，避免不必要的重新加载

  // 处理个体选择变化
  const handleIndividualChange = useCallback((individualKey: string) => {
    if (!paretoData) return;

    const individual = paretoData.allIndividuals?.find(
      (ind: any) => ind.key === individualKey
    );

    if (individual) {
      setSelectedIndividual(individual);
      // 生成几何 PDF
      generateGeometryPdf(individual);

      // 计算CSV路径并加载所有CSV文件到缓存
      const individualIndex = individual.individual_index ?? individual.index;
      const currentPath2FEACsv = paretoData?.path2FEACsv || "";
      if (individualIndex !== undefined && individualIndex !== null && currentPath2FEACsv) {
        const normalizedPath = currentPath2FEACsv.replace(/\\/g, '/');
        let newPath: string;
        if (/\/\d+\/?$/.test(normalizedPath)) {
          newPath = normalizedPath.replace(/\/\d+\/?$/, `/${individualIndex}/`);
        } else {
          newPath = normalizedPath.endsWith('/')
            ? `${normalizedPath}${individualIndex}/`
            : `${normalizedPath}/${individualIndex}/`;
        }
        loadCsvDataForIndividual(newPath);
      }
    }
  }, [paretoData, loadCsvDataForIndividual]);

  // 生成几何 PDF
  const generateGeometryPdf = async (individual: any) => {
    if (!selectedFolder || !individual) return;

    const index = individual.individual_index ?? individual.index;
    if (index === undefined || index === null) return;

    setLoadingGeometry(true);
    setGeometryPdfPath("");

    try {
      const response = await axios.get(`${BACKEND_URL}/api/acmopv2/generate-geometry-pdf`, {
        params: {
          folderName: selectedFolder,
          individualIndex: index
        }
      });

      console.log("PDF generation response:", response.data);

      if (response.data?.pdfPath) {
        // 构建完整的 PDF URL（通过后端）
        // 添加 PDF 查看参数：隐藏侧边栏(navpanes=0)，隐藏工具栏(toolbar=0)，默认缩放400%(zoom=400)
        const pdfUrl = `${BACKEND_URL}/api/results/pdf?path=${encodeURIComponent(response.data.pdfPath)}#navpanes=0&toolbar=0&zoom=400`;
        console.log("PDF URL:", pdfUrl);
        setGeometryPdfPath(pdfUrl);
      } else {
        console.error("No pdfPath in response:", response.data);
      }
    } catch (err: any) {
      console.error("Failed to generate geometry PDF", err);
      console.error("Error details:", err.response?.data || err.message);
      setGeometryPdfPath("");
      // 显示错误提示
      if (err.response?.status === 404) {
        setError(`无法生成几何 PDF: ${err.response?.data?.detail || "API 端点未找到"}`);
      } else {
        setError(`生成几何 PDF 失败: ${err.response?.data?.detail || err.message}`);
      }
    } finally {
      setLoadingGeometry(false);
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
    // 处理 Windows 路径和 Unix 路径
    const normalizedPath = path2FEACsv.replace(/\\/g, '/');

    // 如果路径已经以数字结尾，替换它；否则追加 individual_index
    let newPath: string;
    if (/\/\d+\/?$/.test(normalizedPath)) {
      // 路径末尾有数字，替换它
      newPath = normalizedPath.replace(/\/\d+\/?$/, `/${individualIndex}/`);
    } else {
      // 路径末尾没有数字，追加 individual_index
      newPath = normalizedPath.endsWith('/')
        ? `${normalizedPath}${individualIndex}/`
        : `${normalizedPath}/${individualIndex}/`;
    }

    console.log('CSV Path calculation:', {
      original: path2FEACsv,
      normalized: normalizedPath,
      individualIndex,
      newPath,
      selectedIndividualKey: selectedIndividual.key
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
      // 刷新时不清空个体选择
      fetchParetoData(false);
    }
  };

  // 保存选中的个体到 localStorage 并跳转到 fine-tune 页面
  const handleSendToFineTune = () => {
    if (!selectedFullData || !selectedFolder) return;

    // 保存选中的个体和文件夹
    localStorage.setItem(SELECTED_INDIVIDUAL_KEY, JSON.stringify(selectedFullData));
    localStorage.setItem(SELECTED_FOLDER_KEY, selectedFolder);

    // 跳转到 fine-tune 页面
    router.push("/fine-tune");
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
                <CardContent className="space-y-4">
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
                  {selectedFullData && (
                    <Button
                      onClick={handleSendToFineTune}
                      className="w-full"
                      variant="outline"
                    >
                      <Wrench className="h-4 w-4 mr-2" />
                      发送到 Fine-tune
                    </Button>
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

                      {/* 几何图形 PDF */}
                      <div>
                        <h4 className="font-semibold mb-2">几何图形</h4>
                        <div className="border border-border rounded-lg bg-muted/30 p-4 min-h-[400px] flex items-center justify-center">
                          {loadingGeometry ? (
                            <div className="flex flex-col items-center space-y-2">
                              <Loader2 className="w-6 h-6 animate-spin text-primary" />
                              <span className="text-sm text-muted-foreground">正在生成几何图形...</span>
                            </div>
                          ) : geometryPdfPath ? (
                            <iframe
                              src={geometryPdfPath}
                              className="w-full h-[600px] border-0 rounded"
                              title="几何图形 PDF"
                            />
                          ) : (
                            <div className="text-center text-muted-foreground">
                              <p className="text-sm">点击 Pareto 前沿上的标记以生成几何图形</p>
                            </div>
                          )}
                        </div>
                      </div>
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
                  <CardContent className="h-[600px]">
                    <CsvVisualizer
                      key={`${csvPath}-${selectedIndividual.key}`}
                      path2FEACsv={csvPath}
                      selectedFile={selectedCsvFileName}
                      onFileSelect={setSelectedCsvFileName}
                      onCurrentFileChange={setCurrentCsvFilePath}
                      csvDataCache={csvDataCache.get(csvPath)}
                      loadingCache={loadingCsvCache}
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
                      highlightedKeys={highlightedKeys}
                    />
                  </CardContent>
                </Card>
              )}

              {/* Parameter Histogram */}
              {paretoData?.paretoFront && paretoData.paretoFront.length > 0 && (
                <ParameterHistogram
                  individuals={paretoData.paretoFront}
                  parameterName="stator_tooth_span_angle"
                  f3Threshold={5}
                  onBinSelect={setHighlightedKeys}
                />
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
                ? "未找到优化文件夹，请确保 backend_v2/output/_default 目录下有包含 SwarmData.json 的文件夹"
                : "请从上方选择器中选择一个优化文件夹"}
            </CardContent>
          </Card>
        )}
      </div>
    </div>
  );
}

export const dynamic = 'force-dynamic';
