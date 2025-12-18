"use client"

import { useState, useEffect, useCallback, useRef } from "react";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/alert";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Loader2, AlertCircle, ArrowLeft, RefreshCw, Copy, Check } from "lucide-react";
import { useRouter } from "next/navigation";
import axios from "axios";
import CsvVisualizer from "@/components/CsvVisualizer";
import { ScrollArea } from "@/components/ui/scroll-area";

const BACKEND_URL = process.env.NEXT_PUBLIC_BACKEND_URL || "http://localhost:8000";
const SELECTED_INDIVIDUAL_KEY = "fine-tune-selected-individual";
const SELECTED_FOLDER_KEY = "fine-tune-selected-folder";

interface ParameterInfo {
  name: string; // 参数名（Parameter.name）
  fieldName?: string; // 字段名（用于区分）
  displayName?: string;
  type: "fixed" | "free" | "derived";
  value: number;
  bounds: [number, number] | null;
  unit: string;
  comment?: string;
}

export default function FineTunePage() {
  const router = useRouter();
  const [selectedIndividual, setSelectedIndividual] = useState<any>(null);
  const [selectedFolder, setSelectedFolder] = useState<string>("");
  const [paretoData, setParetoData] = useState<any>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [geometryPdfPath, setGeometryPdfPath] = useState<string>("");
  const [loadingGeometry, setLoadingGeometry] = useState<boolean>(false);
  const [csvPath, setCsvPath] = useState<string>("");
  const [selectedCsvFileName, setSelectedCsvFileName] = useState<string>("");
  const [currentCsvFilePath, setCurrentCsvFilePath] = useState<string>("");
  
  // 参数编辑相关状态
  const [parametersInfo, setParametersInfo] = useState<Record<string, ParameterInfo>>({});
  const [editableParameters, setEditableParameters] = useState<Record<string, number>>({});
  const [isLoadingParams, setIsLoadingParams] = useState(false);
  const updateTimeoutRef = useRef<NodeJS.Timeout | null>(null);
  
  // 复制状态
  const [copied, setCopied] = useState(false);

  // 从 localStorage 加载选中的个体
  useEffect(() => {
    if (typeof window !== "undefined") {
      const savedIndividual = localStorage.getItem(SELECTED_INDIVIDUAL_KEY);
      const savedFolder = localStorage.getItem(SELECTED_FOLDER_KEY);
      
      if (savedIndividual) {
        try {
          const individual = JSON.parse(savedIndividual);
          setSelectedIndividual(individual);
          console.log("Loaded individual from localStorage:", individual);
        } catch (e) {
          console.error("Failed to parse saved individual:", e);
        }
      }
      
      if (savedFolder) {
        setSelectedFolder(savedFolder);
        console.log("Loaded folder from localStorage:", savedFolder);
      }
    }
  }, []);

  // 加载参数信息
  useEffect(() => {
    const loadParametersInfo = async () => {
      setIsLoadingParams(true);
      try {
        const response = await axios.get(`${BACKEND_URL}/api/acmopv2/all-parameters`);
        const apiData = response.data;
        
        // 构建参数信息映射
        const paramsMap: Record<string, ParameterInfo> = {};
        
        // 合并所有类型的参数
        const allParams = [
          ...(apiData.parameters.fixed || []),
          ...(apiData.parameters.free || []),
          ...(apiData.parameters.derived || [])
        ];
        
        allParams.forEach((param: any) => {
          // API 返回的 param.name 是字段名，param.displayName 是参数名（Parameter.name）
          // x_denorm_dict 使用参数名（param.name，即 Parameter.name）作为键
          const fieldName = param.name; // 字段名（如 'mm_r_so'）
          const paramName = param.displayName || param.name; // 参数名（Parameter.name）
          
          const paramInfo: ParameterInfo = {
            name: paramName,
            fieldName: fieldName,
            displayName: param.displayName || param.name,
            type: param.type,
            unit: param.unit || "mm",
            value: param.value,
            bounds: param.bounds || null,
            comment: param.comment || ""
          };
          
          // 同时使用字段名和参数名作为 key，以便匹配
          paramsMap[fieldName] = paramInfo;
          if (paramName !== fieldName) {
            paramsMap[paramName] = paramInfo;
          }
        });
        
        setParametersInfo(paramsMap);
      } catch (err: any) {
        console.error("加载参数信息失败:", err);
      } finally {
        setIsLoadingParams(false);
      }
    };

    loadParametersInfo();
  }, []);

  // 初始化可编辑参数（从选中个体的参数中提取 free 类型参数）
  useEffect(() => {
    console.log("Initializing editable parameters:", {
      hasSelectedIndividual: !!selectedIndividual,
      hasParameters: !!selectedIndividual?.parameters,
      parametersCount: selectedIndividual?.parameters ? Object.keys(selectedIndividual.parameters).length : 0,
      parametersInfoCount: Object.keys(parametersInfo).length
    });
    
    if (selectedIndividual?.parameters && Object.keys(parametersInfo).length > 0) {
      const editable: Record<string, number> = {};
      const unmatchedParams: string[] = [];
      
      console.log("Checking parameters:", {
        individualParams: Object.keys(selectedIndividual.parameters),
        parametersInfoKeys: Object.keys(parametersInfo).slice(0, 10) // 只显示前10个
      });
      
      // 只包含 free 类型的参数
      Object.entries(selectedIndividual.parameters).forEach(([paramName, value]: [string, any]) => {
        // 尝试多种匹配方式
        let paramInfo = parametersInfo[paramName]; // 直接匹配
        
        // 如果直接匹配失败，尝试通过 displayName 匹配
        if (!paramInfo) {
          paramInfo = Object.values(parametersInfo).find(
            (info: any) => info.displayName === paramName || info.name === paramName
          ) as ParameterInfo | undefined;
        }
        
        if (paramInfo && paramInfo.type === "free") {
          editable[paramName] = typeof value === 'number' ? value : parseFloat(value) || 0;
          console.log(`Matched parameter: ${paramName} -> ${paramInfo.displayName || paramInfo.name} (free)`);
        } else {
          unmatchedParams.push(paramName);
          if (paramInfo) {
            console.log(`Parameter ${paramName} is not free type, it's ${paramInfo.type}`);
          } else {
            console.warn(`Parameter ${paramName} not found in parametersInfo`);
          }
        }
      });
      
      console.log("Initialized editable parameters:", Object.keys(editable).length, editable);
      if (unmatchedParams.length > 0) {
        console.warn("Unmatched parameters:", unmatchedParams);
      }
      setEditableParameters(editable);
    } else {
      console.warn("Cannot initialize editable parameters:", {
        hasSelectedIndividual: !!selectedIndividual,
        hasParameters: !!selectedIndividual?.parameters,
        parametersInfoCount: Object.keys(parametersInfo).length
      });
      // 即使没有参数，也设置空对象，以便后续生成 PDF
      setEditableParameters({});
    }
  }, [selectedIndividual, parametersInfo]);

  // 加载 Pareto 数据
  useEffect(() => {
    const fetchParetoData = async () => {
      if (!selectedFolder) return;
      
      setLoading(true);
      setError(null);
      
      try {
        const response = await axios.get(`${BACKEND_URL}/api/acmopv2/pareto-front`, {
          params: {
            folderName: selectedFolder
          }
        });
        
        setParetoData(response.data);
        
        // 如果已保存的个体存在，尝试在数据中找到完整信息
        // 优先从 localStorage 读取，如果没有则使用当前的 selectedIndividual
        const savedIndividualStr = typeof window !== "undefined" 
          ? localStorage.getItem(SELECTED_INDIVIDUAL_KEY) 
          : null;
        
        let individualToMatch = selectedIndividual;
        if (savedIndividualStr) {
          try {
            individualToMatch = JSON.parse(savedIndividualStr);
          } catch (e) {
            console.error("Failed to parse saved individual:", e);
          }
        }
        
        if (individualToMatch?.key) {
          // 优先从 paretoFront 中查找（包含完整的 parameters）
          let fullData = response.data.paretoFront?.find(
            (ind: any) => ind.key === individualToMatch.key
          );
          
          // 如果 paretoFront 中没有，从 allIndividuals 中查找
          if (!fullData) {
            fullData = response.data.allIndividuals?.find(
              (ind: any) => ind.key === individualToMatch.key
            );
          }
          
          // 如果还是没有找到，尝试使用保存的个体数据（可能包含 parameters）
          if (!fullData && individualToMatch.parameters) {
            fullData = individualToMatch;
          }
          
          if (fullData) {
            console.log("Found full data for individual:", fullData);
            console.log("Individual has parameters:", !!fullData.parameters, "Count:", fullData.parameters ? Object.keys(fullData.parameters).length : 0);
            setSelectedIndividual(fullData);
            // 更新 localStorage 中的完整数据
            if (typeof window !== "undefined") {
              localStorage.setItem(SELECTED_INDIVIDUAL_KEY, JSON.stringify(fullData));
            }
          } else {
            console.warn("Could not find individual in pareto data:", individualToMatch.key);
          }
        }
      } catch (err: any) {
        console.error("Failed to fetch Pareto front data", err);
        setError(err.response?.data?.detail || err.message || "获取数据失败");
      } finally {
        setLoading(false);
      }
    };

    if (selectedFolder) {
      fetchParetoData();
    }
  }, [selectedFolder]);

  // 生成几何 PDF（使用自定义参数）
  const generateGeometryPdfWithParams = useCallback(async (params: Record<string, number>) => {
    if (!selectedFolder || !selectedIndividual) return;
    
    const index = selectedIndividual.individual_index ?? selectedIndividual.index;
    if (index === undefined || index === null) return;
    
    setLoadingGeometry(true);
    setError(null);
    
    try {
      const response = await axios.post(`${BACKEND_URL}/api/acmopv2/generate-geometry-pdf-with-params`, {
        folderName: selectedFolder,
        individualIndex: index,
        parameters: params
      });
      
      if (response.data?.pdfPath) {
        const pdfUrl = `${BACKEND_URL}/api/results/pdf?path=${encodeURIComponent(response.data.pdfPath)}#navpanes=0&toolbar=0&zoom=400`;
        setGeometryPdfPath(pdfUrl);
      }
    } catch (err: any) {
      console.error("Failed to generate geometry PDF", err);
      setError(err.response?.data?.detail || err.message || "生成几何 PDF 失败");
      setGeometryPdfPath("");
    } finally {
      setLoadingGeometry(false);
    }
  }, [selectedFolder, selectedIndividual]);

  // 初始生成几何 PDF（当个体和参数都准备好时）
  useEffect(() => {
    console.log("Checking if should generate PDF:", {
      hasSelectedIndividual: !!selectedIndividual,
      hasSelectedFolder: !!selectedFolder,
      hasParetoData: !!paretoData,
      editableParamsCount: Object.keys(editableParameters).length,
      hasParameters: !!selectedIndividual?.parameters
    });
    
    if (selectedIndividual && selectedFolder && paretoData) {
      // 如果有可编辑参数，使用它们；否则使用个体的原始参数
      const paramsToUse = Object.keys(editableParameters).length > 0 
        ? editableParameters 
        : (selectedIndividual.parameters || {});
      
      console.log("Generating PDF with parameters:", Object.keys(paramsToUse).length, paramsToUse);
      
      // 使用参数生成 PDF（即使参数为空，后端也会使用默认值）
      if (Object.keys(paramsToUse).length > 0 || selectedIndividual.parameters) {
        generateGeometryPdfWithParams(paramsToUse);
      }
      
      // 计算 CSV 路径
      const individualIndex = selectedIndividual.individual_index ?? selectedIndividual.index;
      const path2FEACsv = paretoData?.path2FEACsv || "";
      
      if (individualIndex !== undefined && individualIndex !== null && path2FEACsv) {
        const normalizedPath = path2FEACsv.replace(/\\/g, '/');
        let newPath: string;
        if (/\/\d+\/?$/.test(normalizedPath)) {
          newPath = normalizedPath.replace(/\/\d+\/?$/, `/${individualIndex}/`);
        } else {
          newPath = normalizedPath.endsWith('/') 
            ? `${normalizedPath}${individualIndex}/`
            : `${normalizedPath}/${individualIndex}/`;
        }
        setCsvPath(newPath);
      }
    }
  }, [selectedIndividual, selectedFolder, paretoData, editableParameters, generateGeometryPdfWithParams]); // 初始加载时生成

  // 处理参数值变化（带防抖）
  const handleParameterChange = useCallback((paramName: string, value: string) => {
    const numValue = parseFloat(value);
    if (isNaN(numValue)) return;
    
    // 更新本地状态
    setEditableParameters(prev => ({
      ...prev,
      [paramName]: numValue
    }));
    
    // 清除之前的定时器
    if (updateTimeoutRef.current) {
      clearTimeout(updateTimeoutRef.current);
    }
    
    // 设置新的定时器（防抖：500ms）
    updateTimeoutRef.current = setTimeout(() => {
      const updatedParams = {
        ...editableParameters,
        [paramName]: numValue
      };
      generateGeometryPdfWithParams(updatedParams);
    }, 500);
  }, [editableParameters, generateGeometryPdfWithParams]);

  // 清理定时器
  useEffect(() => {
    return () => {
      if (updateTimeoutRef.current) {
        clearTimeout(updateTimeoutRef.current);
      }
    };
  }, []);

  // 获取选中个体的完整数据
  const getSelectedIndividualFullData = () => {
    if (!selectedIndividual || !paretoData) return selectedIndividual;

    const fullData = paretoData.paretoFront?.find(
      (ind: any) => ind.key === selectedIndividual.key
    ) || paretoData.allIndividuals?.find(
      (ind: any) => ind.key === selectedIndividual.key
    );

    return fullData || selectedIndividual;
  };

  const selectedFullData = getSelectedIndividualFullData();

  // 生成参数字典字符串（Python 格式）
  const generateParametersDict = useCallback(() => {
    // 使用当前的可编辑参数，如果没有则使用原始参数
    const paramsToUse = Object.keys(editableParameters).length > 0 
      ? editableParameters 
      : (selectedIndividual?.parameters || {});
    
    if (Object.keys(paramsToUse).length === 0) {
      return "prev_params = { # No parameters available\n}";
    }
    
    // 格式化参数字典
    const lines: string[] = [];
    lines.push("prev_params = { # Fine-tuned parameters");
    
    // 计算最长的键名长度，用于对齐注释
    const maxKeyLength = Math.max(...Object.keys(paramsToUse).map(k => k.length));
    
    Object.entries(paramsToUse).forEach(([key, value]) => {
      const paramInfo = parametersInfo[key];
      // 构建注释：如果有参数信息，显示类型；如果有注释，也显示
      let comment = "";
      if (paramInfo) {
        if (paramInfo.type === "free") {
          comment = " # free variable";
        } else if (paramInfo.comment) {
          comment = ` # ${paramInfo.comment}`;
        }
      }
      
      // 使用参数名作为键，值保留适当的小数位数
      const formattedValue = typeof value === 'number' 
        ? (value % 1 === 0 ? value.toFixed(0) : value.toFixed(5))
        : value;
      
      // 对齐格式：键名 + 空格对齐 + 值 + 逗号 + 注释
      const keyWithQuotes = `"${key}"`;
      const padding = " ".repeat(Math.max(1, maxKeyLength - key.length + 2));
      lines.push(`    ${keyWithQuotes}:${padding}${formattedValue},${comment}`);
    });
    
    lines.push("}");
    return lines.join("\n");
  }, [editableParameters, selectedIndividual, parametersInfo]);

  const parametersDictString = generateParametersDict();

  // 复制参数字典到剪贴板
  const handleCopyParameters = useCallback(async () => {
    try {
      await navigator.clipboard.writeText(parametersDictString);
      setCopied(true);
      setTimeout(() => setCopied(false), 2000);
    } catch (err) {
      console.error("Failed to copy:", err);
    }
  }, [parametersDictString]);

  // 调试信息
  useEffect(() => {
    console.log("Fine-tune page state:", {
      selectedIndividual: selectedIndividual ? { key: selectedIndividual.key, generation: selectedIndividual.generation, individual_index: selectedIndividual.individual_index } : null,
      selectedFolder,
      paretoData: paretoData ? "loaded" : "not loaded",
      editableParametersCount: Object.keys(editableParameters).length,
      parametersInfoCount: Object.keys(parametersInfo).length
    });
  }, [selectedIndividual, selectedFolder, paretoData, editableParameters, parametersInfo]);

  if (!selectedIndividual) {
    return (
      <div className="container mx-auto py-8">
        <Card>
          <CardHeader>
            <CardTitle>Fine-tune</CardTitle>
            <CardDescription>微调选中的个体</CardDescription>
          </CardHeader>
          <CardContent>
            <Alert>
              <AlertCircle className="h-4 w-4" />
              <AlertTitle>未选择个体</AlertTitle>
              <AlertDescription>
                请先在优化页面选择一个个体，然后点击"发送到 Fine-tune"按钮。
              </AlertDescription>
            </Alert>
            <div className="mt-4">
              <Button onClick={() => router.push("/optimization")}>
                <ArrowLeft className="h-4 w-4 mr-2" />
                前往优化页面
              </Button>
            </div>
          </CardContent>
        </Card>
      </div>
    );
  }

  // 获取可编辑参数的列表
  const editableParamsList = Object.entries(editableParameters)
    .map(([name, value]) => ({
      name,
      value,
      info: parametersInfo[name]
    }))
    .filter(item => item.info); // 只显示有信息的参数

  return (
    <div className="container mx-auto py-8 space-y-6">
      {/* Header */}
      <div className="flex items-center justify-between">
        <div>
          <h1 className="text-3xl font-bold tracking-tight">Fine-tune</h1>
          <p className="text-muted-foreground mt-1">
            微调选中的个体: {selectedFullData.generation !== undefined 
              ? `Gen${selectedFullData.generation}-Ind${selectedFullData.individual_index}`
              : selectedFullData.key}
          </p>
        </div>
        <div className="flex items-center gap-4">
          <Button variant="outline" onClick={() => router.push("/optimization")}>
            <ArrowLeft className="h-4 w-4 mr-2" />
            返回优化页面
          </Button>
          <Button onClick={() => {
            if (selectedFolder) {
              window.location.reload();
            }
          }} disabled={loading || !selectedFolder}>
            <RefreshCw className={`h-4 w-4 mr-2 ${loading ? 'animate-spin' : ''}`} />
            刷新数据
          </Button>
        </div>
      </div>

      {/* Error Display */}
      {error && (
        <Alert variant="destructive">
          <AlertCircle className="h-4 w-4" />
          <AlertTitle>错误</AlertTitle>
          <AlertDescription>{error}</AlertDescription>
        </Alert>
      )}

      {/* Main Content */}
      <div className="grid grid-cols-3 gap-6">
        {/* Left Column: Parameter Editing */}
        <div className="space-y-4">
          <Card>
            <CardHeader>
              <CardTitle>设计参数编辑</CardTitle>
              <CardDescription>修改参数值，几何图形将实时更新</CardDescription>
            </CardHeader>
            <CardContent>
              {isLoadingParams ? (
                <div className="text-center py-4">
                  <Loader2 className="h-5 w-5 animate-spin mx-auto mb-2" />
                  <span className="text-sm text-muted-foreground">加载参数信息...</span>
                </div>
              ) : editableParamsList.length === 0 ? (
                <div className="text-center py-4 space-y-2">
                  <p className="text-sm text-muted-foreground">没有可编辑的参数</p>
                  {selectedIndividual?.parameters && Object.keys(selectedIndividual.parameters).length > 0 ? (
                    <div className="text-xs text-muted-foreground">
                      <p>个体包含 {Object.keys(selectedIndividual.parameters).length} 个参数，</p>
                      <p>但没有 free 类型的参数可编辑。</p>
                      <p className="mt-2">几何图形将使用原始参数值生成。</p>
                    </div>
                  ) : (
                    <div className="text-xs text-muted-foreground">
                      <p>个体数据中没有参数信息。</p>
                    </div>
                  )}
                </div>
              ) : (
                <ScrollArea className="h-[600px] pr-4">
                  <div className="space-y-4">
                    {editableParamsList.map(({ name, value, info }) => (
                      <div key={name} className="space-y-2">
                        <Label htmlFor={name} className="text-sm">
                          {info.displayName || name}
                          {info.unit && <span className="text-muted-foreground ml-1">({info.unit})</span>}
                        </Label>
                        <Input
                          id={name}
                          type="number"
                          value={value}
                          onChange={(e) => handleParameterChange(name, e.target.value)}
                          min={info.bounds?.[0]}
                          max={info.bounds?.[1]}
                          step="any"
                          className="w-full"
                        />
                        {info.bounds && (
                          <p className="text-xs text-muted-foreground">
                            范围: [{info.bounds[0].toFixed(4)}, {info.bounds[1].toFixed(4)}]
                          </p>
                        )}
                        {info.comment && (
                          <p className="text-xs text-muted-foreground">{info.comment}</p>
                        )}
                      </div>
                    ))}
                  </div>
                </ScrollArea>
              )}
            </CardContent>
          </Card>
        </div>

        {/* Middle Column: Geometry PDF */}
        <div className="space-y-4">
          <Card>
            <CardHeader>
              <CardTitle>几何图形</CardTitle>
              <CardDescription>实时显示更新后的几何尺寸</CardDescription>
            </CardHeader>
            <CardContent>
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
                    key={geometryPdfPath} // 强制重新加载
                  />
                ) : (
                  <div className="text-center text-muted-foreground">
                    <p className="text-sm">几何图形未生成</p>
                  </div>
                )}
              </div>
            </CardContent>
          </Card>
        </div>

        {/* Right Column: Performance Metrics and FEA Results */}
        <div className="space-y-4">
          {/* Performance Metrics */}
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
              <ScrollArea className="h-[300px]">
                <div className="space-y-4">
                  {/* 目标函数值 */}
                  {selectedFullData.objectives && (
                    <div>
                      <h4 className="font-semibold mb-2 text-sm">目标函数值</h4>
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
                </div>
              </ScrollArea>
            </CardContent>
          </Card>

          {/* FEA Results */}
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
              <CardContent className="h-[400px]">
                <CsvVisualizer 
                  key={`${csvPath}-${selectedIndividual.key}`} 
                  path2FEACsv={csvPath}
                  selectedFile={selectedCsvFileName}
                  onFileSelect={setSelectedCsvFileName}
                  onCurrentFileChange={setCurrentCsvFilePath}
                />
              </CardContent>
            </Card>
          )}
        </div>
      </div>

      {/* Parameters Dictionary Output */}
      <Card>
        <CardHeader>
          <div className="flex items-center justify-between">
            <div>
              <CardTitle>参数字典</CardTitle>
              <CardDescription>调整后的参数字典，可直接复制粘贴使用</CardDescription>
            </div>
            <Button
              variant="outline"
              size="sm"
              onClick={handleCopyParameters}
              className="flex items-center gap-2"
            >
              {copied ? (
                <>
                  <Check className="h-4 w-4" />
                  已复制
                </>
              ) : (
                <>
                  <Copy className="h-4 w-4" />
                  复制
                </>
              )}
            </Button>
          </div>
        </CardHeader>
        <CardContent>
          <div className="bg-muted rounded-lg p-4 font-mono text-sm overflow-x-auto">
            <pre className="whitespace-pre-wrap break-words">{parametersDictString}</pre>
          </div>
        </CardContent>
      </Card>
    </div>
  );
}
