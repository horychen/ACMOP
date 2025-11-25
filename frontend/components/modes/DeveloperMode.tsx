"use client"

import { useState, useEffect } from "react";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Separator } from "@/components/ui/separator";
import { ScrollArea } from "@/components/ui/scroll-area";
import { Loader2, Download, Play, RefreshCw } from "lucide-react";
import { ParameterConfig } from "./ParameterEditor";
import { CompactParameterTable } from "./CompactParameterTable";

export default function DeveloperMode() {
  const [parameters, setParameters] = useState<ParameterConfig[]>([]);
  const [isLoadingDefaults, setIsLoadingDefaults] = useState(true);
  const [isGenerating, setIsGenerating] = useState(false);
  const [generatedJson, setGeneratedJson] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);

  // 加载默认参数
  useEffect(() => {
    loadDefaultParameters();
  }, []);

  const loadDefaultParameters = async () => {
    setIsLoadingDefaults(true);
    setError(null);
    try {
      const response = await fetch('/api/acmopv2/all-parameters');
      if (!response.ok) {
        throw new Error('加载默认参数失败');
      }
      const data = await response.json();
      
      // 转换后端数据为前端格式
      const allParams: ParameterConfig[] = [];
      
      // Fixed parameters
      data.parameters.fixed?.forEach((p: any) => {
        allParams.push({
          name: p.name,
          type: "fixed",
          value: p.value,
          unit: p.unit || "mm",
          comment: p.comment || "",
          sensitivityAnalysis: false,
        });
      });
      
      // Free parameters
      data.parameters.free?.forEach((p: any) => {
        allParams.push({
          name: p.name,
          type: "free",
          value: p.value,
          bounds: p.bounds || [0, 100],
          calc_bounds: p.calc_bounds || null,
          args: p.args || [],
          unit: p.unit || "mm",
          comment: p.comment || "",
          sensitivityAnalysis: false,
        });
      });
      
      // Derived parameters
      data.parameters.derived?.forEach((p: any) => {
        allParams.push({
          name: p.name,
          type: "derived",
          value: p.value,
          calc: p.calc || `lambda ${p.name}: ${p.name} * 1.0`,
          args: p.args || [],
          unit: p.unit || "mm",
          comment: p.comment || "",
          sensitivityAnalysis: false,
        });
      });
      
      setParameters(allParams);
    } catch (err: any) {
      setError(err.message || "加载默认参数时出错");
    } finally {
      setIsLoadingDefaults(false);
    }
  };

  const updateParameter = (index: number, param: ParameterConfig) => {
    const updated = [...parameters];
    updated[index] = param;
    setParameters(updated);
  };

  // 获取所有参数名称，用于 derived 类型的 calc 函数
  const getAvailableParameterNames = (currentIndex: number): string[] => {
    return parameters
      .filter((_, i) => i !== currentIndex && parameters[i].name)
      .map((p) => p.name);
  };

  const validateParameters = (): string[] => {
    const errors: string[] = [];
    
    parameters.forEach((param, index) => {
      if (!param.name || param.name.trim() === "") {
        errors.push(`参数 ${index + 1}: 参数名称不能为空`);
      }
      
      if (param.type === "free") {
        if (!param.bounds || param.bounds.length !== 2) {
          errors.push(`参数 ${param.name}: Free 类型必须指定边界 (bounds)`);
        } else if (param.bounds[0] >= param.bounds[1]) {
          errors.push(`参数 ${param.name}: 下界必须小于上界`);
        }
      }
      
      if (param.type === "derived") {
        if (!param.calc || param.calc.trim() === "") {
          errors.push(`参数 ${param.name}: Derived 类型必须指定计算函数 (calc)`);
        }
      }
      
      if (param.type === "fixed" && (param.value === null || param.value === undefined)) {
        errors.push(`参数 ${param.name}: Fixed 类型必须指定值`);
      }
    });
    
    // 检查重复的参数名称
    const names = parameters.map(p => p.name).filter(n => n);
    const duplicates = names.filter((name, index) => names.indexOf(name) !== index);
    if (duplicates.length > 0) {
      errors.push(`存在重复的参数名称: ${[...new Set(duplicates)].join(", ")}`);
    }
    
    return errors;
  };

  const handleGenerate = async () => {
    setError(null);
    setGeneratedJson(null);

    // 验证参数
    const validationErrors = validateParameters();
    if (validationErrors.length > 0) {
      setError(validationErrors.join("\n"));
      return;
    }

    setIsGenerating(true);

    try {
      const response = await fetch('/api/acmopv2/generate', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({ 
          parameters: parameters.map(p => ({
            name: p.name,
            type: p.type,
            value: p.value,
            bounds: p.bounds,
            calc: p.calc,
            unit: p.unit || "mm",
            comment: p.comment,
            sensitivityAnalysis: p.sensitivityAnalysis
          }))
        }),
      });

      if (!response.ok) {
        const errorData = await response.json();
        throw new Error(errorData.error || errorData.detail || '生成JSON时出错');
      }

      const data = await response.json();
      setGeneratedJson(JSON.stringify(data, null, 2));
    } catch (err: any) {
      setError(err.message || "生成JSON时出错");
    } finally {
      setIsGenerating(false);
    }
  };

  const handleDownload = () => {
    if (!generatedJson) return;
    
    const blob = new Blob([generatedJson], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `design-config-${Date.now()}.json`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  };

  // 获取需要敏感性分析的参数
  const sensitivityParams = parameters.filter(p => p.sensitivityAnalysis && p.type === "free");

  return (
    <div className="space-y-6">
      <Card>
        <CardHeader>
          <div className="flex items-center justify-between">
            <div>
              <CardTitle>Modern_Machine_Designer 参数配置</CardTitle>
              <CardDescription>
                紧凑的代码编辑器风格界面，所有默认值已预填充
              </CardDescription>
            </div>
            <Button 
              onClick={loadDefaultParameters} 
              variant="outline" 
              size="sm"
              disabled={isLoadingDefaults}
            >
              <RefreshCw className={`mr-2 h-4 w-4 ${isLoadingDefaults ? 'animate-spin' : ''}`} />
              重新加载默认值
            </Button>
          </div>
        </CardHeader>
        <CardContent className="space-y-4">
          {isLoadingDefaults ? (
            <div className="flex items-center justify-center py-12">
              <Loader2 className="h-6 w-6 animate-spin text-muted-foreground" />
              <span className="ml-2 text-sm text-muted-foreground">正在加载默认参数...</span>
            </div>
          ) : parameters.length === 0 ? (
            <div className="text-center py-8 text-muted-foreground">
              <p>未找到默认参数，请检查后端连接</p>
            </div>
          ) : (
            <ScrollArea className="max-h-[600px] pr-4">
              <div className="bg-slate-50 dark:bg-slate-900 rounded-lg p-4 font-mono text-xs">
                <div className="mb-2 text-muted-foreground">
                  <span className="text-blue-600 dark:text-blue-400">@dataclass</span>
                  <br />
                  <span className="text-purple-600 dark:text-purple-400">class</span>{" "}
                  <span className="text-green-600 dark:text-green-400">Modern_Machine_Designer</span>
                  <span className="text-muted-foreground">(object):</span>
                </div>
                <CompactParameterTable
                  parameters={parameters}
                  onParameterChange={updateParameter}
                  availableParameters={getAvailableParameterNames}
                />
              </div>
            </ScrollArea>
          )}

          <Separator />

          {/* 敏感性分析摘要 */}
          {sensitivityParams.length > 0 && (
            <div className="p-4 bg-blue-50 dark:bg-blue-950/30 rounded-md border border-blue-200 dark:border-blue-800">
              <p className="text-sm font-medium text-blue-900 dark:text-blue-100 mb-2">
                敏感性分析参数 ({sensitivityParams.length} 个)
              </p>
              <p className="text-xs text-blue-700 dark:text-blue-300">
                {sensitivityParams.map(p => p.name).join(", ")}
              </p>
            </div>
          )}

          <div className="flex items-center justify-between">
            <div className="text-sm text-muted-foreground">
              {parameters.length > 0 
                ? `已配置 ${parameters.length} 个参数 (Fixed: ${parameters.filter(p => p.type === "fixed").length}, Free: ${parameters.filter(p => p.type === "free").length}, Derived: ${parameters.filter(p => p.type === "derived").length})` 
                : "配置完成后，点击生成按钮创建JSON配置文件"}
            </div>
            <Button 
              onClick={handleGenerate} 
              disabled={isGenerating || parameters.length === 0}
              className="min-w-[120px]"
            >
              {isGenerating ? (
                <>
                  <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                  生成中...
                </>
              ) : (
                <>
                  <Play className="mr-2 h-4 w-4" />
                  生成JSON
                </>
              )}
            </Button>
          </div>

          {error && (
            <div className="p-4 bg-destructive/10 border border-destructive/20 rounded-md">
              <p className="text-destructive text-sm font-medium mb-1">错误:</p>
              <pre className="text-xs text-destructive whitespace-pre-wrap">{error}</pre>
            </div>
          )}
        </CardContent>
      </Card>

      {generatedJson && (
        <Card>
          <CardHeader>
            <div className="flex items-center justify-between">
              <div>
                <CardTitle>生成的JSON配置</CardTitle>
                <CardDescription>
                  生成的配置文件可用于后续的敏感性分析和优化
                </CardDescription>
              </div>
              <Button onClick={handleDownload} variant="outline" size="sm">
                <Download className="mr-2 h-4 w-4" />
                下载JSON
              </Button>
            </div>
          </CardHeader>
          <CardContent>
            <ScrollArea className="h-[500px] w-full rounded-md border p-4">
              <pre className="text-xs font-mono">
                {generatedJson}
              </pre>
            </ScrollArea>
          </CardContent>
        </Card>
      )}
    </div>
  );
}
