"use client"

import { useState, useEffect } from "react";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { ScrollArea } from "@/components/ui/scroll-area";
import { Loader2, Upload, FileText } from "lucide-react";
import { SensitivityAnalysisViewer } from "./SensitivityAnalysisViewer";
import { OptimizationResultsViewer } from "./OptimizationResultsViewer";

export default function VisualizationMode() {
  const [jsonFile, setJsonFile] = useState<File | null>(null);
  const [jsonContent, setJsonContent] = useState<any>(null);
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [analysisData, setAnalysisData] = useState<any>(null);
  const [optimizationData, setOptimizationData] = useState<any>(null);
  const [activeView, setActiveView] = useState<"sensitivity" | "optimization">("sensitivity");

  const handleFileUpload = (event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0];
    if (file && file.type === "application/json") {
      setJsonFile(file);
      const reader = new FileReader();
      reader.onload = (e) => {
        try {
          const content = JSON.parse(e.target?.result as string);
          setJsonContent(content);
          setError(null);
        } catch (err) {
          setError("无效的JSON文件");
        }
      };
      reader.readAsText(file);
    } else {
      setError("请选择有效的JSON文件");
    }
  };

  const handleLoadAnalysis = async () => {
    if (!jsonContent) {
      setError("请先上传JSON文件");
      return;
    }

    setIsLoading(true);
    setError(null);

    try {
      // 加载敏感性分析数据
      const sensitivityResponse = await fetch('/api/acmopv2/sensitivity-analysis', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({ config: jsonContent }),
      });

      if (!sensitivityResponse.ok) {
        const errorData = await sensitivityResponse.json();
        throw new Error(errorData.error || '加载敏感性分析数据时出错');
      }

      const sensitivityData = await sensitivityResponse.json();
      setAnalysisData(sensitivityData);

      // 加载优化结果数据
      const optimizationResponse = await fetch('/api/acmopv2/optimization-results', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({ config: jsonContent }),
      });

      if (!optimizationResponse.ok) {
        const errorData = await optimizationResponse.json();
        throw new Error(errorData.error || '加载优化结果数据时出错');
      }

      const optimizationData = await optimizationResponse.json();
      setOptimizationData(optimizationData);
    } catch (err: any) {
      setError(err.message || "加载分析数据时出错");
    } finally {
      setIsLoading(false);
    }
  };

  return (
    <div className="space-y-6">
      <Card>
        <CardHeader>
          <CardTitle>加载分析数据</CardTitle>
          <CardDescription>
            上传JSON配置文件以加载敏感性分析和优化结果
          </CardDescription>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="space-y-2">
            <Label htmlFor="json-upload">JSON配置文件</Label>
            <div className="flex items-center gap-4">
              <Input
                id="json-upload"
                type="file"
                accept=".json"
                onChange={handleFileUpload}
                className="flex-1"
              />
              <Button 
                onClick={handleLoadAnalysis} 
                disabled={!jsonContent || isLoading}
                className="min-w-[120px]"
              >
                {isLoading ? (
                  <>
                    <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                    加载中...
                  </>
                ) : (
                  <>
                    <FileText className="mr-2 h-4 w-4" />
                    加载分析
                  </>
                )}
              </Button>
            </div>
          </div>

          {jsonFile && (
            <div className="p-3 bg-muted rounded-md text-sm">
              <div className="flex items-center gap-2">
                <FileText className="h-4 w-4" />
                <span className="font-medium">{jsonFile.name}</span>
              </div>
            </div>
          )}

          {error && (
            <div className="p-4 bg-destructive/10 border border-destructive/20 rounded-md text-destructive text-sm">
              {error}
            </div>
          )}
        </CardContent>
      </Card>

      {(analysisData || optimizationData) && (
        <Card>
          <CardHeader>
            <div className="flex items-center justify-between">
              <div>
                <CardTitle>分析结果可视化</CardTitle>
                <CardDescription>
                  查看敏感性分析和优化结果
                </CardDescription>
              </div>
              <div className="flex gap-2">
                <Button
                  variant={activeView === "sensitivity" ? "default" : "outline"}
                  size="sm"
                  onClick={() => setActiveView("sensitivity")}
                >
                  敏感性分析
                </Button>
                <Button
                  variant={activeView === "optimization" ? "default" : "outline"}
                  size="sm"
                  onClick={() => setActiveView("optimization")}
                >
                  优化结果
                </Button>
              </div>
            </div>
          </CardHeader>
          <CardContent>
            {activeView === "sensitivity" && analysisData && (
              <SensitivityAnalysisViewer data={analysisData} />
            )}
            {activeView === "optimization" && optimizationData && (
              <OptimizationResultsViewer data={optimizationData} />
            )}
          </CardContent>
        </Card>
      )}
    </div>
  );
}

