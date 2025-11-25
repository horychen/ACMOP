"use client"

import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer } from "recharts";

interface SensitivityAnalysisViewerProps {
  data: any;
}

export function SensitivityAnalysisViewer({ data }: SensitivityAnalysisViewerProps) {
  // 这里可以根据实际的数据结构来渲染敏感性分析结果
  // 示例：假设数据包含参数变化对性能指标的影响

  const chartData = data?.sensitivityData || [];
  const parameters = data?.parameters || [];

  return (
    <div className="space-y-6">
      <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
        {parameters.map((param: string, index: number) => {
          const paramData = chartData.filter((d: any) => d.parameter === param);
          
          return (
            <Card key={index}>
              <CardHeader>
                <CardTitle className="text-lg">{param}</CardTitle>
                <CardDescription>参数敏感性分析</CardDescription>
              </CardHeader>
              <CardContent>
                <ResponsiveContainer width="100%" height={300}>
                  <LineChart data={paramData}>
                    <CartesianGrid strokeDasharray="3 3" />
                    <XAxis dataKey="value" />
                    <YAxis />
                    <Tooltip />
                    <Legend />
                    <Line 
                      type="monotone" 
                      dataKey="efficiency" 
                      stroke="#8884d8" 
                      name="效率 (%)"
                    />
                    <Line 
                      type="monotone" 
                      dataKey="torque" 
                      stroke="#82ca9d" 
                      name="转矩 (Nm)"
                    />
                  </LineChart>
                </ResponsiveContainer>
              </CardContent>
            </Card>
          );
        })}
      </div>

      {data?.summary && (
        <Card>
          <CardHeader>
            <CardTitle>敏感性分析摘要</CardTitle>
          </CardHeader>
          <CardContent>
            <div className="space-y-2">
              {Object.entries(data.summary).map(([key, value]: [string, any]) => (
                <div key={key} className="flex justify-between items-center p-2 bg-muted rounded">
                  <span className="font-medium">{key}</span>
                  <span className="text-sm text-muted-foreground">{JSON.stringify(value)}</span>
                </div>
              ))}
            </div>
          </CardContent>
        </Card>
      )}
    </div>
  );
}

