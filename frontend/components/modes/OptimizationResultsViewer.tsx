"use client"

import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { ScatterChart, Scatter, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer, ZAxis } from "recharts";

interface OptimizationResultsViewerProps {
  data: any;
}

export function OptimizationResultsViewer({ data }: OptimizationResultsViewerProps) {
  // 这里可以根据实际的数据结构来渲染优化结果
  // 示例：Pareto前沿、优化历史等

  const paretoData = data?.paretoFront || [];
  const optimizationHistory = data?.history || [];
  const bestSolutions = data?.bestSolutions || [];

  return (
    <div className="space-y-6">
      {/* Pareto前沿图 */}
      {paretoData.length > 0 && (
        <Card>
          <CardHeader>
            <CardTitle>Pareto前沿</CardTitle>
            <CardDescription>多目标优化的Pareto最优解集</CardDescription>
          </CardHeader>
          <CardContent>
            <ResponsiveContainer width="100%" height={400}>
              <ScatterChart data={paretoData}>
                <CartesianGrid strokeDasharray="3 3" />
                <XAxis 
                  dataKey="objective1" 
                  name={data?.objectives?.[0] || "目标1"}
                  type="number"
                />
                <YAxis 
                  dataKey="objective2" 
                  name={data?.objectives?.[1] || "目标2"}
                  type="number"
                />
                <ZAxis 
                  dataKey="objective3" 
                  name={data?.objectives?.[2] || "目标3"}
                  type="number"
                  range={[50, 400]}
                />
                <Tooltip cursor={{ strokeDasharray: '3 3' }} />
                <Legend />
                <Scatter 
                  name="Pareto解" 
                  data={paretoData} 
                  fill="#8884d8"
                />
              </ScatterChart>
            </ResponsiveContainer>
          </CardContent>
        </Card>
      )}

      {/* 优化历史 */}
      {optimizationHistory.length > 0 && (
        <Card>
          <CardHeader>
            <CardTitle>优化历史</CardTitle>
            <CardDescription>优化过程中的目标函数值变化</CardDescription>
          </CardHeader>
          <CardContent>
            <ResponsiveContainer width="100%" height={300}>
              <ScatterChart data={optimizationHistory}>
                <CartesianGrid strokeDasharray="3 3" />
                <XAxis dataKey="iteration" name="迭代次数" type="number" />
                <YAxis dataKey="bestValue" name="最佳值" type="number" />
                <Tooltip cursor={{ strokeDasharray: '3 3' }} />
                <Legend />
                <Scatter 
                  name="最佳值" 
                  data={optimizationHistory} 
                  fill="#82ca9d"
                />
              </ScatterChart>
            </ResponsiveContainer>
          </CardContent>
        </Card>
      )}

      {/* 最佳解列表 */}
      {bestSolutions.length > 0 && (
        <Card>
          <CardHeader>
            <CardTitle>最佳解决方案</CardTitle>
            <CardDescription>优化得到的最佳设计参数和性能指标</CardDescription>
          </CardHeader>
          <CardContent>
            <div className="space-y-4">
              {bestSolutions.map((solution: any, index: number) => (
                <Card key={index} className="border-l-4 border-l-primary">
                  <CardHeader>
                    <CardTitle className="text-lg">方案 {index + 1}</CardTitle>
                  </CardHeader>
                  <CardContent>
                    <div className="grid grid-cols-2 gap-4">
                      <div>
                        <h4 className="font-semibold mb-2">设计参数</h4>
                        <div className="space-y-1 text-sm">
                          {Object.entries(solution.parameters || {}).map(([key, value]: [string, any]) => (
                            <div key={key} className="flex justify-between">
                              <span className="text-muted-foreground">{key}:</span>
                              <span className="font-medium">{value}</span>
                            </div>
                          ))}
                        </div>
                      </div>
                      <div>
                        <h4 className="font-semibold mb-2">性能指标</h4>
                        <div className="space-y-1 text-sm">
                          {Object.entries(solution.performance || {}).map(([key, value]: [string, any]) => (
                            <div key={key} className="flex justify-between">
                              <span className="text-muted-foreground">{key}:</span>
                              <span className="font-medium">{value}</span>
                            </div>
                          ))}
                        </div>
                      </div>
                    </div>
                  </CardContent>
                </Card>
              ))}
            </div>
          </CardContent>
        </Card>
      )}

      {paretoData.length === 0 && optimizationHistory.length === 0 && bestSolutions.length === 0 && (
        <Card>
          <CardContent className="py-8 text-center text-muted-foreground">
            暂无优化结果数据
          </CardContent>
        </Card>
      )}
    </div>
  );
}

