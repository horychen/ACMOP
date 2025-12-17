"use client"

import { useMemo, useState } from "react";
import { BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer, Cell } from "recharts";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";

interface ParameterHistogramProps {
  individuals: any[]; // Pareto前沿个体列表
  parameterName: string; // 参数名称，如 "stator_tooth_span_angle"
  f3Threshold?: number; // f3阈值，默认5
  onBinSelect?: (selectedKeys: string[]) => void; // 选中分段时的回调，返回该分段内的个体key列表
}

interface BinData {
  binIndex: number;
  binStart: number;
  binEnd: number;
  count: number;
  countLowF3: number; // f3 < threshold 的个体数量
  countHighF3: number; // f3 >= threshold 的个体数量
  individualKeys: string[]; // 该分段内的个体key列表
  individualKeysLowF3: string[]; // f3 < threshold 的个体key列表
  individualKeysHighF3: string[]; // f3 >= threshold 的个体key列表
}

/**
 * 创建固定数量的均匀分段：将参数范围均匀分成10段
 */
function createSmartBins(
  individuals: any[],
  parameterName: string,
  numBins: number = 10
): BinData[] {
  // 提取参数值
  const values: Array<{ value: number; key: string; f3: number }> = individuals
    .map(ind => {
      const paramValue = ind.parameters?.[parameterName];
      if (paramValue === undefined || paramValue === null) return null;
      return {
        value: Number(paramValue),
        key: ind.key,
        f3: ind.objectives?.f3 ?? ind.f3 ?? 0
      };
    })
    .filter((v): v is { value: number; key: string; f3: number } => v !== null);

  if (values.length === 0) return [];

  // 找到参数的最小值和最大值
  const paramValues = values.map(v => v.value);
  const minValue = Math.min(...paramValues);
  const maxValue = Math.max(...paramValues);

  // 如果所有值相同，创建一个包含所有值的分段
  if (minValue === maxValue) {
    const countLowF3 = values.filter(v => v.f3 < 5).length;
    const countHighF3 = values.filter(v => v.f3 >= 5).length;
    return [{
      binIndex: 0,
      binStart: minValue,
      binEnd: maxValue,
      count: values.length,
      countLowF3,
      countHighF3,
      individualKeys: values.map(v => v.key),
      individualKeysLowF3: values.filter(v => v.f3 < 5).map(v => v.key),
      individualKeysHighF3: values.filter(v => v.f3 >= 5).map(v => v.key)
    }];
  }

  // 计算每个分段的宽度
  const binWidth = (maxValue - minValue) / numBins;

  // 创建分段
  const bins: BinData[] = [];
  for (let i = 0; i < numBins; i++) {
    const binStart = minValue + i * binWidth;
    const binEnd = i === numBins - 1 ? maxValue : minValue + (i + 1) * binWidth; // 最后一个分段包含最大值

    // 找出属于当前分段的个体
    // 注意：最后一个分段包含右边界，其他分段不包含右边界
    const binValues = values.filter(v => {
      if (i === numBins - 1) {
        // 最后一个分段：包含右边界
        return v.value >= binStart && v.value <= binEnd;
      } else {
        // 其他分段：不包含右边界
        return v.value >= binStart && v.value < binEnd;
      }
    });

    const countLowF3 = binValues.filter(v => v.f3 < 5).length;
    const countHighF3 = binValues.filter(v => v.f3 >= 5).length;

    bins.push({
      binIndex: i,
      binStart: binStart,
      binEnd: binEnd,
      count: binValues.length,
      countLowF3,
      countHighF3,
      individualKeys: binValues.map(v => v.key),
      individualKeysLowF3: binValues.filter(v => v.f3 < 5).map(v => v.key),
      individualKeysHighF3: binValues.filter(v => v.f3 >= 5).map(v => v.key)
    });
  }

  return bins;
}

export function ParameterHistogram({
  individuals,
  parameterName,
  f3Threshold = 5,
  onBinSelect
}: ParameterHistogramProps) {
  const [selectedBinIndex, setSelectedBinIndex] = useState<number | null>(null);

  // 创建均匀分段（10段）
  const bins = useMemo(() => {
    return createSmartBins(individuals, parameterName, 10);
  }, [individuals, parameterName]);

  // 准备图表数据（需要为每个分段创建两个数据点：一个用于f3<5，一个用于f3>=5）
  const chartData = useMemo(() => {
    return bins.map(bin => ({
      binIndex: bin.binIndex,
      binLabel: `${bin.binStart.toFixed(2)}-${bin.binEnd.toFixed(2)}`,
      binStart: bin.binStart,
      binEnd: bin.binEnd,
      countLowF3: bin.countLowF3,
      countHighF3: bin.countHighF3,
      totalCount: bin.count,
      individualKeys: bin.individualKeys,
      individualKeysLowF3: bin.individualKeysLowF3,
      individualKeysHighF3: bin.individualKeysHighF3
    }));
  }, [bins]);

  // 处理分段点击
  const handleBinClick = (data: any) => {
    if (!data) return;
    
    const binIndex = data.binIndex;
    if (selectedBinIndex === binIndex) {
      // 取消选中
      setSelectedBinIndex(null);
      onBinSelect?.([]);
    } else {
      // 选中新分段
      setSelectedBinIndex(binIndex);
      const bin = bins[binIndex];
      onBinSelect?.(bin.individualKeys);
    }
  };

  // 获取选中分段的个体列表
  const selectedBin = selectedBinIndex !== null ? bins[selectedBinIndex] : null;

  if (individuals.length === 0) {
    return (
      <Card>
        <CardHeader>
          <CardTitle>参数分布直方图</CardTitle>
          <CardDescription>暂无数据</CardDescription>
        </CardHeader>
      </Card>
    );
  }

  return (
    <Card>
      <CardHeader>
        <CardTitle>参数分布直方图</CardTitle>
        <CardDescription>
          {parameterName} 的分布情况（f3 &lt; {f3Threshold}: 番茄红，f3 ≥ {f3Threshold}: 灰色）
        </CardDescription>
      </CardHeader>
      <CardContent>
        <div className="space-y-4">
          {/* 直方图 */}
          <div style={{ height: '400px' }}>
            <ResponsiveContainer width="100%" height="100%">
              <BarChart
                data={chartData}
                margin={{ top: 20, right: 30, left: 20, bottom: 60 }}
              >
                <CartesianGrid strokeDasharray="3 3" stroke="rgba(0, 0, 0, 0.1)" />
                <XAxis
                  dataKey="binLabel"
                  angle={-45}
                  textAnchor="end"
                  height={80}
                  tick={{ fontSize: 10 }}
                />
                <YAxis
                  label={{ value: 'Count [1]', angle: -90, position: 'insideLeft' }}
                  tick={{ fontSize: 12 }}
                />
                <Tooltip
                  content={({ active, payload }) => {
                    if (!active || !payload || payload.length === 0) return null;
                    const data = payload[0].payload;
                    return (
                      <div className="bg-background border border-border rounded-lg p-3 shadow-lg">
                        <p className="font-semibold mb-2">{data.binLabel}</p>
                        <p className="text-sm">总个体数: {data.totalCount}</p>
                        <p className="text-sm" style={{ color: '#ff6347' }}>
                          f3 &lt; {f3Threshold}: {data.countLowF3}
                        </p>
                        <p className="text-sm" style={{ color: '#808080' }}>
                          f3 ≥ {f3Threshold}: {data.countHighF3}
                        </p>
                      </div>
                    );
                  }}
                />
                <Legend />
                {/* f3 < threshold 的柱状图（番茄红，带透明度） */}
                <Bar
                  dataKey="countLowF3"
                  name={`f3 < ${f3Threshold}`}
                  fill="#ff6347"
                  opacity={0.7}
                  style={{ cursor: 'pointer' }}
                >
                  {chartData.map((entry, index) => (
                    <Cell
                      key={`cell-low-${index}`}
                      fill={selectedBinIndex === entry.binIndex ? "#ff4500" : "#ff6347"}
                      opacity={selectedBinIndex === entry.binIndex ? 0.9 : 0.7}
                      onClick={() => handleBinClick(entry)}
                    />
                  ))}
                </Bar>
                {/* f3 >= threshold 的柱状图（灰色，带透明度） */}
                <Bar
                  dataKey="countHighF3"
                  name={`f3 ≥ ${f3Threshold}`}
                  fill="#808080"
                  opacity={0.7}
                  style={{ cursor: 'pointer' }}
                >
                  {chartData.map((entry, index) => (
                    <Cell
                      key={`cell-high-${index}`}
                      fill={selectedBinIndex === entry.binIndex ? "#696969" : "#808080"}
                      opacity={selectedBinIndex === entry.binIndex ? 0.9 : 0.7}
                      onClick={() => handleBinClick(entry)}
                    />
                  ))}
                </Bar>
              </BarChart>
            </ResponsiveContainer>
          </div>

          {/* 选中分段的个体列表 */}
          {selectedBin && (
            <div className="border border-border rounded-lg p-4 bg-muted/30">
              <h4 className="font-semibold mb-2">
                选中分段: {selectedBin.binStart.toFixed(2)} - {selectedBin.binEnd.toFixed(2)}
              </h4>
              <div className="space-y-2">
                <div>
                  <p className="text-sm text-muted-foreground mb-1">
                    总个体数: {selectedBin.count} 个
                  </p>
                  <div className="flex flex-wrap gap-2">
                    {selectedBin.individualKeys.map((key) => {
                      const individual = individuals.find(ind => ind.key === key);
                      const f3 = individual?.objectives?.f3 ?? individual?.f3 ?? 0;
                      const isLowF3 = f3 < f3Threshold;
                      return (
                        <Badge
                          key={key}
                          variant="outline"
                          style={{
                            backgroundColor: isLowF3 ? 'rgba(255, 99, 71, 0.2)' : 'rgba(128, 128, 128, 0.2)',
                            borderColor: isLowF3 ? '#ff6347' : '#808080',
                            color: isLowF3 ? '#ff6347' : '#808080'
                          }}
                        >
                          {key} (f3: {f3.toFixed(2)})
                        </Badge>
                      );
                    })}
                  </div>
                </div>
              </div>
            </div>
          )}

          {/* Pareto前沿个体数量统计 */}
          <div className="border-t border-border pt-4">
            <p className="text-sm text-muted-foreground">
              Pareto前沿个体总数: <span className="font-semibold text-foreground">{individuals.length}</span> 个
            </p>
          </div>
        </div>
      </CardContent>
    </Card>
  );
}

