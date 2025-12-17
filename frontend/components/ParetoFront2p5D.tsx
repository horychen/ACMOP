"use client"

import { useMemo, useState, useEffect } from "react";
import { ScatterChart, Scatter, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer } from "recharts";
import { fastNonDominatedSorting, extractObjectives, filterByF3, Point } from "@/utils/paretoUtils";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Button } from "@/components/ui/button";
import { Maximize2 } from "lucide-react";
import { Dialog, DialogContent, DialogHeader, DialogTitle, DialogDescription } from "@/components/ui/dialog";

interface ParetoFront2p5DProps {
  individuals: Point[];
  objectives?: string[]; // [f1名称, f2名称, f3名称]
  comp?: [number, number]; // [x轴目标索引, y轴目标索引]，默认 [0, 1]，z轴自动为剩余的一个
  upToRankNo?: number; // 显示到第几层前沿，默认 1
  zFilter?: number; // f3 过滤阈值（例如 20）
  onZFilterChange?: (value: number | undefined) => void;
  onIndividualSelect?: (individualKey: string) => void; // 选择个体时的回调
  highlightedKeys?: string[]; // 高亮的个体key列表，未在此列表中的个体将变灰
}

/**
 * 将数值映射到颜色（使用 plasma colormap 类似的效果）
 * plasma colormap: 从深紫色到黄色
 */
function valueToColor(value: number, min: number, max: number): string {
  if (max === min) return "#0d0887"; // 单一值，返回深紫色
  
  // 归一化到 [0, 1]
  const normalized = (value - min) / (max - min);
  
  // Plasma colormap 的近似实现
  // 从深紫色 (#0d0887) 到黄色 (#f0f921)
  const colors = [
    { r: 13, g: 8, b: 135 },    // #0d0887
    { r: 75, g: 10, b: 161 },   // #4b0aa1
    { r: 125, g: 13, b: 135 },  // #7d0d87
    { r: 168, g: 42, b: 99 },   // #a82a63
    { r: 203, g: 71, b: 119 },  // #cb4777
    { r: 229, g: 107, b: 93 },  // #e56b5d
    { r: 248, g: 148, b: 65 },  // #f89441
    { r: 253, g: 195, b: 40 },  // #fdc328
    { r: 240, g: 249, b: 33 }   // #f0f921
  ];
  
  const index = Math.floor(normalized * (colors.length - 1));
  const t = (normalized * (colors.length - 1)) - index;
  
  const c1 = colors[Math.min(index, colors.length - 1)];
  const c2 = colors[Math.min(index + 1, colors.length - 1)];
  
  const r = Math.round(c1.r + t * (c2.r - c1.r));
  const g = Math.round(c1.g + t * (c2.g - c1.g));
  const b = Math.round(c1.b + t * (c2.b - c1.b));
  
  return `rgb(${r}, ${g}, ${b})`;
}

/**
 * 计算合理的坐标轴刻度值
 * @param min 最小值（数据范围的最小值，已包含 padding）
 * @param max 最大值（数据范围的最大值，已包含 padding）
 * @returns 刻度值数组，确保第一个值 <= min，最后一个值 >= max
 */
function calculateTicks(min: number, max: number): number[] {
  const range = max - min;
  
  // 处理范围很小或为 0 的情况
  if (range === 0 || range < 0.01) {
    const center = (min + max) / 2;
    const padding = Math.max(Math.abs(center) * 0.1, 0.5);
    return [
      Math.floor((center - padding) * 10) / 10,
      Math.floor((center - padding / 2) * 10) / 10,
      Math.round(center * 10) / 10,
      Math.ceil((center + padding / 2) * 10) / 10,
      Math.ceil((center + padding) * 10) / 10
    ];
  }
  
  // 计算合适的间隔，使得总共有 5-7 个刻度
  let bestInterval = range / 6; // 默认 6 个间隔（7 个刻度）
  
  // 将间隔调整为"友好"的数字（1, 2, 5, 10, 20, 50, 100, 200, 500, 1000 等）
  const magnitude = Math.pow(10, Math.floor(Math.log10(bestInterval)));
  const normalized = bestInterval / magnitude;
  
  let friendlyInterval: number;
  if (normalized <= 1) {
    friendlyInterval = 1 * magnitude;
  } else if (normalized <= 2) {
    friendlyInterval = 2 * magnitude;
  } else if (normalized <= 5) {
    friendlyInterval = 5 * magnitude;
  } else {
    friendlyInterval = 10 * magnitude;
  }
  
  // 确保刻度范围完全覆盖数据范围
  // 起点应该 <= min，终点应该 >= max
  const actualStart = Math.floor(min / friendlyInterval) * friendlyInterval;
  const actualEnd = Math.ceil(max / friendlyInterval) * friendlyInterval;
  
  // 生成刻度值
  const ticks: number[] = [];
  let current = actualStart;
  while (current <= actualEnd + friendlyInterval * 0.001) { // 添加小的容差避免浮点误差
    ticks.push(current);
    current += friendlyInterval;
    if (ticks.length > 15) break; // 防止无限循环
  }
  
  // 确保刻度范围覆盖数据范围
  if (ticks.length > 0) {
    if (ticks[0] > min) {
      // 如果第一个刻度大于最小值，添加更小的刻度
      let extraTick = ticks[0] - friendlyInterval;
      while (extraTick >= min - friendlyInterval * 0.001) {
        ticks.unshift(extraTick);
        extraTick -= friendlyInterval;
        if (ticks.length > 15) break;
      }
    }
    if (ticks[ticks.length - 1] < max) {
      // 如果最后一个刻度小于最大值，添加更大的刻度
      let extraTick = ticks[ticks.length - 1] + friendlyInterval;
      while (extraTick <= max + friendlyInterval * 0.001) {
        ticks.push(extraTick);
        extraTick += friendlyInterval;
        if (ticks.length > 15) break;
      }
    }
  }
  
  // 确保至少有 5 个刻度，最多 7 个
  if (ticks.length < 5) {
    // 如果刻度太少，减小间隔
    const smallerInterval = friendlyInterval / 2;
    const newStart = Math.floor(min / smallerInterval) * smallerInterval;
    const newEnd = Math.ceil(max / smallerInterval) * smallerInterval;
    const newTicks: number[] = [];
    current = newStart;
    while (current <= newEnd + smallerInterval * 0.001 && newTicks.length < 10) {
      newTicks.push(current);
      current += smallerInterval;
    }
    // 确保覆盖范围
    if (newTicks.length > 0) {
      if (newTicks[0] > min) {
        let extraTick = newTicks[0] - smallerInterval;
        while (extraTick >= min - smallerInterval * 0.001 && newTicks.length < 10) {
          newTicks.unshift(extraTick);
          extraTick -= smallerInterval;
        }
      }
      if (newTicks[newTicks.length - 1] < max) {
        let extraTick = newTicks[newTicks.length - 1] + smallerInterval;
        while (extraTick <= max + smallerInterval * 0.001 && newTicks.length < 10) {
          newTicks.push(extraTick);
          extraTick += smallerInterval;
        }
      }
    }
    return newTicks.length >= 5 ? (newTicks.length <= 7 ? newTicks : newTicks.slice(0, 7)) : newTicks;
  } else if (ticks.length > 7) {
    // 如果刻度太多，尝试使用更大的间隔，但确保覆盖范围
    const largerInterval = friendlyInterval * 2;
    const newStart = Math.floor(min / largerInterval) * largerInterval;
    const newEnd = Math.ceil(max / largerInterval) * largerInterval;
    const newTicks: number[] = [];
    current = newStart;
    while (current <= newEnd + largerInterval * 0.001 && newTicks.length < 10) {
      newTicks.push(current);
      current += largerInterval;
    }
    // 确保覆盖范围
    if (newTicks.length > 0) {
      if (newTicks[0] > min) {
        let extraTick = newTicks[0] - largerInterval;
        while (extraTick >= min - largerInterval * 0.001 && newTicks.length < 10) {
          newTicks.unshift(extraTick);
          extraTick -= largerInterval;
        }
      }
      if (newTicks[newTicks.length - 1] < max) {
        let extraTick = newTicks[newTicks.length - 1] + largerInterval;
        while (extraTick <= max + largerInterval * 0.001 && newTicks.length < 10) {
          newTicks.push(extraTick);
          extraTick += largerInterval;
        }
      }
    }
    if (newTicks.length >= 5 && newTicks.length <= 7) {
      return newTicks;
    }
    // 如果更大的间隔导致刻度太少或太多，返回原来的，但确保覆盖范围
    // 如果原来的刻度太多，只保留必要的部分以确保覆盖
    if (ticks.length > 7) {
      // 保留第一个和最后一个，中间均匀选择 5 个
      const result = [ticks[0]];
      const step = Math.floor((ticks.length - 1) / 6);
      for (let i = 1; i < ticks.length - 1; i += step) {
        if (result.length < 6) {
          result.push(ticks[i]);
        }
      }
      result.push(ticks[ticks.length - 1]);
      return result;
    }
  }
  
  return ticks;
}

/**
 * 生成颜色条渐变
 */
function ColorBar({ min, max, label, height = 200 }: { min: number; max: number; label: string; height?: number | string }) {
  const steps = 100;
  const gradientId = `colorbar-gradient-${Math.random().toString(36).substr(2, 9)}`;
  
  // 将高度转换为数字（如果是字符串，尝试解析）
  const heightNum = typeof height === 'string' ? parseFloat(height) || 200 : height;
  
  const stops = [];
  for (let i = 0; i <= steps; i++) {
    const value = min + (max - min) * (i / steps);
    const color = valueToColor(value, min, max);
    stops.push(
      <stop key={i} offset={`${(i / steps) * 100}%`} stopColor={color} />
    );
  }
  
  // 生成刻度值（更多刻度）
  const numTicks = 8; // 刻度数量
  const tickValues: number[] = [];
  for (let i = 0; i <= numTicks; i++) {
    const value = min + (max - min) * (i / numTicks);
    tickValues.push(value);
  }
  
  // 计算刻度位置（从底部到顶部）
  // 最小值在底部，最大值在顶部
  const tickPositions = tickValues.map((_, i) => (1 - i / numTicks) * heightNum);
  
  return (
    <div className="flex items-center gap-2">
      <div className="flex flex-col items-end">
        <div className="text-xs text-muted-foreground mb-1">{label}</div>
        <div className="flex items-start gap-2">
          {/* 刻度标签（左侧，有足够空间） */}
          <div className="relative text-xs text-muted-foreground text-right" style={{ height: heightNum }}>
            {tickValues.map((value, i) => (
              <span
                key={i}
                className="block"
                style={{
                  position: 'absolute',
                  top: `${tickPositions[i]}px`,
                  transform: 'translateY(-50%)',
                  whiteSpace: 'nowrap',
                  right: '4px',
                  minWidth: '50px'
                }}
              >
                {value.toFixed(2)}
              </span>
            ))}
          </div>
          {/* 颜色条 */}
          <div className="relative flex-shrink-0">
            <svg width="20" height={heightNum} className="border border-border rounded" style={{ height: heightNum }}>
              <defs>
                <linearGradient id={gradientId} x1="0%" y1="100%" x2="0%" y2="0%">
                  {stops}
                </linearGradient>
              </defs>
              <rect width="20" height={heightNum} fill={`url(#${gradientId})`} />
              {/* 绘制刻度线 */}
              {tickPositions.map((y, i) => (
                <line
                  key={i}
                  x1={0}
                  y1={y}
                  x2={20}
                  y2={y}
                  stroke="rgba(0, 0, 0, 0.3)"
                  strokeWidth={0.5}
                />
              ))}
            </svg>
          </div>
        </div>
      </div>
    </div>
  );
}

// 全屏图表组件
function FullscreenChart({
  height,
  chartData,
  coloredData,
  xLabel,
  yLabel,
  zLabel,
  xTicks,
  yTicks,
  onIndividualSelect,
  handlePointClick
}: {
  height: number;
  chartData: any;
  coloredData: any[];
  xLabel: string;
  yLabel: string;
  zLabel: string;
  xTicks: number[];
  yTicks: number[];
  onIndividualSelect?: (key: string) => void;
  handlePointClick: (data: any) => void;
}) {
  // 过滤掉可能重叠的刻度（如果 X 轴和 Y 轴的最小刻度都接近 0，则移除 Y 轴的最小刻度）
  const xMinTick = xTicks[0];
  const yMinTick = yTicks[0];
  const threshold = Math.max(Math.abs(chartData.xMax - chartData.xMin), Math.abs(chartData.yMax - chartData.yMin)) * 0.01;
  
  // 如果两个轴的最小刻度都接近 0，移除 Y 轴的最小刻度以避免重叠
  const filteredXTicks = xTicks;
  const filteredYTicks = (Math.abs(xMinTick) < threshold && Math.abs(yMinTick) < threshold && yTicks.length > 1)
    ? yTicks.slice(1)
    : yTicks;
  
  return (
    <div className="flex gap-4 items-start h-full">
      <div className="flex-1 h-full">
        <ResponsiveContainer width="100%" height={height}>
          <ScatterChart margin={{ top: 20, right: 20, bottom: 20, left: 20 }}>
            <CartesianGrid 
              strokeDasharray="3 3" 
              stroke="rgba(0, 0, 0, 0.08)"
              strokeOpacity={0.3}
            />
            <XAxis
              type="number"
              dataKey="x"
              name={xLabel}
              label={{ value: xLabel, position: "insideBottom", offset: -5 }}
              domain={[chartData.xMin, chartData.xMax]}
              ticks={filteredXTicks}
              tickFormatter={(value) => value.toFixed(1)}
              allowDataOverflow={false}
              tick={{ fontSize: 12 }}
              tickMargin={8}
            />
            <YAxis
              type="number"
              dataKey="y"
              name={yLabel}
              label={{ value: yLabel, angle: -90, position: "insideLeft" }}
              domain={[chartData.yMin, chartData.yMax]}
              ticks={filteredYTicks}
              tickFormatter={(value) => value.toFixed(1)}
              allowDataOverflow={false}
              tick={{ fontSize: 12 }}
              tickMargin={8}
              width={60}
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
                        {xLabel}: {data.x?.toFixed(4)}
                      </p>
                      <p className="text-sm">
                        {yLabel}: {data.y?.toFixed(4)}
                      </p>
                      <p className="text-sm">
                        {zLabel}: {data.z?.toFixed(4)}
                      </p>
                      <p className="text-xs text-muted-foreground mt-1">
                        前沿等级: Rank {data.rank}
                      </p>
                      {onIndividualSelect && (
                        <p className="text-xs text-primary mt-2 font-medium">
                          点击选择此个体
                        </p>
                      )}
                    </div>
                  );
                }
                return null;
              }}
            />
            <Scatter
              name="Pareto解"
              data={coloredData}
              onClick={(data: any) => {
                // Recharts 的 onClick 事件参数格式：{ payload: data, ... }
                const pointData = data?.payload || data;
                if (onIndividualSelect && pointData?.key) {
                  handlePointClick(pointData);
                }
              }}
              shape={(props: any) => {
                const { cx, cy, payload } = props;
                if (!payload || !payload.fill || cx === undefined || cy === undefined) {
                  return <circle cx={0} cy={0} r={0} />;
                }
                return (
                  <circle
                    cx={cx}
                    cy={cy}
                    r={6}
                    fill={payload.fill}
                    stroke="rgba(0,0,0,0.3)"
                    strokeWidth={1}
                    opacity={payload.opacity ?? 0.8}
                    style={{ cursor: onIndividualSelect ? 'pointer' : 'default' }}
                    onClick={(e) => {
                      e.stopPropagation();
                      if (onIndividualSelect && payload?.key) {
                        handlePointClick(payload);
                      }
                    }}
                  />
                );
              }}
            />
          </ScatterChart>
        </ResponsiveContainer>
      </div>
      
      {/* 颜色条 */}
      <ColorBar
        min={chartData.zMin}
        max={chartData.zMax}
        label={zLabel}
        height={height}
      />
    </div>
  );
}

export function ParetoFront2p5D({
  individuals,
  objectives = ["f1", "f2", "f3"],
  comp = [0, 1],
  upToRankNo = 1,
  zFilter,
  onZFilterChange,
  onIndividualSelect,
  highlightedKeys = undefined
}: ParetoFront2p5DProps) {
  const [isMaximized, setIsMaximized] = useState(false);
  const [fullscreenHeight, setFullscreenHeight] = useState(600);

  // 计算全屏高度
  useEffect(() => {
    if (isMaximized && typeof window !== 'undefined') {
      const updateHeight = () => {
        setFullscreenHeight(window.innerHeight * 0.95 - 140);
      };
      updateHeight();
      window.addEventListener('resize', updateHeight);
      return () => window.removeEventListener('resize', updateHeight);
    }
  }, [isMaximized]);

  // 确定 z 轴目标索引（第三个目标）
  const zComp = useMemo(() => {
    const fullComp = [0, 1, 2];
    fullComp.splice(fullComp.indexOf(comp[0]), 1);
    fullComp.splice(fullComp.indexOf(comp[1]), 1);
    return fullComp[0];
  }, [comp]);

  // 应用 z 轴目标过滤（根据 zComp 确定要过滤的目标函数）
  const filteredIndividuals = useMemo(() => {
    if (zFilter === undefined || zFilter === null) {
      return individuals;
    }
    return individuals.filter(ind => {
      // 根据 zComp 获取对应的目标函数值
      const zValue = zComp === 0 ? ind.f1 ?? ind.objectives?.f1 ?? 0
            : zComp === 1 ? ind.f2 ?? ind.objectives?.f2 ?? 0
            : ind.f3 ?? ind.objectives?.f3 ?? 0;
      
      // 应用缩放（如果 zComp 是效率，需要乘以100）
      const zScale = zComp === 1 ? 100 : 1;
      const scaledZValue = zValue * zScale;
      
      return scaledZValue < zFilter;
    });
  }, [individuals, zFilter, zComp]);

  // 提取目标函数值
  const points = useMemo(() => {
    return extractObjectives(filteredIndividuals);
  }, [filteredIndividuals]);

  // 非支配排序
  const fronts = useMemo(() => {
    if (points.length === 0) return [];
    return fastNonDominatedSorting(points);
  }, [points]);

  // 准备图表数据
  const chartData = useMemo(() => {
    const data: any[] = [];
    const allZValues: number[] = [];
    const allXValues: number[] = [];
    const allYValues: number[] = [];

    // 处理每个前沿（最多到 upToRankNo）
    for (let rank = 0; rank < Math.min(fronts.length, upToRankNo); rank++) {
      const front = fronts[rank];
      
      for (const idx of front) {
        if (idx >= filteredIndividuals.length) continue;
        
        const ind = filteredIndividuals[idx];
        const x = comp[0] === 0 ? ind.f1 ?? ind.objectives?.f1 ?? 0
              : comp[0] === 1 ? ind.f2 ?? ind.objectives?.f2 ?? 0
              : ind.f3 ?? ind.objectives?.f3 ?? 0;
        const y = comp[1] === 0 ? ind.f1 ?? ind.objectives?.f1 ?? 0
              : comp[1] === 1 ? ind.f2 ?? ind.objectives?.f2 ?? 0
              : ind.f3 ?? ind.objectives?.f3 ?? 0;
        const z = zComp === 0 ? ind.f1 ?? ind.objectives?.f1 ?? 0
              : zComp === 1 ? ind.f2 ?? ind.objectives?.f2 ?? 0
              : ind.f3 ?? ind.objectives?.f3 ?? 0;

        // 应用缩放（效率需要乘以100）
        const xScale = comp[0] === 1 ? 100 : 1;
        const yScale = comp[1] === 1 ? 100 : 1;
        const zScale = zComp === 1 ? 100 : 1;

        const scaledX = x * xScale;
        const scaledY = y * yScale;
        const scaledZ = z * zScale;

        data.push({
          x: scaledX,
          y: scaledY,
          z: scaledZ,
          name: `Gen${ind.generation ?? '?'}-Ind${ind.individual_index ?? ind.index ?? '?'}`,
          key: ind.key,
          originalIndex: idx,
          rank: rank + 1
        });

        allXValues.push(scaledX);
        allYValues.push(scaledY);
        allZValues.push(scaledZ);
      }
    }

    if (allZValues.length === 0) {
      return { 
        data, 
        zMin: 0, 
        zMax: 1,
        xMin: 0,
        xMax: 1,
        yMin: 0,
        yMax: 1
      };
    }

    // 计算范围并添加 padding（5%）
    const calculateRange = (values: number[]) => {
      const min = Math.min(...values);
      const max = Math.max(...values);
      const range = max - min;
      
      // 如果所有值相同，添加一个小的默认范围
      if (range === 0) {
        const center = min;
        const defaultRange = Math.abs(center) * 0.1 || 1; // 使用中心值的10%或默认1
        return {
          min: center - defaultRange,
          max: center + defaultRange
        };
      }
      
      const padding = range * 0.05; // 5% padding
      return {
        min: min - padding,
        max: max + padding
      };
    };

    const xRange = calculateRange(allXValues);
    const yRange = calculateRange(allYValues);

    return { 
      data, 
      zMin: Math.min(...allZValues), 
      zMax: Math.max(...allZValues),
      xMin: xRange.min,
      xMax: xRange.max,
      yMin: yRange.min,
      yMax: yRange.max
    };
  }, [fronts, filteredIndividuals, comp, zComp, upToRankNo]);

  // 为每个点添加颜色（支持高亮）
  const coloredData = useMemo(() => {
    if (chartData.data.length === 0) return [];
    const { zMin, zMax } = chartData;
    
    return chartData.data.map(point => {
      const baseColor = valueToColor(point.z, zMin, zMax);
      
      // 如果指定了高亮列表，未选中的个体变灰
      if (highlightedKeys !== undefined && highlightedKeys.length > 0) {
        const isHighlighted = highlightedKeys.includes(point.key);
        return {
          ...point,
          fill: isHighlighted ? baseColor : '#d3d3d3', // 未选中的变灰
          opacity: isHighlighted ? 0.8 : 0.3 // 未选中的降低透明度
        };
      }
      
      return {
        ...point,
        fill: baseColor,
        opacity: 0.8
      };
    });
  }, [chartData, highlightedKeys]);

  // 注意：Pareto 前沿轮廓（step plot）在 ScatterChart 中较难实现
  // 如果需要，可以使用叠加的 LineChart 或自定义 SVG 路径

  if (coloredData.length === 0) {
    return (
      <div className="text-center py-8 text-muted-foreground">
        暂无数据可显示
      </div>
    );
  }

  const xLabel = objectives[comp[0]] || `目标${comp[0] + 1}`;
  const yLabel = objectives[comp[1]] || `目标${comp[1] + 1}`;
  const zLabel = objectives[zComp] || `目标${zComp + 1}`;

  // 计算坐标轴刻度
  const xTicks = calculateTicks(chartData.xMin, chartData.xMax);
  const yTicks = calculateTicks(chartData.yMin, chartData.yMax);

  // 处理点击选择个体
  const handlePointClick = (data: any) => {
    if (onIndividualSelect && data?.key) {
      onIndividualSelect(data.key);
      setIsMaximized(false); // 选择后关闭全屏视图
    }
  };

  // 渲染图表的函数（可在普通视图和全屏视图中复用）
  const renderChart = (height: number | string | undefined = 400, isFullscreen: boolean = false) => {
    // 对于全屏模式，确保有明确的高度值（转换为数字）
    const chartHeight = typeof height === 'string' 
      ? (isFullscreen ? parseFloat(height) || 400 : 400)
      : (height || 400);
    
    // 在函数内部计算过滤后的刻度，避免重叠
    const xMinTick = xTicks[0];
    const yMinTick = yTicks[0];
    const threshold = Math.max(Math.abs(chartData.xMax - chartData.xMin), Math.abs(chartData.yMax - chartData.yMin)) * 0.01;
    
    // 如果两个轴的最小刻度都接近 0，移除 Y 轴的最小刻度以避免重叠
    const filteredXTicks = xTicks;
    const filteredYTicks = (Math.abs(xMinTick) < threshold && Math.abs(yMinTick) < threshold && yTicks.length > 1)
      ? yTicks.slice(1)
      : yTicks;
    
    return (
    <div className="flex gap-4 items-start" style={isFullscreen ? { height: chartHeight } : {}}>
      <div className="flex-1" style={isFullscreen ? { height: chartHeight } : {}}>
        <ResponsiveContainer width="100%" height={chartHeight} key={isFullscreen ? "fullscreen" : "normal"}>
          <ScatterChart margin={{ top: 20, right: 20, bottom: 20, left: 20 }}>
            <CartesianGrid 
              strokeDasharray="3 3" 
              stroke="rgba(0, 0, 0, 0.08)"
              strokeOpacity={0.3}
            />
            <XAxis
              type="number"
              dataKey="x"
              name={xLabel}
              label={{ value: xLabel, position: "insideBottom", offset: -5 }}
              domain={[chartData.xMin, chartData.xMax]}
              ticks={filteredXTicks}
              tickFormatter={(value) => value.toFixed(1)}
              allowDataOverflow={false}
              tick={{ fontSize: 12 }}
              tickMargin={8}
            />
            <YAxis
              type="number"
              dataKey="y"
              name={yLabel}
              label={{ value: yLabel, angle: -90, position: "insideLeft" }}
              domain={[chartData.yMin, chartData.yMax]}
              ticks={filteredYTicks}
              tickFormatter={(value) => value.toFixed(1)}
              allowDataOverflow={false}
              tick={{ fontSize: 12 }}
              tickMargin={8}
              width={60}
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
                        {xLabel}: {data.x?.toFixed(4)}
                      </p>
                      <p className="text-sm">
                        {yLabel}: {data.y?.toFixed(4)}
                      </p>
                      <p className="text-sm">
                        {zLabel}: {data.z?.toFixed(4)}
                      </p>
                      <p className="text-xs text-muted-foreground mt-1">
                        前沿等级: Rank {data.rank}
                      </p>
                      {onIndividualSelect && (
                        <p className="text-xs text-primary mt-2 font-medium">
                          点击选择此个体
                        </p>
                      )}
                    </div>
                  );
                }
                return null;
              }}
            />
            <Scatter
              name="Pareto解"
              data={coloredData}
              onClick={(data: any) => {
                // Recharts 的 onClick 事件参数格式：{ payload: data, ... }
                const pointData = data?.payload || data;
                if (onIndividualSelect && pointData?.key) {
                  handlePointClick(pointData);
                }
              }}
              shape={(props: any) => {
                const { cx, cy, payload } = props;
                if (!payload || !payload.fill || cx === undefined || cy === undefined) {
                  return <circle cx={0} cy={0} r={0} />;
                }
                return (
                  <circle
                    cx={cx}
                    cy={cy}
                    r={isFullscreen ? 6 : 5}
                    fill={payload.fill}
                    stroke={isFullscreen ? "rgba(0,0,0,0.3)" : "rgba(0,0,0,0.1)"}
                    strokeWidth={isFullscreen ? 1 : 0.5}
                    opacity={payload.opacity ?? 0.8}
                    style={{ cursor: onIndividualSelect ? 'pointer' : 'default' }}
                    onClick={(e) => {
                      e.stopPropagation();
                      if (onIndividualSelect && payload?.key) {
                        handlePointClick(payload);
                      }
                    }}
                  />
                );
              }}
            />
          </ScatterChart>
        </ResponsiveContainer>
      </div>
      
      {/* 颜色条 */}
      <ColorBar
        min={chartData.zMin}
        max={chartData.zMax}
        label={zLabel}
        height={isFullscreen ? (typeof chartHeight === 'string' ? parseFloat(chartHeight) || 400 : chartHeight) : 200}
      />
    </div>
    );
  };

  return (
    <div className="space-y-4">
      {/* 过滤控制 */}
      {onZFilterChange && (
        <div className="flex items-center gap-4">
          <Label htmlFor="z-filter" className="whitespace-nowrap">
            {zLabel} 过滤阈值:
          </Label>
          <Input
            id="z-filter"
            type="number"
            placeholder="例如: 20"
            value={zFilter ?? ""}
            onChange={(e) => {
              const value = e.target.value === "" ? undefined : parseFloat(e.target.value);
              onZFilterChange(isNaN(value as number) ? undefined : value);
            }}
            className="w-32"
          />
          {zFilter !== undefined && (
            <span className="text-sm text-muted-foreground">
              显示 {zLabel} &lt; {zFilter} 的个体
            </span>
          )}
        </div>
      )}

      {/* 图表区域 */}
      <div className="relative">
        {renderChart(400, false)}
        {/* 最大化按钮 */}
        <Button
          variant="outline"
          size="icon"
          className="absolute top-2 right-2 z-10"
          onClick={() => setIsMaximized(true)}
          title="最大化查看"
        >
          <Maximize2 className="h-4 w-4" />
        </Button>
      </div>

      {/* 全屏Dialog */}
      <Dialog open={isMaximized} onOpenChange={setIsMaximized}>
        <DialogContent className="max-w-none w-[95vw] h-[95vh] p-6 flex flex-col" style={{ maxWidth: '95vw', width: '95vw', height: '95vh' }}>
          <DialogHeader className="flex-shrink-0 mb-4">
            <DialogTitle>Pareto前沿 - 全屏查看</DialogTitle>
            <DialogDescription>
              点击图表中的点来选择个体。点击后会自动关闭此窗口并更新个体选择。
            </DialogDescription>
          </DialogHeader>
          <div className="flex-1 overflow-hidden min-h-0" style={{ height: 'calc(95vh - 140px)' }}>
            <FullscreenChart 
              height={fullscreenHeight}
              chartData={chartData}
              coloredData={coloredData}
              xLabel={xLabel}
              yLabel={yLabel}
              zLabel={zLabel}
              xTicks={xTicks}
              yTicks={yTicks}
              onIndividualSelect={onIndividualSelect}
              handlePointClick={handlePointClick}
            />
          </div>
        </DialogContent>
      </Dialog>
    </div>
  );
}

