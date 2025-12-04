"use client"

import { useMemo } from "react";
import { ScatterChart, Scatter, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer } from "recharts";
import { fastNonDominatedSorting, extractObjectives, filterByF3, Point } from "@/utils/paretoUtils";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";

interface ParetoFront2p5DProps {
  individuals: Point[];
  objectives?: string[]; // [f1名称, f2名称, f3名称]
  comp?: [number, number]; // [x轴目标索引, y轴目标索引]，默认 [0, 1]，z轴自动为剩余的一个
  upToRankNo?: number; // 显示到第几层前沿，默认 1
  zFilter?: number; // f3 过滤阈值（例如 20）
  onZFilterChange?: (value: number | undefined) => void;
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
 * 生成颜色条渐变
 */
function ColorBar({ min, max, label }: { min: number; max: number; label: string }) {
  const steps = 100;
  const gradientId = `colorbar-gradient-${Math.random().toString(36).substr(2, 9)}`;
  
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
  const tickPositions = tickValues.map((_, i) => (i / numTicks) * 200);
  
  return (
    <div className="flex items-center gap-2">
      <div className="flex flex-col items-end">
        <div className="text-xs text-muted-foreground mb-1">{label}</div>
        <div className="flex items-start gap-2">
          {/* 刻度标签（左侧，有足够空间） */}
          <div className="relative h-[200px] text-xs text-muted-foreground text-right">
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
            <svg width="20" height="200" className="border border-border rounded">
              <defs>
                <linearGradient id={gradientId} x1="0%" y1="100%" x2="0%" y2="0%">
                  {stops}
                </linearGradient>
              </defs>
              <rect width="20" height="200" fill={`url(#${gradientId})`} />
              {/* 绘制刻度线 */}
              {tickPositions.map((y, i) => (
                <line
                  key={i}
                  x1={0}
                  y1={200 - y}
                  x2={20}
                  y2={200 - y}
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

export function ParetoFront2p5D({
  individuals,
  objectives = ["f1", "f2", "f3"],
  comp = [0, 1],
  upToRankNo = 1,
  zFilter,
  onZFilterChange
}: ParetoFront2p5DProps) {
  // 确定 z 轴目标索引（第三个目标）
  const zComp = useMemo(() => {
    const fullComp = [0, 1, 2];
    fullComp.splice(fullComp.indexOf(comp[0]), 1);
    fullComp.splice(fullComp.indexOf(comp[1]), 1);
    return fullComp[0];
  }, [comp]);

  // 应用 f3 过滤
  const filteredIndividuals = useMemo(() => {
    if (zFilter === undefined || zFilter === null) {
      return individuals;
    }
    return filterByF3(individuals, zFilter);
  }, [individuals, zFilter]);

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

        data.push({
          x: x * xScale,
          y: y * yScale,
          z: z * zScale,
          name: `Gen${ind.generation ?? '?'}-Ind${ind.individual_index ?? ind.index ?? '?'}`,
          key: ind.key,
          originalIndex: idx,
          rank: rank + 1
        });

        allZValues.push(z * zScale);
      }
    }

    if (allZValues.length === 0) {
      return { data, zMin: 0, zMax: 1 };
    }
    return { data, zMin: Math.min(...allZValues), zMax: Math.max(...allZValues) };
  }, [fronts, filteredIndividuals, comp, zComp, upToRankNo]);

  // 为每个点添加颜色
  const coloredData = useMemo(() => {
    if (chartData.data.length === 0) return [];
    const { zMin, zMax } = chartData;
    
    return chartData.data.map(point => ({
      ...point,
      fill: valueToColor(point.z, zMin, zMax)
    }));
  }, [chartData]);

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

      {/* 图表 */}
      <div className="flex gap-4 items-start">
        <div className="flex-1">
          <ResponsiveContainer width="100%" height={400}>
            <ScatterChart margin={{ top: 20, right: 20, bottom: 20, left: 20 }}>
              <CartesianGrid strokeDasharray="3 3" />
              <XAxis
                type="number"
                dataKey="x"
                name={xLabel}
                label={{ value: xLabel, position: "insideBottom", offset: -5 }}
              />
              <YAxis
                type="number"
                dataKey="y"
                name={yLabel}
                label={{ value: yLabel, angle: -90, position: "insideLeft" }}
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
                      </div>
                    );
                  }
                  return null;
                }}
              />
              <Scatter
                name="Pareto解"
                data={coloredData}
                shape={(props: any) => {
                  const { cx, cy, payload } = props;
                  if (!payload || !payload.fill || cx === undefined || cy === undefined) return null;
                  return (
                    <circle
                      cx={cx}
                      cy={cy}
                      r={5}
                      fill={payload.fill}
                      stroke="rgba(0,0,0,0.1)"
                      strokeWidth={0.5}
                      opacity={0.8}
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
        />
      </div>
    </div>
  );
}

