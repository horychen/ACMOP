/**
 * 非支配排序工具函数
 * 实现类似 pygmo.fast_non_dominated_sorting 的功能
 */

export interface Point {
  f1: number;
  f2: number;
  f3: number;
  [key: string]: any;
}

/**
 * 判断个体 i 是否支配个体 j
 * 对于最小化问题，如果 i 的所有目标值都 <= j，且至少有一个 < j，则 i 支配 j
 */
function dominates(pointI: number[], pointJ: number[]): boolean {
  const allLessEqual = pointI.every((val, idx) => val <= pointJ[idx]);
  const atLeastOneLess = pointI.some((val, idx) => val < pointJ[idx]);
  return allLessEqual && atLeastOneLess;
}

/**
 * 快速非支配排序
 * 返回所有前沿的索引数组
 */
export function fastNonDominatedSorting(points: number[][]): number[][] {
  const n = points.length;
  if (n === 0) {
    return [];
  }

  // 存储每个个体被哪些个体支配
  const dominatedBy: number[][] = Array(n).fill(null).map(() => []);
  // 存储每个个体支配多少个其他个体
  const dominationCount: number[] = Array(n).fill(0);
  // 存储每个个体属于哪个前沿
  const rank: number[] = Array(n).fill(-1);
  // 存储每个前沿包含的个体索引
  const fronts: number[][] = [];

  // 第一遍：计算支配关系
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      if (i === j) continue;
      
      if (dominates(points[i], points[j])) {
        dominatedBy[i].push(j);
        dominationCount[j]++;
      }
    }
  }

  // 找到第一前沿（Rank 1）：所有 dominationCount 为 0 的个体
  let currentFront: number[] = [];
  for (let i = 0; i < n; i++) {
    if (dominationCount[i] === 0) {
      rank[i] = 0;
      currentFront.push(i);
    }
  }

  let frontIndex = 0;
  while (currentFront.length > 0) {
    fronts.push([...currentFront]);
    const nextFront: number[] = [];

    // 对于当前前沿的每个个体，减少被它支配的个体的 dominationCount
    for (const i of currentFront) {
      for (const j of dominatedBy[i]) {
        dominationCount[j]--;
        if (dominationCount[j] === 0 && rank[j] === -1) {
          rank[j] = frontIndex + 1;
          nextFront.push(j);
        }
      }
    }

    currentFront = nextFront;
    frontIndex++;
  }

  return fronts;
}

/**
 * 从个体数据中提取目标函数值数组
 */
export function extractObjectives(individuals: Point[]): number[][] {
  return individuals.map(ind => [
    ind.f1 ?? ind.objectives?.f1 ?? 0,
    ind.f2 ?? ind.objectives?.f2 ?? 0,
    ind.f3 ?? ind.objectives?.f3 ?? 0
  ]);
}

/**
 * 根据 f3 值过滤个体
 */
export function filterByF3(individuals: Point[], maxF3?: number): Point[] {
  if (maxF3 === undefined || maxF3 === null) {
    return individuals;
  }
  return individuals.filter(ind => {
    const f3 = ind.f3 ?? ind.objectives?.f3;
    return f3 !== undefined && f3 < maxF3;
  });
}

