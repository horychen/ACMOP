"use client"

import { useState } from "react";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { Checkbox } from "@/components/ui/checkbox";
import Editor from "@monaco-editor/react";
import { Trash2, ChevronDown, ChevronUp } from "lucide-react";
import { Button } from "@/components/ui/button";

export interface ParameterConfig {
  name: string;
  type: "fixed" | "free" | "derived";
  value: any;
  bounds?: [number, number] | null;
  calc?: string | null; // Python code as string
  calc_bounds?: string | null; // Python code for calculating bounds
  args?: string[]; // Arguments for calc/calc_bounds
  unit?: string;
  comment?: string;
  sensitivityAnalysis?: boolean;
}

interface ParameterEditorProps {
  parameter: ParameterConfig;
  index: number;
  onChange: (index: number, param: ParameterConfig) => void;
  onDelete: (index: number) => void;
  availableParameters: string[]; // For derived calc function parameters
}

const getDefaultCalcCode = (paramName: string, availableParams: string[]): string => {
  return `# Python lambda function
# Example: lambda x, y: x * y
# Available parameters: ${availableParams.length > 0 ? availableParams.join(", ") : "无"}

lambda ${paramName}: ${paramName} * 1.0`;
};

export function ParameterEditor({
  parameter,
  index,
  onChange,
  onDelete,
  availableParameters
}: ParameterEditorProps) {
  const [isExpanded, setIsExpanded] = useState(true);

  const handleChange = (field: keyof ParameterConfig, value: any) => {
    const updated = { ...parameter, [field]: value };
    
    // Reset dependent fields when type changes
    if (field === "type") {
      if (value === "fixed") {
        updated.bounds = null;
        updated.calc = null;
      } else if (value === "free") {
        updated.calc = null;
        if (!updated.bounds) {
          updated.bounds = [0, 100];
        }
      } else if (value === "derived") {
        updated.bounds = null;
        if (!updated.calc) {
          // Generate default calc code
          const paramName = parameter.name || "x";
          updated.calc = getDefaultCalcCode(paramName, availableParameters);
        }
      }
    }
    
    onChange(index, updated);
  };

  const handleBoundsChange = (index: 0 | 1, value: string) => {
    const numValue = parseFloat(value);
    if (!isNaN(numValue)) {
      const newBounds: [number, number] = parameter.bounds || [0, 100];
      newBounds[index] = numValue;
      handleChange("bounds", newBounds);
    }
  };

  return (
    <Card className="border-l-4 border-l-primary">
      <CardHeader className="pb-3">
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-2 flex-1">
            <Button
              variant="ghost"
              size="sm"
              onClick={() => setIsExpanded(!isExpanded)}
              className="h-6 w-6 p-0"
            >
              {isExpanded ? (
                <ChevronDown className="h-4 w-4" />
              ) : (
                <ChevronUp className="h-4 w-4" />
              )}
            </Button>
            <CardTitle className="text-base font-medium">
              {parameter.name || `参数 ${index + 1}`}
            </CardTitle>
          </div>
          <Button
            variant="ghost"
            size="sm"
            onClick={() => onDelete(index)}
            className="text-destructive hover:text-destructive"
          >
            <Trash2 className="h-4 w-4" />
          </Button>
        </div>
      </CardHeader>
      
      {isExpanded && (
        <CardContent className="space-y-4">
          {/* 参数名称 */}
          <div className="space-y-2">
            <Label htmlFor={`param-name-${index}`}>参数名称</Label>
            <Input
              id={`param-name-${index}`}
              value={parameter.name || ""}
              onChange={(e) => handleChange("name", e.target.value)}
              placeholder="例如: mm_r_so"
            />
          </div>

          {/* 参数类型 */}
          <div className="space-y-2">
            <Label htmlFor={`param-type-${index}`}>参数类型</Label>
            <Select
              value={parameter.type}
              onValueChange={(value: "fixed" | "free" | "derived") => handleChange("type", value)}
            >
              <SelectTrigger id={`param-type-${index}`}>
                <SelectValue />
              </SelectTrigger>
              <SelectContent>
                <SelectItem value="fixed">Fixed (固定值)</SelectItem>
                <SelectItem value="free">Free (自由变量，用于优化)</SelectItem>
                <SelectItem value="derived">Derived (派生值，由计算得出)</SelectItem>
              </SelectContent>
            </Select>
          </div>

          {/* 参数值 - 仅当类型为 fixed 或 free 时显示 */}
          {(parameter.type === "fixed" || parameter.type === "free") && (
            <div className="space-y-2">
              <Label htmlFor={`param-value-${index}`}>
                参数值 {parameter.unit && `(${parameter.unit})`}
              </Label>
              <Input
                id={`param-value-${index}`}
                type="number"
                step="any"
                value={parameter.value !== null && parameter.value !== undefined ? String(parameter.value) : ""}
                onChange={(e) => {
                  const val = e.target.value === "" ? null : parseFloat(e.target.value);
                  handleChange("value", isNaN(val as number) ? null : val);
                }}
                placeholder="输入数值"
              />
            </div>
          )}

          {/* Bounds - 仅当类型为 free 时显示 */}
          {parameter.type === "free" && (
            <div className="space-y-2">
              <Label>优化边界 (Bounds)</Label>
              <div className="grid grid-cols-2 gap-2">
                <div className="space-y-1">
                  <Label htmlFor={`param-bound-lower-${index}`} className="text-xs text-muted-foreground">
                    下界 (Lower)
                  </Label>
                  <Input
                    id={`param-bound-lower-${index}`}
                    type="number"
                    step="any"
                    value={parameter.bounds ? String(parameter.bounds[0]) : ""}
                    onChange={(e) => handleBoundsChange(0, e.target.value)}
                    placeholder="最小值"
                  />
                </div>
                <div className="space-y-1">
                  <Label htmlFor={`param-bound-upper-${index}`} className="text-xs text-muted-foreground">
                    上界 (Upper)
                  </Label>
                  <Input
                    id={`param-bound-upper-${index}`}
                    type="number"
                    step="any"
                    value={parameter.bounds ? String(parameter.bounds[1]) : ""}
                    onChange={(e) => handleBoundsChange(1, e.target.value)}
                    placeholder="最大值"
                  />
                </div>
              </div>
            </div>
          )}

          {/* Calc 方法编辑器 - 仅当类型为 derived 时显示 */}
          {parameter.type === "derived" && (
            <div className="space-y-2">
              <Label>计算函数 (Python Lambda)</Label>
              <div className="border rounded-lg overflow-hidden" style={{ height: "200px" }}>
                <Editor
                  height="200px"
                  defaultLanguage="python"
                  value={parameter.calc || ""}
                  onChange={(value) => handleChange("calc", value || "")}
                  theme="vs-dark"
                  options={{
                    minimap: { enabled: false },
                    fontSize: 12,
                    wordWrap: "on",
                    automaticLayout: true,
                  }}
                />
              </div>
              <p className="text-xs text-muted-foreground">
                提示: 使用 lambda 函数格式，例如: <code>lambda x, y: x * y</code>
                <br />
                可用参数: {availableParameters.length > 0 ? availableParameters.join(", ") : "无"}
              </p>
            </div>
          )}

          {/* 单位 */}
          <div className="space-y-2">
            <Label htmlFor={`param-unit-${index}`}>单位</Label>
            <Input
              id={`param-unit-${index}`}
              value={parameter.unit || ""}
              onChange={(e) => handleChange("unit", e.target.value)}
              placeholder="例如: mm, deg, A"
            />
          </div>

          {/* 注释 */}
          <div className="space-y-2">
            <Label htmlFor={`param-comment-${index}`}>注释</Label>
            <Input
              id={`param-comment-${index}`}
              value={parameter.comment || ""}
              onChange={(e) => handleChange("comment", e.target.value)}
              placeholder="参数说明"
            />
          </div>

          {/* 敏感性分析选项 */}
          {parameter.type === "free" && (
            <div className="flex items-center space-x-2">
              <Checkbox
                id={`param-sensitivity-${index}`}
                checked={parameter.sensitivityAnalysis || false}
                onCheckedChange={(checked) => handleChange("sensitivityAnalysis", checked)}
              />
              <Label
                htmlFor={`param-sensitivity-${index}`}
                className="text-sm font-normal cursor-pointer"
              >
                请求敏感性分析
              </Label>
            </div>
          )}
        </CardContent>
      )}
    </Card>
  );
}

