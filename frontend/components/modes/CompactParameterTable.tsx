"use client"

import { useState } from "react";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { Checkbox } from "@/components/ui/checkbox";
import Editor from "@monaco-editor/react";
import { ParameterConfig } from "./ParameterEditor";

interface CompactParameterTableProps {
  parameters: ParameterConfig[];
  onParameterChange: (index: number, param: ParameterConfig) => void;
  availableParameters: (index: number) => string[];
}

export function CompactParameterTable({
  parameters,
  onParameterChange,
  availableParameters
}: CompactParameterTableProps) {
  const [editingCalcIndex, setEditingCalcIndex] = useState<number | null>(null);
  const [editingCalcBoundsIndex, setEditingCalcBoundsIndex] = useState<number | null>(null);

  const handleChange = (index: number, field: keyof ParameterConfig, value: any) => {
    const param = parameters[index];
    const updated = { ...param, [field]: value };
    
    // 当类型改变时，重置相关字段
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
          const paramName = param.name || "x";
          updated.calc = `lambda ${paramName}: ${paramName} * 1.0`;
        }
      }
    }
    
    onParameterChange(index, updated);
  };

  return (
    <div className="space-y-4">
      {/* Fixed Variables Section */}
      <div>
        <h3 className="text-sm font-semibold mb-2 text-muted-foreground"># Fixed variables</h3>
        <div className="space-y-1">
          {parameters
            .map((param, originalIndex) => ({ param, originalIndex }))
            .filter(({ param }) => param.type === "fixed")
            .map(({ param, originalIndex }) => {
              return (
                <div key={originalIndex} className="grid grid-cols-12 gap-2 items-center text-xs font-mono py-1 hover:bg-muted/50 rounded px-2">
                  <div className="col-span-3">
                    <Input
                      value={param.name}
                      onChange={(e) => handleChange(originalIndex, "name", e.target.value)}
                      placeholder="param_name"
                      className="h-7 text-xs font-mono"
                    />
                  </div>
                  <div className="col-span-1">
                    <Select
                      value={param.type}
                      onValueChange={(v: "fixed" | "free" | "derived") => handleChange(originalIndex, "type", v)}
                    >
                      <SelectTrigger className="h-7 text-xs">
                        <SelectValue />
                      </SelectTrigger>
                      <SelectContent>
                        <SelectItem value="fixed">fixed</SelectItem>
                        <SelectItem value="free">free</SelectItem>
                        <SelectItem value="derived">derived</SelectItem>
                      </SelectContent>
                    </Select>
                  </div>
                  <div className="col-span-2">
                    <Input
                      type="number"
                      step="any"
                      value={param.value !== null && param.value !== undefined ? String(param.value) : ""}
                      onChange={(e) => handleChange(originalIndex, "value", e.target.value === "" ? null : parseFloat(e.target.value))}
                      placeholder="value"
                      className="h-7 text-xs"
                    />
                  </div>
                  <div className="col-span-1">
                    <Input
                      value={param.unit || "mm"}
                      onChange={(e) => handleChange(originalIndex, "unit", e.target.value)}
                      placeholder="unit"
                      className="h-7 text-xs"
                    />
                  </div>
                  <div className="col-span-5 text-xs text-muted-foreground truncate">
                    {param.comment || ""}
                  </div>
                </div>
              );
            })}
        </div>
      </div>

      {/* Free Variables Section */}
      <div>
        <h3 className="text-sm font-semibold mb-2 text-muted-foreground"># Free variables</h3>
        <div className="space-y-1">
          {parameters
            .map((param, originalIndex) => ({ param, originalIndex }))
            .filter(({ param }) => param.type === "free")
            .map(({ param, originalIndex }) => {
              const isEditingCalcBounds = editingCalcBoundsIndex === originalIndex;
              return (
                <div key={originalIndex} className="space-y-1">
                  <div className="grid grid-cols-12 gap-2 items-center text-xs font-mono py-1 hover:bg-muted/50 rounded px-2">
                    <div className="col-span-2">
                      <Input
                        value={param.name}
                        onChange={(e) => handleChange(originalIndex, "name", e.target.value)}
                        placeholder="param_name"
                        className="h-7 text-xs font-mono"
                      />
                    </div>
                    <div className="col-span-1">
                      <Select
                        value={param.type}
                        onValueChange={(v: "fixed" | "free" | "derived") => handleChange(originalIndex, "type", v)}
                      >
                        <SelectTrigger className="h-7 text-xs">
                          <SelectValue />
                        </SelectTrigger>
                        <SelectContent>
                          <SelectItem value="fixed">fixed</SelectItem>
                          <SelectItem value="free">free</SelectItem>
                          <SelectItem value="derived">derived</SelectItem>
                        </SelectContent>
                      </Select>
                    </div>
                    <div className="col-span-1">
                      <Input
                        type="number"
                        step="any"
                        value={param.value !== null && param.value !== undefined ? String(param.value) : ""}
                        onChange={(e) => handleChange(originalIndex, "value", e.target.value === "" ? null : parseFloat(e.target.value))}
                        placeholder="value"
                        className="h-7 text-xs"
                      />
                    </div>
                    <div className="col-span-2 flex gap-1">
                      <Input
                        type="number"
                        step="any"
                        value={param.bounds ? String(param.bounds[0]) : ""}
                        onChange={(e) => {
                          const bounds = param.bounds || [0, 100];
                          bounds[0] = parseFloat(e.target.value) || 0;
                          handleChange(originalIndex, "bounds", bounds);
                        }}
                        placeholder="min"
                        className="h-7 text-xs"
                      />
                      <Input
                        type="number"
                        step="any"
                        value={param.bounds ? String(param.bounds[1]) : ""}
                        onChange={(e) => {
                          const bounds = param.bounds || [0, 100];
                          bounds[1] = parseFloat(e.target.value) || 100;
                          handleChange(originalIndex, "bounds", bounds);
                        }}
                        placeholder="max"
                        className="h-7 text-xs"
                      />
                    </div>
                    <div className="col-span-1">
                      <Input
                        value={param.unit || "mm"}
                        onChange={(e) => handleChange(originalIndex, "unit", e.target.value)}
                        placeholder="unit"
                        className="h-7 text-xs"
                      />
                    </div>
                    <div className="col-span-3 flex items-center gap-2">
                      {param.calc_bounds ? (
                        <button
                          onClick={() => setEditingCalcBoundsIndex(isEditingCalcBounds ? null : originalIndex)}
                          className="text-xs px-2 py-1 bg-blue-100 dark:bg-blue-900 text-blue-700 dark:text-blue-300 rounded hover:bg-blue-200 dark:hover:bg-blue-800"
                          title="点击查看 calc_bounds 计算方法代码"
                        >
                          {isEditingCalcBounds ? "收起代码" : "查看 calc_bounds"}
                        </button>
                      ) : (
                        <span className="text-xs text-muted-foreground">无 calc_bounds</span>
                      )}
                      <Checkbox
                        checked={param.sensitivityAnalysis || false}
                        onCheckedChange={(checked) => handleChange(originalIndex, "sensitivityAnalysis", checked)}
                      />
                      <span className="text-xs text-muted-foreground">Sens.</span>
                    </div>
                    <div className="col-span-2 text-xs text-muted-foreground truncate">
                      {param.comment || ""}
                    </div>
                  </div>
                  {isEditingCalcBounds && param.calc_bounds && (
                    <div className="ml-4 border rounded-lg overflow-hidden bg-slate-50 dark:bg-slate-900" style={{ height: "120px" }}>
                      <div className="px-2 py-1 text-xs font-semibold bg-slate-200 dark:bg-slate-800">
                        calc_bounds (args: {param.args?.join(", ") || "无"})
                      </div>
                      <Editor
                        height="90px"
                        defaultLanguage="python"
                        value={param.calc_bounds}
                        onChange={(value) => handleChange(originalIndex, "calc_bounds", value || "")}
                        theme="vs-dark"
                        options={{
                          minimap: { enabled: false },
                          fontSize: 11,
                          wordWrap: "on",
                          automaticLayout: true,
                          readOnly: true, // 只读，用于检查代码
                        }}
                      />
                    </div>
                  )}
                </div>
              );
            })}
        </div>
      </div>

      {/* Derived Variables Section */}
      <div>
        <h3 className="text-sm font-semibold mb-2 text-muted-foreground"># Derived variables</h3>
        <div className="space-y-1">
          {parameters
            .map((param, originalIndex) => ({ param, originalIndex }))
            .filter(({ param }) => param.type === "derived")
            .map(({ param, originalIndex }) => {
              const isEditing = editingCalcIndex === originalIndex;
              return (
                <div key={originalIndex} className="space-y-1">
                  <div className="grid grid-cols-12 gap-2 items-center text-xs font-mono py-1 hover:bg-muted/50 rounded px-2">
                    <div className="col-span-3">
                      <Input
                        value={param.name}
                        onChange={(e) => handleChange(originalIndex, "name", e.target.value)}
                        placeholder="param_name"
                        className="h-7 text-xs font-mono"
                      />
                    </div>
                    <div className="col-span-1">
                      <Select
                        value={param.type}
                        onValueChange={(v: "fixed" | "free" | "derived") => handleChange(originalIndex, "type", v)}
                      >
                        <SelectTrigger className="h-7 text-xs">
                          <SelectValue />
                        </SelectTrigger>
                        <SelectContent>
                          <SelectItem value="fixed">fixed</SelectItem>
                          <SelectItem value="free">free</SelectItem>
                          <SelectItem value="derived">derived</SelectItem>
                        </SelectContent>
                      </Select>
                    </div>
                    <div className="col-span-1">
                      <Input
                        value={param.unit || "mm"}
                        onChange={(e) => handleChange(originalIndex, "unit", e.target.value)}
                        placeholder="unit"
                        className="h-7 text-xs"
                      />
                    </div>
                    <div className="col-span-7 flex items-center gap-2">
                      {param.calc ? (
                        <>
                          <button
                            onClick={() => setEditingCalcIndex(isEditing ? null : originalIndex)}
                            className="text-xs px-2 py-1 bg-green-100 dark:bg-green-900 text-green-700 dark:text-green-300 rounded hover:bg-green-200 dark:hover:bg-green-800"
                            title="点击查看 calc 计算方法代码"
                          >
                            {isEditing ? "收起代码" : "查看 calc"}
                          </button>
                          <span className="text-xs text-muted-foreground truncate" title={param.calc}>
                            {param.calc.length > 50 ? param.calc.substring(0, 50) + "..." : param.calc}
                          </span>
                        </>
                      ) : (
                        <span className="text-xs text-muted-foreground">无 calc</span>
                      )}
                    </div>
                  </div>
                  {isEditing && param.calc && (
                    <div className="ml-4 border-2 border-green-300 dark:border-green-700 rounded-lg overflow-hidden bg-slate-50 dark:bg-slate-900 shadow-md" style={{ height: "180px" }}>
                      <div className="px-3 py-2 text-xs font-semibold bg-green-200 dark:bg-green-800 text-green-900 dark:text-green-100 flex items-center justify-between">
                        <span>calc 计算方法</span>
                        <span className="text-xs font-normal text-green-700 dark:text-green-300">
                          args: [{param.args?.join(", ") || "无"}]
                        </span>
                      </div>
                      <div style={{ height: "140px" }}>
                        <Editor
                          height="140px"
                          defaultLanguage="python"
                          value={param.calc}
                          onChange={(value) => handleChange(originalIndex, "calc", value || "")}
                          theme="vs-dark"
                          options={{
                            minimap: { enabled: false },
                            fontSize: 12,
                            wordWrap: "on",
                            automaticLayout: true,
                            readOnly: true, // 只读，用于检查代码
                            lineNumbers: "on",
                          }}
                        />
                      </div>
                    </div>
                  )}
                </div>
              );
            })}
        </div>
      </div>
    </div>
  );
}

