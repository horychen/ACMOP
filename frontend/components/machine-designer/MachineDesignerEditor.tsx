"use client"

import React, { useState, useEffect } from "react";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Switch } from "@/components/ui/switch";
import { Button } from "@/components/ui/button";
import { Separator } from "@/components/ui/separator";
import { ScrollArea } from "@/components/ui/scroll-area";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { Save, Upload, Download, RefreshCw, Plus, Trash2, ChevronDown, ChevronUp, Code } from "lucide-react";
import { Checkbox } from "@/components/ui/checkbox";
import Editor from "@monaco-editor/react";
import { Badge } from "@/components/ui/badge";
import { Loader2 } from "lucide-react";

interface ParameterData {
  name: string;
  type: "fixed" | "free" | "derived";
  value: any;
  bounds: [number, number] | null;
  unit: string;
  comment: string | null;
  calc_dependencies: any;
}

interface MachineDesignerData {
  machine_class: string;
  bool_PermanentMagnet: boolean;
  bool_StatorSlotClosed: boolean;
  bool_RotorNotched: boolean;
  [key: string]: any; // 允许动态参数
  machineGeometry: {
    rotorCore: {
      color: string | null;
      GP: Record<string, string>; // 从数组改为对象
      _needs_rebuild: boolean;
    };
    shaft: {
      color: string | null;
      GP: Record<string, string>;
      _needs_rebuild: boolean;
    };
    rotorMagnet: {
      color: string | null;
      GP: Record<string, string>;
      _needs_rebuild: boolean;
    };
    statorCore: {
      color: string | null;
      GP: Record<string, string>;
      _needs_rebuild: boolean;
    };
    coils: {
      color: string | null;
      GP: Record<string, string>;
      _needs_rebuild: boolean;
    };
  };
  wily: {
    phase_number_m: number;
    stator_slot_number_Qs: number;
    pole_pair_number_p: number;
    suspension_pole_pair_number_ps: number;
  };
}

const DEFAULT_DATA: MachineDesignerData = {
  machine_class: "bearingless_spmsm_heart.bearingless_spmsm_design_variant",
  bool_PermanentMagnet: true,
  bool_StatorSlotClosed: false,
  bool_RotorNotched: true,
  wily: {
    phase_number_m: 3,
    stator_slot_number_Qs: 12,
    pole_pair_number_p: 4,
    suspension_pole_pair_number_ps: 5
  },
  machineGeometry: {
    rotorCore: {
      color: null,
      GP: {
        "mm_r_ro": "rotor_outer_radius",
        "mm_d_ri": "rotor_iron (back iron) depth",
        "mm_d_pm": "magnet_depth",
        "mm_d_rp": "inter_polar_iron_thickness",
        "mm_d_rs": "inter_segment_iron_thickness",
        "p": "pole_pair_number_p",
        "s": "number_of_magnet_segments_per_pole"
      },
      _needs_rebuild: true
    },
    shaft: {
      color: null,
      GP: {
        "mm_r_ri": "rotor_inner_radius"
      },
      _needs_rebuild: true
    },
    rotorMagnet: {
      color: null,
      GP: {
        "mm_d_pm": "magnet_depth",
        "mm_d_ri": "rotor_iron (back iron) depth",
        "mm_r_ri": "rotor_inner_radius"
      },
      _needs_rebuild: true
    },
    statorCore: {
      color: null,
      GP: {
        "mm_r_si": "stator_inner_radius",
        "mm_d_sto": "stator_tooth_open_depth",
        "mm_d_sts": "stator_tooth_shoe_depth",
        "mm_d_st": "stator_tooth_depth",
        "mm_d_sy": "stator_yoke_depth",
        "mm_w_st": "stator_tooth_width",
        "deg_alpha_st": "stator_tooth_span_angle",
        "deg_alpha_sto": "stator_tooth_open_angle",
        "Q": "stator_slot_number_Qs"
      },
      _needs_rebuild: true
    },
    coils: {
      color: null,
      GP: {
        "mm_r_so": "stator_outer_radius",
        "mm_d_sy": "stator_yoke_depth",
        "mm_w_st": "stator_tooth_width",
        "mm_d_st": "stator_tooth_depth"
      },
      _needs_rebuild: true
    }
  }
};

interface ParameterInfo {
  name: string;
  displayName: string;
  type: string;
  unit: string;
  value: any;
  bounds: [number, number] | null;
  calc: string | null;
  calc_bounds: string | null;
  args: string[];
  comment: string;
}

export default function MachineDesignerEditor() {
  const [data, setData] = useState<MachineDesignerData>(DEFAULT_DATA);
  const [activeTab, setActiveTab] = useState<string>("general");
  const [parametersInfo, setParametersInfo] = useState<Record<string, ParameterInfo>>({});
  const [isLoadingParams, setIsLoadingParams] = useState(false);
  const [expandedGP, setExpandedGP] = useState<Record<string, boolean>>({});
  const [editingCalc, setEditingCalc] = useState<Record<string, "calc" | "calc_bounds" | null>>({});

  useEffect(() => {
    // 可以在这里加载保存的数据
    loadSavedData();
    loadParametersInfo();
  }, []);

  const loadParametersInfo = async () => {
    setIsLoadingParams(true);
    try {
      const response = await fetch('/api/acmopv2/all-parameters');
      if (!response.ok) {
        throw new Error('加载参数信息失败');
      }
      const apiData = await response.json();
      
      // 构建参数信息映射
      const paramsMap: Record<string, ParameterInfo> = {};
      
      // 合并所有类型的参数
      const allParams = [
        ...(apiData.parameters.fixed || []),
        ...(apiData.parameters.free || []),
        ...(apiData.parameters.derived || [])
      ];
      
      allParams.forEach((param: any) => {
        paramsMap[param.name] = {
          name: param.name,
          displayName: param.displayName || param.name,
          type: param.type,
          unit: param.unit || "mm",
          value: param.value,
          bounds: param.bounds || null,
          calc: param.calc || null,
          calc_bounds: param.calc_bounds || null,
          args: param.args || [],
          comment: param.comment || ""
        };
      });
      
      setParametersInfo(paramsMap);
    } catch (err) {
      console.error("加载参数信息失败:", err);
    } finally {
      setIsLoadingParams(false);
    }
  };

  const loadSavedData = async () => {
    try {
      // 首先尝试从本地存储加载
      const saved = localStorage.getItem("machine_designer_data");
      if (saved) {
        const parsed = JSON.parse(saved);
        // 验证数据结构，确保GP是对象而不是数组
        if (parsed.machineGeometry) {
          Object.keys(parsed.machineGeometry).forEach(component => {
            if (Array.isArray(parsed.machineGeometry[component].GP)) {
              // 如果是数组，转换为对象格式
              const gpArray = parsed.machineGeometry[component].GP;
              parsed.machineGeometry[component].GP = {};
              gpArray.forEach((gp: string) => {
                parsed.machineGeometry[component].GP[gp] = gp;
              });
            }
          });
        }
        setData(parsed);
        return;
      }
      
      // 如果没有本地数据，尝试从后端加载默认JSON
      try {
        const response = await fetch('/api/machine-designer/default');
        if (response.ok) {
          const defaultData = await response.json();
          setData(defaultData);
        }
      } catch (err) {
        console.log("无法从后端加载默认数据，使用本地默认值");
      }
    } catch (e) {
      console.error("加载保存的数据失败:", e);
    }
  };

  const saveData = () => {
    try {
      localStorage.setItem("machine_designer_data", JSON.stringify(data, null, 2));
      alert("数据已保存到本地存储");
    } catch (e) {
      alert("保存失败: " + e);
    }
  };

  const downloadJSON = () => {
    const blob = new Blob([JSON.stringify(data, null, 2)], { type: "application/json" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = "machine_designer.json";
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  };

  const handleFileUpload = (event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0];
    if (file) {
      const reader = new FileReader();
      reader.onload = (e) => {
        try {
          const jsonData = JSON.parse(e.target?.result as string);
          setData(jsonData);
          alert("文件加载成功");
        } catch (err) {
          alert("文件解析失败: " + err);
        }
      };
      reader.readAsText(file);
    }
  };

  const updateWinding = (field: string, value: number) => {
    setData(prev => ({
      ...prev,
      wily: {
        ...prev.wily,
        [field]: value
      }
    }));
  };

  const updateGeometryGP = (component: string, gpObject: Record<string, string>) => {
    setData(prev => ({
      ...prev,
      machineGeometry: {
        ...prev.machineGeometry,
        [component]: {
          ...prev.machineGeometry[component as keyof typeof prev.machineGeometry],
          GP: gpObject
        }
      }
    }));
  };

  const addGPToComponent = (component: string, gpKey: string, gpDisplayName: string) => {
    const currentGP = data.machineGeometry[component as keyof typeof data.machineGeometry].GP;
    if (!(gpKey in currentGP)) {
      updateGeometryGP(component, { ...currentGP, [gpKey]: gpDisplayName });
    }
  };

  const removeGPFromComponent = (component: string, gpKey: string) => {
    const currentGP = data.machineGeometry[component as keyof typeof data.machineGeometry].GP;
    const newGP = { ...currentGP };
    delete newGP[gpKey];
    updateGeometryGP(component, newGP);
  };

  const updateGPKey = (component: string, oldKey: string, newKey: string) => {
    const currentGP = data.machineGeometry[component as keyof typeof data.machineGeometry].GP;
    if (oldKey in currentGP) {
      const displayName = currentGP[oldKey];
      const newGP = { ...currentGP };
      delete newGP[oldKey];
      newGP[newKey] = displayName;
      updateGeometryGP(component, newGP);
    }
  };

  const updateGPDisplayName = (component: string, gpKey: string, displayName: string) => {
    const currentGP = data.machineGeometry[component as keyof typeof data.machineGeometry].GP;
    updateGeometryGP(component, { ...currentGP, [gpKey]: displayName });
  };

  const resetToDefault = () => {
    if (confirm("确定要重置为默认值吗？")) {
      setData(DEFAULT_DATA);
    }
  };

  return (
    <div className="space-y-6">
      <Card>
        <CardHeader>
          <div className="flex items-center justify-between">
            <div>
              <CardTitle>Machine Designer 配置</CardTitle>
              <CardDescription>
                编辑机器设计器的完整配置，对应 machine_designer.json 结构
              </CardDescription>
            </div>
            <div className="flex gap-2">
              <Button onClick={resetToDefault} variant="outline" size="sm">
                <RefreshCw className="mr-2 h-4 w-4" />
                重置
              </Button>
              <Button onClick={saveData} variant="outline" size="sm">
                <Save className="mr-2 h-4 w-4" />
                保存
              </Button>
              <Button onClick={downloadJSON} variant="outline" size="sm">
                <Download className="mr-2 h-4 w-4" />
                下载JSON
              </Button>
              <label className="cursor-pointer">
                <Button variant="outline" size="sm" asChild>
                  <span>
                    <Upload className="mr-2 h-4 w-4" />
                    上传JSON
                  </span>
                </Button>
                <input
                  type="file"
                  accept=".json"
                  onChange={handleFileUpload}
                  className="hidden"
                />
              </label>
            </div>
          </div>
        </CardHeader>
        <CardContent>
          <Tabs value={activeTab} onValueChange={setActiveTab} className="w-full">
            <TabsList className="grid w-full grid-cols-4">
              <TabsTrigger value="general">基本设置</TabsTrigger>
              <TabsTrigger value="winding">绕组参数</TabsTrigger>
              <TabsTrigger value="geometry">几何参数</TabsTrigger>
              <TabsTrigger value="json">JSON视图</TabsTrigger>
            </TabsList>

            {/* 基本设置 */}
            <TabsContent value="general" className="space-y-4 mt-4">
              <div className="space-y-4">
                <div className="space-y-2">
                  <Label htmlFor="machine_class">机器类型 (machine_class)</Label>
                  <Input
                    id="machine_class"
                    value={data.machine_class}
                    onChange={(e) => setData(prev => ({ ...prev, machine_class: e.target.value }))}
                    placeholder="例如: bearingless_spmsm_heart.bearingless_spmsm_design_variant"
                  />
                </div>

                <Separator />

                <div className="space-y-4">
                  <h3 className="text-lg font-semibold">布尔标志</h3>
                  <div className="space-y-3">
                    <div className="flex items-center justify-between">
                      <Label htmlFor="bool_PermanentMagnet">永久磁铁 (bool_PermanentMagnet)</Label>
                      <Switch
                        id="bool_PermanentMagnet"
                        checked={data.bool_PermanentMagnet}
                        onCheckedChange={(checked) => setData(prev => ({ ...prev, bool_PermanentMagnet: checked }))}
                      />
                    </div>
                    <div className="flex items-center justify-between">
                      <Label htmlFor="bool_StatorSlotClosed">定子槽闭合 (bool_StatorSlotClosed)</Label>
                      <Switch
                        id="bool_StatorSlotClosed"
                        checked={data.bool_StatorSlotClosed}
                        onCheckedChange={(checked) => setData(prev => ({ ...prev, bool_StatorSlotClosed: checked }))}
                      />
                    </div>
                    <div className="flex items-center justify-between">
                      <Label htmlFor="bool_RotorNotched">转子开槽 (bool_RotorNotched)</Label>
                      <Switch
                        id="bool_RotorNotched"
                        checked={data.bool_RotorNotched}
                        onCheckedChange={(checked) => setData(prev => ({ ...prev, bool_RotorNotched: checked }))}
                      />
                    </div>
                  </div>
                </div>
              </div>
            </TabsContent>

            {/* 绕组参数 */}
            <TabsContent value="winding" className="space-y-4 mt-4">
              <div className="grid grid-cols-2 gap-4">
                <div className="space-y-2">
                  <Label htmlFor="phase_number_m">相数 (phase_number_m)</Label>
                  <Input
                    id="phase_number_m"
                    type="number"
                    value={data.wily.phase_number_m}
                    onChange={(e) => updateWinding("phase_number_m", parseInt(e.target.value) || 0)}
                  />
                </div>
                <div className="space-y-2">
                  <Label htmlFor="stator_slot_number_Qs">定子槽数 (stator_slot_number_Qs)</Label>
                  <Input
                    id="stator_slot_number_Qs"
                    type="number"
                    value={data.wily.stator_slot_number_Qs}
                    onChange={(e) => updateWinding("stator_slot_number_Qs", parseInt(e.target.value) || 0)}
                  />
                </div>
                <div className="space-y-2">
                  <Label htmlFor="pole_pair_number_p">极对数 (pole_pair_number_p)</Label>
                  <Input
                    id="pole_pair_number_p"
                    type="number"
                    value={data.wily.pole_pair_number_p}
                    onChange={(e) => updateWinding("pole_pair_number_p", parseInt(e.target.value) || 0)}
                  />
                </div>
                <div className="space-y-2">
                  <Label htmlFor="suspension_pole_pair_number_ps">悬浮极对数 (suspension_pole_pair_number_ps)</Label>
                  <Input
                    id="suspension_pole_pair_number_ps"
                    type="number"
                    value={data.wily.suspension_pole_pair_number_ps}
                    onChange={(e) => updateWinding("suspension_pole_pair_number_ps", parseInt(e.target.value) || 0)}
                  />
                </div>
              </div>
            </TabsContent>

            {/* 几何参数 */}
            <TabsContent value="geometry" className="space-y-4 mt-4">
              <div className="flex items-center justify-between mb-4">
                <div className="text-sm text-muted-foreground">
                  {isLoadingParams ? "正在加载参数信息..." : "紧凑显示所有几何参数的详细信息"}
                </div>
                <Button onClick={loadParametersInfo} variant="outline" size="sm" disabled={isLoadingParams}>
                  <RefreshCw className={`mr-2 h-4 w-4 ${isLoadingParams ? 'animate-spin' : ''}`} />
                  刷新参数信息
                </Button>
              </div>
              <ScrollArea className="h-[700px] pr-4">
                {isLoadingParams ? (
                  <div className="flex items-center justify-center py-12">
                    <Loader2 className="h-6 w-6 animate-spin text-muted-foreground" />
                  </div>
                ) : (
                  <div className="space-y-4">
                    {Object.entries(data.machineGeometry).map(([componentName, componentData]) => {
                      const gpEntries = Object.entries(componentData.GP);
                      return (
                        <Card key={componentName} className="overflow-hidden">
                          <CardHeader className="pb-3">
                            <div className="flex items-center justify-between">
                              <div>
                                <CardTitle className="text-base capitalize">{componentName}</CardTitle>
                                <CardDescription className="text-xs">
                                  {gpEntries.length} 个几何参数
                                </CardDescription>
                              </div>
                              <div className="flex items-center gap-2">
                                <Label className="text-xs">需要重建</Label>
                                <Switch
                                  checked={componentData._needs_rebuild}
                                  onCheckedChange={(checked) => {
                                    setData(prev => ({
                                      ...prev,
                                      machineGeometry: {
                                        ...prev.machineGeometry,
                                        [componentName]: {
                                          ...prev.machineGeometry[componentName as keyof typeof prev.machineGeometry],
                                          _needs_rebuild: checked
                                        }
                                      }
                                    }));
                                  }}
                                  className="scale-75"
                                />
                              </div>
                            </div>
                          </CardHeader>
                          <CardContent className="p-0">
                            <div className="overflow-x-auto">
                              <table className="w-full text-xs border-collapse">
                                <thead>
                                  <tr className="bg-muted/50 border-b">
                                    <th className="text-left p-2 font-semibold">参数键</th>
                                    <th className="text-left p-2 font-semibold">显示名称</th>
                                    <th className="text-center p-2 font-semibold">类型</th>
                                    <th className="text-right p-2 font-semibold">值</th>
                                    <th className="text-right p-2 font-semibold">范围</th>
                                    <th className="text-center p-2 font-semibold">单位</th>
                                    <th className="text-center p-2 font-semibold">代码</th>
                                    <th className="text-center p-2 font-semibold w-12">操作</th>
                                  </tr>
                                </thead>
                                <tbody>
                                  {gpEntries.map(([gpKey, gpDisplayName]) => {
                                    // 从顶层 data 读取参数定义（包含 type, unit, bounds 等）
                                    const paramData = data[gpKey] as ParameterData | undefined;
                                    
                                    // 从组件内部读取实际值（优先级更高，因为这是组件特定的值）
                                    const componentValue = (componentData as any)[gpKey];
                                    
                                    // 从 API 获取的参数信息中查找（使用 name 字段匹配）
                                    let apiParamInfo: ParameterInfo | null = null;
                                    if (paramData?.name) {
                                      apiParamInfo = parametersInfo[paramData.name] || null;
                                    } else {
                                      // 如果顶层没有定义，尝试直接用 gpKey 查找
                                      apiParamInfo = parametersInfo[gpKey] || null;
                                    }
                                    
                                    // 构建完整的参数信息
                                    // 优先使用组件内部的值，如果没有则使用顶层定义的值，最后使用 API 的值
                                    const finalValue = componentValue !== undefined 
                                      ? componentValue 
                                      : (paramData?.value !== undefined ? paramData.value : (apiParamInfo?.value));
                                    
                                    // 优先使用顶层定义的类型和单位，如果没有则使用 API 的信息
                                    const finalType = paramData?.type || apiParamInfo?.type || "unknown";
                                    const finalUnit = paramData?.unit || apiParamInfo?.unit || "-";
                                    const finalBounds = paramData?.bounds || apiParamInfo?.bounds || null;
                                    
                                    const paramInfo: ParameterInfo | null = (paramData || componentValue !== undefined || apiParamInfo) ? {
                                      name: paramData?.name || apiParamInfo?.name || gpKey,
                                      displayName: gpDisplayName,
                                      type: finalType,
                                      unit: finalUnit,
                                      value: finalValue,
                                      bounds: finalBounds,
                                      calc: apiParamInfo?.calc || null,
                                      calc_bounds: apiParamInfo?.calc_bounds || null,
                                      args: apiParamInfo?.args || [],
                                      comment: paramData?.comment || apiParamInfo?.comment || ""
                                    } : null;
                                    
                                    const rowKey = `${componentName}-${gpKey}`;
                                    const isExpanded = expandedGP[rowKey] || false;
                                    const editingType = editingCalc[rowKey] || null;
                                    const isEditingCalc = editingType === "calc";
                                    const isEditingCalcBounds = editingType === "calc_bounds";
                                    
                                    return (
                                      <React.Fragment key={gpKey}>
                                        <tr className="border-b hover:bg-muted/30">
                                          <td className="p-2 font-mono">
                                            <Input
                                              value={gpKey}
                                              onChange={(e) => updateGPKey(componentName, gpKey, e.target.value)}
                                              className="h-7 text-xs font-mono w-32"
                                            />
                                          </td>
                                          <td className="p-2">
                                            <Input
                                              value={gpDisplayName}
                                              onChange={(e) => updateGPDisplayName(componentName, gpKey, e.target.value)}
                                              className="h-7 text-xs w-48"
                                            />
                                          </td>
                                          <td className="p-2 text-center">
                                            {paramInfo ? (
                                              <Badge variant={
                                                paramInfo.type === "fixed" ? "default" :
                                                paramInfo.type === "free" ? "secondary" :
                                                "outline"
                                              } className="text-xs">
                                                {paramInfo.type}
                                              </Badge>
                                            ) : (
                                              <span className="text-muted-foreground">-</span>
                                            )}
                                          </td>
                                          <td className="p-2 text-right font-mono">
                                            {(() => {
                                              // 优先显示 paramInfo 的值，如果没有则显示组件内部的值
                                              const displayValue = paramInfo?.value !== null && paramInfo?.value !== undefined
                                                ? paramInfo.value
                                                : (componentValue !== undefined ? componentValue : null);
                                              
                                              if (displayValue !== null && displayValue !== undefined) {
                                                return typeof displayValue === 'number' 
                                                  ? displayValue.toFixed(3) 
                                                  : String(displayValue);
                                              }
                                              return "-";
                                            })()}
                                          </td>
                                          <td className="p-2 text-right font-mono text-xs">
                                            {paramInfo?.bounds 
                                              ? `[${paramInfo.bounds[0]}, ${paramInfo.bounds[1]}]`
                                              : "-"}
                                          </td>
                                          <td className="p-2 text-center text-muted-foreground">
                                            {paramInfo?.unit || "-"}
                                          </td>
                                          <td className="p-2 text-center">
                                            <div className="flex gap-1 justify-center">
                                              {paramInfo?.calc && (
                                                <Button
                                                  variant="ghost"
                                                  size="sm"
                                                  onClick={() => setEditingCalc(prev => ({
                                                    ...prev,
                                                    [rowKey]: isEditingCalc ? null : "calc"
                                                  }))}
                                                  className="h-6 px-2 text-xs"
                                                  title="查看 calc"
                                                >
                                                  <Code className="h-3 w-3" />
                                                </Button>
                                              )}
                                              {paramInfo?.calc_bounds && (
                                                <Button
                                                  variant="ghost"
                                                  size="sm"
                                                  onClick={() => setEditingCalc(prev => ({
                                                    ...prev,
                                                    [rowKey]: isEditingCalcBounds ? null : "calc_bounds"
                                                  }))}
                                                  className="h-6 px-2 text-xs"
                                                  title="查看 calc_bounds"
                                                >
                                                  <Code className="h-3 w-3" />
                                                </Button>
                                              )}
                                              {!paramInfo?.calc && !paramInfo?.calc_bounds && (
                                                <span className="text-muted-foreground text-xs">-</span>
                                              )}
                                            </div>
                                          </td>
                                          <td className="p-2 text-center">
                                            <Button
                                              variant="ghost"
                                              size="sm"
                                              onClick={() => removeGPFromComponent(componentName, gpKey)}
                                              className="h-6 w-6 p-0"
                                            >
                                              <Trash2 className="h-3 w-3" />
                                            </Button>
                                          </td>
                                        </tr>
                                        {/* 展开显示代码编辑器 */}
                                        {(isEditingCalc || isEditingCalcBounds) && paramInfo && (
                                          <tr>
                                            <td colSpan={8} className="p-0">
                                              {isEditingCalc && paramInfo.calc && (
                                                <div className="border-t-2 border-green-300 dark:border-green-700 bg-slate-50 dark:bg-slate-900">
                                                  <div className="px-3 py-1 text-xs font-semibold bg-green-200 dark:bg-green-800 text-green-900 dark:text-green-100">
                                                    calc 计算方法 (args: {paramInfo.args?.join(", ") || "无"})
                                                  </div>
                                                  <div style={{ height: "120px" }}>
                                                    <Editor
                                                      height="120px"
                                                      defaultLanguage="python"
                                                      value={paramInfo.calc}
                                                      theme="vs-dark"
                                                      options={{
                                                        minimap: { enabled: false },
                                                        fontSize: 11,
                                                        wordWrap: "on",
                                                        automaticLayout: true,
                                                        readOnly: true,
                                                        lineNumbers: "on",
                                                      }}
                                                    />
                                                  </div>
                                                </div>
                                              )}
                                              {isEditingCalcBounds && paramInfo.calc_bounds && (
                                                <div className="border-t-2 border-blue-300 dark:border-blue-700 bg-slate-50 dark:bg-slate-900">
                                                  <div className="px-3 py-1 text-xs font-semibold bg-blue-200 dark:bg-blue-800 text-blue-900 dark:text-blue-100">
                                                    calc_bounds 计算方法 (args: {paramInfo.args?.join(", ") || "无"})
                                                  </div>
                                                  <div style={{ height: "120px" }}>
                                                    <Editor
                                                      height="120px"
                                                      defaultLanguage="python"
                                                      value={paramInfo.calc_bounds}
                                                      theme="vs-dark"
                                                      options={{
                                                        minimap: { enabled: false },
                                                        fontSize: 11,
                                                        wordWrap: "on",
                                                        automaticLayout: true,
                                                        readOnly: true,
                                                        lineNumbers: "on",
                                                      }}
                                                    />
                                                  </div>
                                                </div>
                                              )}
                                            </td>
                                          </tr>
                                        )}
                                      </React.Fragment>
                                    );
                                  })}
                                </tbody>
                              </table>
                            </div>
                            <div className="p-3 border-t">
                              <Button
                                variant="outline"
                                size="sm"
                                onClick={() => {
                                  const gpKey = prompt("输入新的GP参数键:");
                                  const gpDisplayName = prompt("输入显示名称:");
                                  if (gpKey && gpDisplayName) {
                                    addGPToComponent(componentName, gpKey, gpDisplayName);
                                  }
                                }}
                                className="w-full"
                              >
                                <Plus className="mr-2 h-4 w-4" />
                                添加GP参数
                              </Button>
                            </div>
                          </CardContent>
                        </Card>
                      );
                    })}
                  </div>
                )}
              </ScrollArea>
            </TabsContent>

            {/* JSON视图 */}
            <TabsContent value="json" className="mt-4">
              <ScrollArea className="h-[600px] w-full rounded-md border p-4">
                <pre className="text-xs font-mono">
                  {JSON.stringify(data, null, 2)}
                </pre>
              </ScrollArea>
            </TabsContent>
          </Tabs>
        </CardContent>
      </Card>
    </div>
  );
}

