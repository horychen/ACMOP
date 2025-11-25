"use client"

import { useState, useEffect } from "react";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Switch } from "@/components/ui/switch";
import { Button } from "@/components/ui/button";
import { Separator } from "@/components/ui/separator";
import { ScrollArea } from "@/components/ui/scroll-area";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { Save, Upload, Download, RefreshCw, Plus, Trash2 } from "lucide-react";
import { Checkbox } from "@/components/ui/checkbox";

interface MachineDesignerData {
  machine_class: string;
  bool_PermanentMagnet: boolean;
  bool_StatorSlotClosed: boolean;
  bool_RotorNotched: boolean;
  parameters: Record<string, any>;
  winding: {
    phase_number_m: number;
    stator_slot_number_Qs: number;
    pole_pair_number_p: number;
    suspension_pole_pair_number_ps: number;
  };
  machineGeometry: {
    rotorCore: {
      color: string | null;
      GP: string[];
      _needs_rebuild: boolean;
    };
    shaft: {
      color: string | null;
      GP: string[];
      _needs_rebuild: boolean;
    };
    rotorMagnet: {
      color: string | null;
      GP: string[];
      _needs_rebuild: boolean;
    };
    statorCore: {
      color: string | null;
      GP: string[];
      _needs_rebuild: boolean;
    };
    coils: {
      color: string | null;
      GP: string[];
      _needs_rebuild: boolean;
    };
  };
}

const DEFAULT_DATA: MachineDesignerData = {
  machine_class: "bearingless_spmsm_heart.bearingless_spmsm_design_variant",
  bool_PermanentMagnet: true,
  bool_StatorSlotClosed: false,
  bool_RotorNotched: true,
  parameters: {},
  winding: {
    phase_number_m: 3,
    stator_slot_number_Qs: 12,
    pole_pair_number_p: 4,
    suspension_pole_pair_number_ps: 5
  },
  machineGeometry: {
    rotorCore: {
      color: null,
      GP: ["mm_r_ro", "mm_d_ri", "mm_d_pm", "mm_d_rp", "mm_d_rs", "p", "s"],
      _needs_rebuild: true
    },
    shaft: {
      color: null,
      GP: ["mm_r_ri"],
      _needs_rebuild: true
    },
    rotorMagnet: {
      color: null,
      GP: ["mm_d_pm", "mm_d_ri", "mm_r_ri"],
      _needs_rebuild: true
    },
    statorCore: {
      color: null,
      GP: ["mm_r_si", "mm_d_sto", "mm_d_sts", "mm_d_st", "mm_d_sy", "mm_w_st", "deg_alpha_st", "deg_alpha_sto", "Q"],
      _needs_rebuild: true
    },
    coils: {
      color: null,
      GP: ["mm_r_so", "mm_d_sy", "mm_w_st", "mm_d_st"],
      _needs_rebuild: true
    }
  }
};

export default function MachineDesignerEditor() {
  const [data, setData] = useState<MachineDesignerData>(DEFAULT_DATA);
  const [activeTab, setActiveTab] = useState<string>("general");

  useEffect(() => {
    // 可以在这里加载保存的数据
    loadSavedData();
  }, []);

  const loadSavedData = () => {
    try {
      const saved = localStorage.getItem("machine_designer_data");
      if (saved) {
        setData(JSON.parse(saved));
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
      winding: {
        ...prev.winding,
        [field]: value
      }
    }));
  };

  const updateGeometryGP = (component: string, gpList: string[]) => {
    setData(prev => ({
      ...prev,
      machineGeometry: {
        ...prev.machineGeometry,
        [component]: {
          ...prev.machineGeometry[component as keyof typeof prev.machineGeometry],
          GP: gpList
        }
      }
    }));
  };

  const addGPToComponent = (component: string, gpName: string) => {
    const currentGP = data.machineGeometry[component as keyof typeof data.machineGeometry].GP;
    if (!currentGP.includes(gpName)) {
      updateGeometryGP(component, [...currentGP, gpName]);
    }
  };

  const removeGPFromComponent = (component: string, index: number) => {
    const currentGP = data.machineGeometry[component as keyof typeof data.machineGeometry].GP;
    updateGeometryGP(component, currentGP.filter((_, i) => i !== index));
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
                    value={data.winding.phase_number_m}
                    onChange={(e) => updateWinding("phase_number_m", parseInt(e.target.value) || 0)}
                  />
                </div>
                <div className="space-y-2">
                  <Label htmlFor="stator_slot_number_Qs">定子槽数 (stator_slot_number_Qs)</Label>
                  <Input
                    id="stator_slot_number_Qs"
                    type="number"
                    value={data.winding.stator_slot_number_Qs}
                    onChange={(e) => updateWinding("stator_slot_number_Qs", parseInt(e.target.value) || 0)}
                  />
                </div>
                <div className="space-y-2">
                  <Label htmlFor="pole_pair_number_p">极对数 (pole_pair_number_p)</Label>
                  <Input
                    id="pole_pair_number_p"
                    type="number"
                    value={data.winding.pole_pair_number_p}
                    onChange={(e) => updateWinding("pole_pair_number_p", parseInt(e.target.value) || 0)}
                  />
                </div>
                <div className="space-y-2">
                  <Label htmlFor="suspension_pole_pair_number_ps">悬浮极对数 (suspension_pole_pair_number_ps)</Label>
                  <Input
                    id="suspension_pole_pair_number_ps"
                    type="number"
                    value={data.winding.suspension_pole_pair_number_ps}
                    onChange={(e) => updateWinding("suspension_pole_pair_number_ps", parseInt(e.target.value) || 0)}
                  />
                </div>
              </div>
            </TabsContent>

            {/* 几何参数 */}
            <TabsContent value="geometry" className="space-y-4 mt-4">
              <ScrollArea className="h-[600px] pr-4">
                <div className="space-y-6">
                  {Object.entries(data.machineGeometry).map(([componentName, componentData]) => (
                    <Card key={componentName}>
                      <CardHeader>
                        <CardTitle className="text-lg capitalize">{componentName}</CardTitle>
                        <CardDescription>
                          几何参数 (GP) 列表
                        </CardDescription>
                      </CardHeader>
                      <CardContent className="space-y-4">
                        <div className="space-y-2">
                          <Label>几何参数 (GP)</Label>
                          <div className="space-y-2">
                            {componentData.GP.map((gp, index) => (
                              <div key={index} className="flex items-center gap-2">
                                <Input
                                  value={gp}
                                  onChange={(e) => {
                                    const newGP = [...componentData.GP];
                                    newGP[index] = e.target.value;
                                    updateGeometryGP(componentName, newGP);
                                  }}
                                  className="flex-1"
                                />
                                <Button
                                  variant="ghost"
                                  size="sm"
                                  onClick={() => removeGPFromComponent(componentName, index)}
                                >
                                  <Trash2 className="h-4 w-4" />
                                </Button>
                              </div>
                            ))}
                            <Button
                              variant="outline"
                              size="sm"
                              onClick={() => {
                                const newGP = prompt("输入新的GP参数名:");
                                if (newGP) {
                                  addGPToComponent(componentName, newGP);
                                }
                              }}
                              className="w-full"
                            >
                              <Plus className="mr-2 h-4 w-4" />
                              添加GP参数
                            </Button>
                          </div>
                        </div>
                        <div className="flex items-center justify-between">
                          <Label>需要重建 (_needs_rebuild)</Label>
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
                          />
                        </div>
                      </CardContent>
                    </Card>
                  ))}
                </div>
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

