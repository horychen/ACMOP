"use client";

import React from 'react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import LinearMachineView from '@/components/LinearMachineView';
import WindingDiagrams from '@/components/WindingDiagrams';
import { DesignData } from '@/lib/DesignData';
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from '@/components/ui/table';

interface DesignVisualizerClientProps {
  data: DesignData;
}

export default function DesignVisualizerClient({ data }: DesignVisualizerClientProps) {
  // Extract winding layout data
  const wily = data["EX-user"]?.wily || {};
  const layer_X_phases = wily.layer_X_phases || null;
  const layer_X_signs = wily.layer_X_signs || null;
  const layer_Y_phases = wily.layer_Y_phases || null;
  const layer_Y_signs = wily.layer_Y_signs || null;
  const grouping_AC = wily.grouping_AC || null;

  // Extract GP parameters
  const gpParams = data.GP || {};

  // Extract EX-user parameters
  const exParams = data["EX-user"] || {};

  return (
    <div className="space-y-6">
      {/* Machine Geometry Visualization */}
      <Card>
        <CardHeader>
          <CardTitle>机器几何视图</CardTitle>
          <CardDescription>电机横截面可视化</CardDescription>
        </CardHeader>
        <CardContent>
          <LinearMachineView
            Qs={data.Qs}
            p={data.p}
            ps={data.ps}
            coilPitchY={data.coil_pitch_y}
            layer_X_phases={layer_X_phases}
            layer_Y_phases={layer_Y_phases}
            layer_X_signs={layer_X_signs}
            layer_Y_signs={layer_Y_signs}
            grouping_AC={grouping_AC}
          />
        </CardContent>
      </Card>

      {/* Winding Diagrams */}
      {(layer_X_phases || layer_Y_phases) && (
        <Card>
          <CardHeader>
            <CardTitle>绕组布局图</CardTitle>
            <CardDescription>绕组相位分配和连接图</CardDescription>
          </CardHeader>
          <CardContent>
            {layer_X_phases && (
              <WindingDiagrams
                Qs={data.Qs}
                p={data.p}
                m={data.m || 3}
                layer_X_phases={layer_X_phases}
                layer_X_signs={layer_X_signs || []}
              />
            )}
          </CardContent>
        </Card>
      )}

      {/* Parameters Tables */}
      <Tabs defaultValue="geometric" className="w-full">
        <TabsList className="grid w-full grid-cols-3">
          <TabsTrigger value="geometric">几何参数</TabsTrigger>
          <TabsTrigger value="excitation">激励参数</TabsTrigger>
          <TabsTrigger value="performance">性能指标</TabsTrigger>
        </TabsList>

        {/* Geometric Parameters */}
        <TabsContent value="geometric" className="mt-4">
          <Card>
            <CardHeader>
              <CardTitle>几何参数 (GP)</CardTitle>
              <CardDescription>设计几何参数列表</CardDescription>
            </CardHeader>
            <CardContent>
              {Object.keys(gpParams).length > 0 ? (
                <div className="max-h-96 overflow-y-auto">
                  <Table>
                    <TableHeader>
                      <TableRow>
                        <TableHead>参数名称</TableHead>
                        <TableHead>类型</TableHead>
                        <TableHead>数值</TableHead>
                        <TableHead>边界</TableHead>
                      </TableRow>
                    </TableHeader>
                    <TableBody>
                      {Object.entries(gpParams).map(([key, param]: [string, any]) => (
                        <TableRow key={key}>
                          <TableCell className="font-medium">{key}</TableCell>
                          <TableCell>
                            <span className={`px-2 py-1 rounded text-xs ${
                              param.type === 'fixed' ? 'bg-blue-100 text-blue-800 dark:bg-blue-900 dark:text-blue-200' :
                              param.type === 'free' ? 'bg-green-100 text-green-800 dark:bg-green-900 dark:text-green-200' :
                              'bg-gray-100 text-gray-800 dark:bg-gray-900 dark:text-gray-200'
                            }`}>
                              {param.type}
                            </span>
                          </TableCell>
                          <TableCell>{typeof param.value === 'number' ? param.value.toFixed(4) : String(param.value)}</TableCell>
                          <TableCell>
                            {param.bounds 
                              ? `[${param.bounds[0].toFixed(4)}, ${param.bounds[1].toFixed(4)}]`
                              : 'N/A'
                            }
                          </TableCell>
                        </TableRow>
                      ))}
                    </TableBody>
                  </Table>
                </div>
              ) : (
                <div className="text-center py-8 text-muted-foreground">
                  暂无几何参数数据
                </div>
              )}
            </CardContent>
          </Card>
        </TabsContent>

        {/* Excitation Parameters */}
        <TabsContent value="excitation" className="mt-4">
          <Card>
            <CardHeader>
              <CardTitle>激励参数 (EX-user)</CardTitle>
              <CardDescription>激励和运行参数</CardDescription>
            </CardHeader>
            <CardContent>
              {Object.keys(exParams).length > 0 ? (
                <div className="max-h-96 overflow-y-auto">
                  <Table>
                    <TableHeader>
                      <TableRow>
                        <TableHead>参数名称</TableHead>
                        <TableHead>数值</TableHead>
                      </TableRow>
                    </TableHeader>
                    <TableBody>
                      {Object.entries(exParams)
                        .filter(([key]) => key !== 'wily') // Exclude wily as it's displayed separately
                        .map(([key, value]: [string, any]) => (
                          <TableRow key={key}>
                            <TableCell className="font-medium">{key}</TableCell>
                            <TableCell>
                              {typeof value === 'number' 
                                ? value.toFixed(4) 
                                : typeof value === 'object' 
                                  ? JSON.stringify(value, null, 2)
                                  : String(value)
                              }
                            </TableCell>
                          </TableRow>
                        ))}
                    </TableBody>
                  </Table>
                </div>
              ) : (
                <div className="text-center py-8 text-muted-foreground">
                  暂无激励参数数据
                </div>
              )}
            </CardContent>
          </Card>
        </TabsContent>

        {/* Performance Metrics */}
        <TabsContent value="performance" className="mt-4">
          <Card>
            <CardHeader>
              <CardTitle>性能指标</CardTitle>
              <CardDescription>FEA评估的性能指标</CardDescription>
            </CardHeader>
            <CardContent>
              {data["FEA_Evaluated_Performance--1-Initial"] ? (
                <div className="max-h-96 overflow-y-auto">
                  <Table>
                    <TableHeader>
                      <TableRow>
                        <TableHead>指标名称</TableHead>
                        <TableHead>数值</TableHead>
                      </TableRow>
                    </TableHeader>
                    <TableBody>
                      {Object.entries(data["FEA_Evaluated_Performance--1-Initial"])
                        .filter(([key]) => !key.startsWith('_'))
                        .map(([key, value]: [string, any]) => (
                          <TableRow key={key}>
                            <TableCell className="font-medium">{key}</TableCell>
                            <TableCell>
                              {typeof value === 'number' 
                                ? value.toFixed(4) 
                                : String(value)
                              }
                            </TableCell>
                          </TableRow>
                        ))}
                    </TableBody>
                  </Table>
                </div>
              ) : (
                <div className="text-center py-8 text-muted-foreground">
                  暂无性能指标数据
                </div>
              )}
            </CardContent>
          </Card>
        </TabsContent>
      </Tabs>

      {/* Basic Machine Info */}
      <Card>
        <CardHeader>
          <CardTitle>机器基本信息</CardTitle>
          <CardDescription>机器类型和基本规格</CardDescription>
        </CardHeader>
        <CardContent>
          <div className="grid grid-cols-2 gap-4">
            <div>
              <div className="text-sm text-muted-foreground">机器类型</div>
              <div className="text-lg font-semibold">{data.machine_type || 'N/A'}</div>
            </div>
            <div>
              <div className="text-sm text-muted-foreground">相数 (m)</div>
              <div className="text-lg font-semibold">{data.m || 'N/A'}</div>
            </div>
            <div>
              <div className="text-sm text-muted-foreground">定子槽数 (Qs)</div>
              <div className="text-lg font-semibold">{data.Qs || 'N/A'}</div>
            </div>
            <div>
              <div className="text-sm text-muted-foreground">极对数 (p)</div>
              <div className="text-lg font-semibold">{data.p || 'N/A'}</div>
            </div>
            <div>
              <div className="text-sm text-muted-foreground">悬浮极对数 (ps)</div>
              <div className="text-lg font-semibold">{data.ps || 'N/A'}</div>
            </div>
            <div>
              <div className="text-sm text-muted-foreground">机械功率 (kW)</div>
              <div className="text-lg font-semibold">{data.mec_power ? `${data.mec_power} kW` : 'N/A'}</div>
            </div>
          </div>
        </CardContent>
      </Card>
    </div>
  );
}

