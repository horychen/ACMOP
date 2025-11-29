"use client";

import React, { useState, useEffect } from 'react';
import { Activity, Settings, Zap, Cpu, Maximize, RotateCw, Wind, BarChart3, Loader2, Sparkles, AlertCircle, Download } from 'lucide-react';
import { DesignSpecs, OptimizationResult } from '../../types';
import { calculateMachineDesign, generateEfficiencyCurve } from '../../services/physicsEngine';
import { suggestSpecsFromDescription, analyzeDesignResult } from '../../services/geminiService';
import LinearMachineView from '../../components/LinearMachineView';
import { EfficiencyChart } from '../../components/Charts';
import { InputGroup, KpiCard, ResultRow } from '../../components/DesignHelpers';
import axios from 'axios';
import CsvVisualizer from '../../components/CsvVisualizer';
import CsvChartVisualizer from '../../components/CsvChartVisualizer';
import PdfViewer from '../../components/PdfViewer';

const DEFAULT_SPECS: DesignSpecs = {
    ratedPower: 5, // 5 kW
    ratedSpeed: 3000, // 3000 RPM
    ratedVoltage: 400,
    outerDiameterLimit: 120,
    axialLengthLimit: 100,
    airGap: 1.0,
    slotCount: 12,
    poleCount: 4,
    currentDensity: 5
};

export default function VisualizerPage() {
    const [specs, setSpecs] = useState<DesignSpecs>(DEFAULT_SPECS);
    const [result, setResult] = useState<OptimizationResult | null>(null);
    const [isAnalyzing, setIsAnalyzing] = useState(false);
    const [aiAnalysis, setAiAnalysis] = useState<string | null>(null);
    const [userPrompt, setUserPrompt] = useState('');
    const [isSuggesting, setIsSuggesting] = useState(false);

    // Data loading state
    const [availableConfigs, setAvailableConfigs] = useState<Record<string, any>>({});
    const [selectedConfigKey, setSelectedConfigKey] = useState<string>('');
    const [isLoadingConfigs, setIsLoadingConfigs] = useState(false);

    // Metadata from machine_designer_full.json
    const [metadata, setMetadata] = useState<Record<string, any> | null>(null);
    const [isLoadingMetadata, setIsLoadingMetadata] = useState(false);
    const [exData, setExData] = useState<Record<string, any> | null>(null);
    const [path2FEACsv, setPath2FEACsv] = useState<string | null>(null);
    const [pdfUrl, setPdfUrl] = useState<string | null>(null);
    const [coilPitchY, setCoilPitchY] = useState<number | null>(null);
    const [wilyData, setWilyData] = useState<Record<string, any> | null>(null);

    // Fetch available configurations on mount
    useEffect(() => {
        const fetchConfigs = async () => {
            setIsLoadingConfigs(true);
            try {
                const response = await axios.get('/api/machine-specs');
                setAvailableConfigs(response.data);
            } catch (error) {
                console.error("Failed to load machine specs", error);
            } finally {
                setIsLoadingConfigs(false);
            }
        };
        fetchConfigs();
    }, []);

    // Fetch metadata from machine_designer_full.json
    useEffect(() => {
        const fetchMetadata = async () => {
            setIsLoadingMetadata(true);
            try {
                const response = await axios.get('/api/machine-designer/full');
                console.log("Full response data:", response.data);
                
                if (response.data) {
                    // Extract metadata fields (lines 2-10 from JSON)
                    const meta = {
                        name: response.data.name,
                        bool_PermanentMagnet: response.data.bool_PermanentMagnet,
                        bool_StatorSlotClosed: response.data.bool_StatorSlotClosed,
                        bool_RotorNotched: response.data.bool_RotorNotched,
                        select_FEA_tool: response.data.select_FEA_tool,
                        select_fea_config_dict: response.data.select_fea_config_dict,
                        bool_jmagDeleteResultsAfterCalculation: response.data.bool_jmagDeleteResultsAfterCalculation,
                        counter: response.data.counter,
                        parameters: response.data.parameters
                    };
                    setMetadata(meta);
                    
                    // Extract EX data (lines 20-53 from JSON)
                    if (response.data.EX && typeof response.data.EX === 'object') {
                        console.log("EX data found:", response.data.EX);
                        setExData(response.data.EX);
                    } else {
                        console.warn("EX data not found or invalid. Available keys:", Object.keys(response.data || {}));
                        setExData(null);
                    }
                    
                    // Extract path2FEACsv
                    if (response.data.path2FEACsv) {
                        setPath2FEACsv(response.data.path2FEACsv);
                    } else {
                        setPath2FEACsv(null);
                    }
                    
                    // Extract wily data (including coil_pitch_y, Qs, p, ps)
                    if (response.data.wily && typeof response.data.wily === 'object') {
                        setWilyData(response.data.wily);
                        if (response.data.wily.coil_pitch_y !== undefined) {
                            setCoilPitchY(response.data.wily.coil_pitch_y);
                        } else {
                            setCoilPitchY(null);
                        }
                    } else {
                        setWilyData(null);
                        setCoilPitchY(null);
                    }
                    
                    // Load PDF URL
                    const backendUrl = process.env.NEXT_PUBLIC_BACKEND_URL || "http://localhost:8000";
                    setPdfUrl(`${backendUrl}/api/results/pdf/machine-geometry`);
                } else {
                    console.error("Response data is empty");
                    setExData(null);
                    setPath2FEACsv(null);
                }
            } catch (error: any) {
                console.error("Failed to load metadata:", error);
                if (error.response) {
                    console.error("Error response:", error.response.data);
                }
                setExData(null);
            } finally {
                setIsLoadingMetadata(false);
            }
        };
        fetchMetadata();
    }, []);

    // Auto-recalculate when specs change
    useEffect(() => {
        const timer = setTimeout(() => {
            const res = calculateMachineDesign(specs);
            setResult(res);
            setAiAnalysis(null);
        }, 200);
        return () => clearTimeout(timer);
    }, [specs]);

    const handleInputChange = (field: keyof DesignSpecs, value: string) => {
        const numVal = parseFloat(value);
        if (!isNaN(numVal)) {
            setSpecs(prev => ({ ...prev, [field]: numVal }));
        }
    };

    const handleConfigSelect = (key: string) => {
        setSelectedConfigKey(key);
        if (!key || !availableConfigs[key]) return;

        const config = availableConfigs[key];

        const newSpecs: DesignSpecs = {
            ratedPower: (config.mec_power || 5000) / 1000,
            ratedSpeed: config.ExcitationFreqSimulated ? (config.ExcitationFreqSimulated * 60) / (config.p || 1) : 3000,
            ratedVoltage: config.VoltageRating || 400,
            outerDiameterLimit: (config["GP-user"]?.mm_r_so?.value * 2) || 120,
            axialLengthLimit: config["EX-user"]?.mm_stack_length || config.mm_stack_length || 100,
            airGap: config["GP-user"]?.mm_d_mech_air_gap?.value || config.minimum_mechanical_air_gap_length_mm || 1.0,
            slotCount: config.Qs || 12,
            poleCount: (config.p || 2) * 2,
            currentDensity: config["EX-user"]?.Js ? config["EX-user"].Js / 1e6 : 5
        };

        if (!newSpecs.ratedSpeed || isNaN(newSpecs.ratedSpeed)) {
            newSpecs.ratedSpeed = 3000;
        }

        setSpecs(newSpecs);
    };

    const handleAiSuggest = async () => {
        if (!userPrompt.trim()) return;
        setIsSuggesting(true);
        const suggested = await suggestSpecsFromDescription(userPrompt);
        if (suggested && Object.keys(suggested).length > 0) {
            setSpecs(prev => ({ ...prev, ...suggested }));
        }
        setIsSuggesting(false);
    };

    const handleAiAnalyze = async () => {
        if (!result) return;
        setIsAnalyzing(true);
        const text = await analyzeDesignResult(result.specs, result.performance);
        setAiAnalysis(text);
        setIsAnalyzing(false);
    };

    const efficiencyData = result ? generateEfficiencyCurve(result.specs) : [];

    return (
        <div className="min-h-screen flex flex-col font-sans bg-background text-foreground">

            {/* Header */}
            <header className="bg-card border-b border-border p-4 flex items-center justify-between sticky top-0 z-50">
                <div className="flex items-center space-x-3">
                    <div className="bg-primary p-2 rounded-lg">
                        <RotateCw className="text-primary-foreground h-6 w-6" />
                    </div>
                    <div>
                        <h1 className="text-xl font-bold tracking-tight">ACMOP <span className="text-primary font-mono text-sm">Visualizer</span></h1>
                        <p className="text-xs text-muted-foreground">Bearingless Machine Designer Integration</p>
                    </div>
                </div>

                {/* Configuration Loader */}
                <div className="flex items-center space-x-4">
                    <div className="flex items-center space-x-2">
                        <span className="text-xs text-muted-foreground">Load Config:</span>
                        <select
                            value={selectedConfigKey}
                            onChange={(e) => handleConfigSelect(e.target.value)}
                            className="bg-card border border-border text-xs rounded px-2 py-1 text-foreground focus:ring-1 focus:ring-ring outline-none max-w-[200px]"
                            disabled={isLoadingConfigs}
                        >
                            <option value="">-- Select Machine --</option>
                            {Object.keys(availableConfigs).map(key => (
                                <option key={key} value={key}>{key}</option>
                            ))}
                        </select>
                        {isLoadingConfigs && <Loader2 className="w-3 h-3 animate-spin text-muted-foreground" />}
                    </div>
                </div>
            </header>

            <main className="flex-1 flex flex-col lg:flex-row overflow-hidden h-[calc(100vh-73px)]">

                {/* Left Sidebar: Controls & AI Input */}
                <aside className="w-full lg:w-96 bg-card border-r border-border flex flex-col h-full overflow-y-auto">

                    {/* AI Design Assistant */}
                    <div className="p-4 border-b border-border bg-muted/50">
                        <h2 className="text-sm font-semibold text-primary mb-3 flex items-center">
                            <Sparkles className="w-4 h-4 mr-2" /> Design Assistant
                        </h2>
                        <div className="space-y-2">
                            <textarea
                                value={userPrompt}
                                onChange={(e) => setUserPrompt(e.target.value)}
                                placeholder="e.g., I need a high-speed motor for a drone, 500W, 12V..."
                                className="w-full bg-background border border-border rounded-md p-3 text-sm text-foreground focus:ring-1 focus:ring-ring outline-none resize-none h-20 placeholder-muted-foreground"
                            />
                            <button
                                onClick={handleAiSuggest}
                                disabled={isSuggesting || !userPrompt}
                                className="w-full bg-primary hover:bg-primary/90 disabled:opacity-50 disabled:cursor-not-allowed text-white text-xs font-semibold py-2 px-4 rounded-md transition-colors flex items-center justify-center"
                            >
                                {isSuggesting ? <Loader2 className="animate-spin w-4 h-4 mr-2" /> : "Generate Specs"}
                            </button>
                        </div>
                    </div>

                    {/* Manual Parameter Inputs */}
                    <div className="p-4 space-y-6">
                        <div className="flex items-center justify-between">
                            <h2 className="text-sm font-semibold text-foreground flex items-center">
                                <Settings className="w-4 h-4 mr-2" /> Parameters
                            </h2>
                        </div>

                        {/* Power & Speed */}
                        <div className="space-y-4">
                            <h3 className="text-xs uppercase tracking-wider text-muted-foreground font-bold">Ratings</h3>
                            <div className="grid grid-cols-2 gap-3">
                                <InputGroup label="Power (kW)" value={specs.ratedPower} onChange={v => handleInputChange('ratedPower', v)} icon={<Zap className="w-3 h-3" />} />
                                <InputGroup label="Speed (RPM)" value={specs.ratedSpeed} onChange={v => handleInputChange('ratedSpeed', v)} icon={<Activity className="w-3 h-3" />} />
                                <InputGroup label="Voltage (V)" value={specs.ratedVoltage} onChange={v => handleInputChange('ratedVoltage', v)} />
                                <InputGroup label="Current Dens. (A/mm²)" value={specs.currentDensity} onChange={v => handleInputChange('currentDensity', v)} />
                            </div>
                        </div>

                        {/* Geometry Limits */}
                        <div className="space-y-4">
                            <h3 className="text-xs uppercase tracking-wider text-muted-foreground font-bold">Constraints</h3>
                            <div className="grid grid-cols-2 gap-3">
                                <InputGroup label="Outer Dia. (mm)" value={specs.outerDiameterLimit} onChange={v => handleInputChange('outerDiameterLimit', v)} icon={<Maximize className="w-3 h-3" />} />
                                <InputGroup label="Length (mm)" value={specs.axialLengthLimit} onChange={v => handleInputChange('axialLengthLimit', v)} />
                                <InputGroup label="Air Gap (mm)" value={specs.airGap} step={0.1} onChange={v => handleInputChange('airGap', v)} />
                            </div>
                        </div>

                        {/* Topology */}
                        <div className="space-y-4">
                            <h3 className="text-xs uppercase tracking-wider text-muted-foreground font-bold">Topology</h3>
                            <div className="grid grid-cols-2 gap-3">
                                <InputGroup label="Slot Count" value={specs.slotCount} step={3} onChange={v => handleInputChange('slotCount', v)} />
                                <InputGroup label="Pole Count" value={specs.poleCount} step={2} onChange={v => handleInputChange('poleCount', v)} />
                            </div>
                        </div>

                        {/* Alternative CSV Chart Visualizer */}
                        {path2FEACsv && (
                            <div className="p-4 border-t border-border mt-6">
                                <h2 className="text-sm font-semibold text-foreground mb-3 flex items-center">
                                    <BarChart3 className="w-4 h-4 mr-2" /> CSV 图表可视化
                                </h2>
                                <div className="h-[400px]">
                                    <CsvChartVisualizer 
                                        path2FEACsv={path2FEACsv}
                                        projectName={metadata?.name}
                                    />
                                </div>
                            </div>
                        )}
                    </div>
                </aside>

                {/* Center: Visualization & Results */}
                <div className="flex-1 bg-background flex flex-col h-full overflow-hidden">

                    <div className="flex-1 overflow-y-auto p-6">
                        {/* Visualization Panel - Full Width */}
                        <div className="mb-8">
                            <div className="flex items-center justify-between mb-4">
                                <h3 className="text-lg font-medium text-foreground">Geometry</h3>
                            </div>
                            
                            <div className="grid grid-cols-1 lg:grid-cols-2 gap-6 mb-6">
                                {/* PDF Cross Section Viewer */}
                                <div className="h-[400px] bg-card rounded-lg border border-border overflow-hidden">
                                    {pdfUrl ? (
                                        <PdfViewer pdfUrl={pdfUrl} />
                                    ) : isLoadingMetadata ? (
                                        <div className="h-full flex items-center justify-center">
                                            <div className="text-center">
                                                <Loader2 className="w-6 h-6 animate-spin text-primary mx-auto mb-2" />
                                                <p className="text-sm text-muted-foreground">加载横截面 PDF...</p>
                                            </div>
                                        </div>
                                    ) : (
                                        <div className="h-full flex items-center justify-center">
                                            <div className="text-center p-6">
                                                <AlertCircle className="w-8 h-8 text-muted-foreground mx-auto mb-2" />
                                                <p className="text-sm text-muted-foreground">PDF 路径未加载</p>
                                            </div>
                                        </div>
                                    )}
                                </div>

                                {/* Stats Summary */}
                                <div className="flex flex-col justify-center">
                                    <div className="grid grid-cols-3 gap-2 text-center">
                                        <div className="bg-card p-3 rounded border border-border">
                                            <div className="text-xs text-muted-foreground">Stator OD</div>
                                            <div className="text-primary font-mono font-bold">{(result?.geometry.statorOuterRadius! * 2).toFixed(1)} mm</div>
                                        </div>
                                        <div className="bg-card p-3 rounded border border-border">
                                            <div className="text-xs text-muted-foreground">Rotor OD</div>
                                            <div className="text-primary font-mono font-bold">{(result?.geometry.rotorOuterRadius! * 2).toFixed(1)} mm</div>
                                        </div>
                                        <div className="bg-card p-3 rounded border border-border">
                                            <div className="text-xs text-muted-foreground">Slot Fill</div>
                                            <div className="text-primary font-mono font-bold">45%</div>
                                        </div>
                                    </div>
                                </div>
                            </div>

                            {/* Linear Geometry View - Full Width with More Space */}
                            {wilyData && (
                                <div className="w-full mb-6">
                                    <h4 className="text-sm font-medium text-foreground mb-4">
                                        Linear View (Auto-scaled)
                                        {coilPitchY !== null && (
                                            <span className="ml-2 text-muted-foreground font-normal">
                                                - coil_pitch_y = {coilPitchY}
                                            </span>
                                        )}
                                    </h4>
                                    <div className="w-full bg-card rounded-lg border border-border p-4" style={{ minHeight: '700px' }}>
                                        <LinearMachineView 
                                            Qs={wilyData.stator_slot_number_Qs || wilyData.Qs || 12}
                                            p={wilyData.pole_pair_number_p || wilyData.p || 2}
                                            ps={wilyData.suspension_pole_pair_number_ps || wilyData.ps || 3}
                                            coilPitchY={coilPitchY} 
                                        />
                                    </div>
                                </div>
                            )}
                        </div>

                        {/* Metadata Section */}
                        <div className="mb-8">
                            <h3 className="text-lg font-medium text-foreground mb-4 flex items-center">
                                <Settings className="w-5 h-5 mr-2" /> 元数据信息
                            </h3>
                            {isLoadingMetadata ? (
                                <div className="bg-card rounded-lg border border-border p-6 flex items-center justify-center">
                                    <Loader2 className="w-5 h-5 animate-spin text-muted-foreground mr-2" />
                                    <span className="text-sm text-muted-foreground">加载元数据中...</span>
                                </div>
                            ) : metadata ? (
                                <div className="bg-card rounded-lg border border-border overflow-hidden">
                                    <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4 p-6">
                                        <div className="space-y-1">
                                            <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">名称</div>
                                            <div className="text-sm font-mono text-foreground">{metadata.name || 'N/A'}</div>
                                        </div>
                                        <div className="space-y-1">
                                            <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">永磁体</div>
                                            <div className="text-sm text-foreground">
                                                {metadata.bool_PermanentMagnet ? (
                                                    <span className="text-emerald-500">是</span>
                                                ) : (
                                                    <span className="text-red-500">否</span>
                                                )}
                                            </div>
                                        </div>
                                        <div className="space-y-1">
                                            <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">定子槽封闭</div>
                                            <div className="text-sm text-foreground">
                                                {metadata.bool_StatorSlotClosed ? (
                                                    <span className="text-emerald-500">是</span>
                                                ) : (
                                                    <span className="text-red-500">否</span>
                                                )}
                                            </div>
                                        </div>
                                        <div className="space-y-1">
                                            <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">转子开槽</div>
                                            <div className="text-sm text-foreground">
                                                {metadata.bool_RotorNotched ? (
                                                    <span className="text-emerald-500">是</span>
                                                ) : (
                                                    <span className="text-red-500">否</span>
                                                )}
                                            </div>
                                        </div>
                                        <div className="space-y-1">
                                            <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">FEA 工具</div>
                                            <div className="text-sm font-mono text-foreground">{metadata.select_FEA_tool || 'N/A'}</div>
                                        </div>
                                        <div className="space-y-1">
                                            <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">FEA 配置</div>
                                            <div className="text-sm text-foreground break-words">{metadata.select_fea_config_dict || 'N/A'}</div>
                                        </div>
                                        <div className="space-y-1">
                                            <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">计算后删除结果</div>
                                            <div className="text-sm text-foreground">
                                                {metadata.bool_jmagDeleteResultsAfterCalculation ? (
                                                    <span className="text-emerald-500">是</span>
                                                ) : (
                                                    <span className="text-red-500">否</span>
                                                )}
                                            </div>
                                        </div>
                                        <div className="space-y-1">
                                            <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">计数器</div>
                                            <div className="text-sm font-mono text-foreground">{metadata.counter ?? 'N/A'}</div>
                                        </div>
                                        <div className="space-y-1">
                                            <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">参数数量</div>
                                            <div className="text-sm font-mono text-foreground">
                                                {metadata.parameters && typeof metadata.parameters === 'object' 
                                                    ? Object.keys(metadata.parameters).length 
                                                    : 0}
                                            </div>
                                        </div>
                                    </div>
                                </div>
                            ) : (
                                <div className="bg-card rounded-lg border border-border p-6">
                                    <div className="flex items-center text-muted-foreground">
                                        <AlertCircle className="w-4 h-4 mr-2" />
                                        <span className="text-sm">无法加载元数据</span>
                                    </div>
                                </div>
                            )}
                        </div>

                        {/* EX Parameters Section */}
                        <div className="mb-8">
                            <h3 className="text-lg font-medium text-foreground mb-4 flex items-center">
                                <Zap className="w-5 h-5 mr-2" /> 激励参数 (EX - Excitation)
                            </h3>
                            {isLoadingMetadata ? (
                                <div className="bg-card rounded-lg border border-border p-6 flex items-center justify-center">
                                    <Loader2 className="w-5 h-5 animate-spin text-muted-foreground mr-2" />
                                    <span className="text-sm text-muted-foreground">加载参数中...</span>
                                </div>
                            ) : exData && typeof exData === 'object' && Object.keys(exData).length > 0 ? (
                                <div className="space-y-6">
                                    {/* 几何参数 */}
                                    <div className="bg-card rounded-lg border border-border overflow-hidden">
                                        <div className="bg-muted/50 px-6 py-3 border-b border-border">
                                            <h4 className="text-sm font-semibold text-foreground flex items-center">
                                                <Maximize className="w-4 h-4 mr-2" /> 几何参数
                                            </h4>
                                        </div>
                                        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4 p-6">
                                            <div className="space-y-1">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">叠片长度</div>
                                                <div className="text-sm font-mono text-foreground">{exData.mm_stack_length_specified?.toFixed(2) || 'N/A'} <span className="text-muted-foreground">mm</span></div>
                                            </div>
                                            <div className="space-y-1">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">槽面积</div>
                                                <div className="text-sm font-mono text-foreground">{exData.mm2_slot_area?.toFixed(2) || 'N/A'} <span className="text-muted-foreground">mm²</span></div>
                                            </div>
                                            <div className="space-y-1">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">磁体面积</div>
                                                <div className="text-sm font-mono text-foreground">{exData.mm2_magnet_area?.toFixed(2) || 'N/A'} <span className="text-muted-foreground">mm²</span></div>
                                            </div>
                                            <div className="space-y-1">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">磁体起始角度</div>
                                                <div className="text-sm font-mono text-foreground">{exData.Magnet_StartAngle?.toFixed(2) || 'N/A'} <span className="text-muted-foreground">°</span></div>
                                            </div>
                                            <div className="space-y-1">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">初始旋转角度</div>
                                                <div className="text-sm font-mono text-foreground">{exData.InitialRotationAngle?.toFixed(2) || 'N/A'} <span className="text-muted-foreground">°</span></div>
                                            </div>
                                        </div>
                                    </div>

                                    {/* 材料参数 */}
                                    <div className="bg-card rounded-lg border border-border overflow-hidden">
                                        <div className="bg-muted/50 px-6 py-3 border-b border-border">
                                            <h4 className="text-sm font-semibold text-foreground flex items-center">
                                                <Cpu className="w-4 h-4 mr-2" /> 材料参数
                                            </h4>
                                        </div>
                                        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4 p-6">
                                            <div className="space-y-1">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">磁体材料</div>
                                                <div className="text-sm font-mono text-foreground">{exData.Magnet_Name || 'N/A'}</div>
                                            </div>
                                            <div className="space-y-1">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">钢材材料</div>
                                                <div className="text-sm font-mono text-foreground">{exData.SteelMaterial || 'N/A'}</div>
                                            </div>
                                            <div className="space-y-1">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">定子铁心材料</div>
                                                <div className="text-sm font-mono text-foreground">{exData.StatorCore_Material || 'N/A'}</div>
                                            </div>
                                            <div className="space-y-1">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">转子铁心材料</div>
                                                <div className="text-sm font-mono text-foreground">{exData.RotorCore_Material || 'N/A'}</div>
                                            </div>
                                            <div className="space-y-1">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">叠片系数</div>
                                                <div className="text-sm font-mono text-foreground">{exData.LaminationFactor?.toFixed(1) || 'N/A'} <span className="text-muted-foreground">%</span></div>
                                            </div>
                                        </div>
                                    </div>

                                    {/* 电气参数 */}
                                    <div className="bg-card rounded-lg border border-border overflow-hidden">
                                        <div className="bg-muted/50 px-6 py-3 border-b border-border">
                                            <h4 className="text-sm font-semibold text-foreground flex items-center">
                                                <Zap className="w-4 h-4 mr-2" /> 电气参数
                                            </h4>
                                        </div>
                                        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4 p-6">
                                            <div className="space-y-1">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">额定功率</div>
                                                <div className="text-sm font-mono text-foreground">{(exData.RatedPower / 1000)?.toFixed(1) || 'N/A'} <span className="text-muted-foreground">kW</span></div>
                                            </div>
                                            <div className="space-y-1">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">额定转速</div>
                                                <div className="text-sm font-mono text-foreground">{exData.RatedSpeed?.toLocaleString() || 'N/A'} <span className="text-muted-foreground">RPM</span></div>
                                            </div>
                                            <div className="space-y-1">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">激励频率</div>
                                                <div className="text-sm font-mono text-foreground">{exData.ExcitationFreqSimulated?.toFixed(1) || 'N/A'} <span className="text-muted-foreground">Hz</span></div>
                                            </div>
                                            <div className="space-y-1">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">直流母线电压</div>
                                                <div className="text-sm font-mono text-foreground">{exData.DCBusVoltage?.toFixed(0) || 'N/A'} <span className="text-muted-foreground">V</span></div>
                                            </div>
                                            <div className="space-y-1">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">电流密度</div>
                                                <div className="text-sm font-mono text-foreground">{(exData.Js / 1e6)?.toFixed(1) || 'N/A'} <span className="text-muted-foreground">A/mm²</span></div>
                                            </div>
                                            <div className="space-y-1">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">驱动绕组电阻</div>
                                                <div className="text-sm font-mono text-foreground">{exData.DriveW_Rs?.toFixed(2) || 'N/A'} <span className="text-muted-foreground">Ω</span></div>
                                            </div>
                                            <div className="space-y-1">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">悬浮绕组电阻</div>
                                                <div className="text-sm font-mono text-foreground">{exData.BeariW_Rs?.toFixed(2) || 'N/A'} <span className="text-muted-foreground">Ω</span></div>
                                            </div>
                                            <div className="space-y-1">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">连接方式</div>
                                                <div className="text-sm text-foreground">
                                                    {exData.bool_WyeConnectOrDeltaConnect ? (
                                                        <span className="text-emerald-500 font-mono">Y型连接</span>
                                                    ) : (
                                                        <span className="text-blue-500 font-mono">Δ型连接</span>
                                                    )}
                                                </div>
                                            </div>
                                        </div>
                                    </div>

                                    {/* 温度参数 */}
                                    <div className="bg-card rounded-lg border border-border overflow-hidden">
                                        <div className="bg-muted/50 px-6 py-3 border-b border-border">
                                            <h4 className="text-sm font-semibold text-foreground flex items-center">
                                                <Activity className="w-4 h-4 mr-2" /> 温度参数
                                            </h4>
                                        </div>
                                        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4 p-6">
                                            <div className="space-y-1">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">工作温度</div>
                                                <div className="text-sm font-mono text-foreground">{exData.Temperature?.toFixed(0) || 'N/A'} <span className="text-muted-foreground">°C</span></div>
                                            </div>
                                            <div className="space-y-1">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">磁体温度</div>
                                                <div className="text-sm font-mono text-foreground">{exData.Magnet_Temperature?.toFixed(0) || 'N/A'} <span className="text-muted-foreground">°C</span></div>
                                            </div>
                                        </div>
                                    </div>

                                    {/* 绕组参数 */}
                                    <div className="bg-card rounded-lg border border-border overflow-hidden">
                                        <div className="bg-muted/50 px-6 py-3 border-b border-border">
                                            <h4 className="text-sm font-semibold text-foreground flex items-center">
                                                <Wind className="w-4 h-4 mr-2" /> 绕组参数
                                            </h4>
                                        </div>
                                        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4 p-6">
                                            <div className="space-y-1">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">串联线圈匝数</div>
                                                <div className="text-sm font-mono text-foreground">{exData.no_series_coil_turns_N?.toFixed(0) || 'N/A'}</div>
                                            </div>
                                            <div className="space-y-1">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">驱动绕组每槽导体数</div>
                                                <div className="text-sm font-mono text-foreground">{exData.DriveW_zQ?.toFixed(0) || 'N/A'}</div>
                                            </div>
                                            <div className="space-y-1">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">悬浮绕组每槽导体数</div>
                                                <div className="text-sm font-mono text-foreground">{exData.BeariW_zQ?.toFixed(0) || 'N/A'}</div>
                                            </div>
                                            <div className="space-y-1">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">绕组填充因子</div>
                                                <div className="text-sm font-mono text-foreground">{(exData.WindingFill * 100)?.toFixed(2) || 'N/A'} <span className="text-muted-foreground">%</span></div>
                                            </div>
                                            <div className="space-y-1">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">转矩电流比例</div>
                                                <div className="text-sm font-mono text-foreground">{(exData.TORQUE_CURRENT_RATIO * 100)?.toFixed(1) || 'N/A'} <span className="text-muted-foreground">%</span></div>
                                            </div>
                                            <div className="space-y-1">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">悬浮电流比例</div>
                                                <div className="text-sm font-mono text-foreground">{(exData.SUSPENSION_CURRENT_RATIO * 100)?.toFixed(1) || 'N/A'} <span className="text-muted-foreground">%</span></div>
                                            </div>
                                            <div className="space-y-1">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">槽电流利用率（转矩）</div>
                                                <div className="text-sm font-mono text-foreground">{(exData.slot_current_utilizing_ratio_for_torque * 100)?.toFixed(2) || 'N/A'} <span className="text-muted-foreground">%</span></div>
                                            </div>
                                        </div>
                                    </div>

                                    {/* 电流参数 */}
                                    <div className="bg-card rounded-lg border border-border overflow-hidden">
                                        <div className="bg-muted/50 px-6 py-3 border-b border-border">
                                            <h4 className="text-sm font-semibold text-foreground flex items-center">
                                                <BarChart3 className="w-4 h-4 mr-2" /> 电流参数
                                            </h4>
                                        </div>
                                        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4 p-6">
                                            <div className="space-y-1">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">槽内电流</div>
                                                <div className="text-sm font-mono text-foreground">{exData.CurrentAmp_in_the_slot?.toFixed(2) || 'N/A'} <span className="text-muted-foreground">A</span></div>
                                            </div>
                                            <div className="space-y-1">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">每导体电流</div>
                                                <div className="text-sm font-mono text-foreground">{exData.CurrentAmp_per_conductor?.toFixed(2) || 'N/A'} <span className="text-muted-foreground">A</span></div>
                                            </div>
                                            <div className="space-y-1">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">每相电流</div>
                                                <div className="text-sm font-mono text-foreground">{exData.CurrentAmp_per_phase?.toFixed(2) || 'N/A'} <span className="text-muted-foreground">A</span></div>
                                            </div>
                                            <div className="space-y-1">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">驱动绕组电流</div>
                                                <div className="text-sm font-mono text-foreground">{exData.DriveW_CurrentAmp?.toFixed(2) || 'N/A'} <span className="text-muted-foreground">A</span></div>
                                            </div>
                                            <div className="space-y-1">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">悬浮绕组电流</div>
                                                <div className="text-sm font-mono text-foreground">{exData.BeariW_CurrentAmp?.toFixed(2) || 'N/A'} <span className="text-muted-foreground">A</span></div>
                                            </div>
                                        </div>
                                    </div>
                                </div>
                            ) : (
                                <div className="bg-card rounded-lg border border-border p-6">
                                    <div className="flex items-center text-muted-foreground">
                                        <AlertCircle className="w-4 h-4 mr-2" />
                                        <span className="text-sm">无法加载激励参数 (EX数据为空或未找到)</span>
                                    </div>
                                </div>
                            )}
                        </div>

                        {/* CSV Visualizer Section */}
                        <div className="mb-8">
                            <h3 className="text-lg font-medium text-foreground mb-4">Simulation Results (CSV)</h3>
                            <div className="h-[500px]">
                                <CsvVisualizer 
                                    path2FEACsv={path2FEACsv || undefined}
                                />
                            </div>
                        </div>

                        {/* Performance Panel - Moved to Bottom */}
                        <div className="mb-8">
                            <div className="flex items-center justify-between mb-4">
                                <h3 className="text-lg font-medium text-foreground">Performance</h3>
                                <button
                                    onClick={handleAiAnalyze}
                                    disabled={isAnalyzing}
                                    className="text-xs flex items-center text-primary hover:text-primary/80 transition-colors"
                                >
                                    {isAnalyzing ? "Analyzing..." : "AI Analysis"} <Sparkles className="w-3 h-3 ml-1" />
                                </button>
                            </div>

                            <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
                                {/* KPI Cards */}
                                <div>
                                    <div className="grid grid-cols-2 gap-4 mb-6">
                                        <KpiCard label="Efficiency" value={result?.performance.efficiency.toFixed(1)} unit="%" icon={<Wind className="w-4 h-4 text-emerald-500" />} />
                                        <KpiCard label="Rated Torque" value={result?.performance.torque.toFixed(2)} unit="Nm" icon={<RotateCw className="w-4 h-4 text-amber-500" />} />
                                        <KpiCard label="Suspension Force" value={result?.performance.suspensionForce.toFixed(1)} unit="N" icon={<Activity className="w-4 h-4 text-purple-500" />} />
                                        <KpiCard label="Cost Est." value={result?.performance.materialCost.toFixed(2)} unit="$" icon={<BarChart3 className="w-4 h-4 text-blue-500" />} />
                                    </div>

                                    {/* AI Analysis Box */}
                                    {aiAnalysis && (
                                        <div className="bg-primary/10/20 border border-primary/30 p-4 rounded-lg animate-in fade-in slide-in-from-bottom-2 duration-500">
                                            <h4 className="text-primary text-xs font-bold uppercase mb-2 flex items-center">
                                                <Cpu className="w-3 h-3 mr-1" /> Engineer's Note
                                            </h4>
                                            <p className="text-sm text-primary-foreground leading-relaxed">
                                                {aiAnalysis}
                                            </p>
                                        </div>
                                    )}
                                </div>

                                {/* Charts */}
                                <div>
                                    <EfficiencyChart data={efficiencyData} />
                                </div>
                            </div>
                        </div>

                        {/* Detailed Table */}
                        <div className="mt-8">
                            <h3 className="text-lg font-medium text-foreground mb-4">Detailed Calculation Results</h3>
                            <div className="bg-card rounded-lg border border-border overflow-hidden">
                                <table className="w-full text-sm text-left text-muted-foreground">
                                    <thead className="text-xs text-muted-foreground uppercase bg-background/50">
                                        <tr>
                                            <th className="px-6 py-3">Parameter</th>
                                            <th className="px-6 py-3">Value</th>
                                            <th className="px-6 py-3">Unit</th>
                                        </tr>
                                    </thead>
                                    <tbody className="divide-y divide-slate-800">
                                        <ResultRow label="Copper Loss" value={result?.performance.copperLoss.toFixed(1)} unit="W" />
                                        <ResultRow label="Iron Loss" value={result?.performance.ironLoss.toFixed(1)} unit="W" />
                                        <ResultRow label="Torque Ripple" value={result?.performance.torqueRipple.toFixed(1)} unit="%" />
                                        <ResultRow label="Power Factor" value={result?.performance.powerFactor.toFixed(2)} unit="-" />
                                        <ResultRow label="Tooth Width" value={result?.geometry.toothWidth.toFixed(2)} unit="mm" />
                                        <ResultRow label="Slot Depth" value={result?.geometry.slotDepth.toFixed(2)} unit="mm" />
                                    </tbody>
                                </table>
                            </div>
                        </div>

                    </div>
                </div>
            </main>
        </div>
    );
}
