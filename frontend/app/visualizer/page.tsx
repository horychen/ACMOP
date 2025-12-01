"use client";

import React, { useState, useEffect } from 'react';
import { Activity, Settings, Zap, Cpu, Maximize, RotateCw, Wind, BarChart3, Loader2, Sparkles, AlertCircle, Download } from 'lucide-react';
import { DesignSpecs, OptimizationResult } from '../../types';
import { calculateMachineDesign, generateEfficiencyCurve } from '../../services/physicsEngine';
import { suggestSpecsFromDescription, analyzeDesignResult } from '../../services/geminiService';
import LinearMachineView from '../../components/LinearMachineView';
import { EfficiencyChart, DonutChart } from '../../components/Charts';
import { InputGroup, KpiCard, ResultRow } from '../../components/DesignHelpers';
import axios from 'axios';
import CsvVisualizer from '../../components/CsvVisualizer';
import CsvChartVisualizer from '../../components/CsvChartVisualizer';
import PdfViewer from '../../components/PdfViewer';
import WindingDiagrams from '../../components/WindingDiagrams';

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
    const [specPerformanceDict, setSpecPerformanceDict] = useState<Record<string, any> | null>(null);
    const [showWindingInfo, setShowWindingInfo] = useState<boolean>(true);

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
                        parameters: response.data.parameters,
                        // Fields from lines 867-870
                        path2SwarmData: response.data.path2SwarmData,
                        project_name: response.data.project_name,
                        expected_project_file: response.data.expected_project_file,
                        path2FEACsv: response.data.path2FEACsv
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

                    // Extract spec_performance_dict
                    if (response.data.spec_performance_dict && typeof response.data.spec_performance_dict === 'object') {
                        console.log("spec_performance_dict found:", response.data.spec_performance_dict);
                        setSpecPerformanceDict(response.data.spec_performance_dict);
                    } else {
                        console.warn("spec_performance_dict not found or invalid");
                        setSpecPerformanceDict(null);
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

            <main className="flex-1 flex flex-col overflow-hidden h-[calc(100vh-73px)]">
                {/* Top Section: Controls & AI Input - Collapsible */}
                <div className="bg-card border-b border-border">
                    <div className="p-4">
                        <div className="grid grid-cols-1 lg:grid-cols-3 gap-4">
                            {/* AI Design Assistant */}
                            <div className="bg-muted/50 rounded-lg p-4 border border-border">
                                <h2 className="text-sm font-semibold text-primary mb-3 flex items-center">
                                    <Sparkles className="w-4 h-4 mr-2" /> Design Assistant
                                </h2>
                                <div className="space-y-2">
                                    <textarea
                                        value={userPrompt}
                                        onChange={(e) => setUserPrompt(e.target.value)}
                                        placeholder="e.g., I need a high-speed motor for a drone, 500W, 12V..."
                                        className="w-full bg-background border border-border rounded-md p-2 text-sm text-foreground focus:ring-1 focus:ring-ring outline-none resize-none h-16 placeholder-muted-foreground"
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
                            <div className="bg-muted/50 rounded-lg p-4 border border-border">
                                <h2 className="text-sm font-semibold text-foreground mb-3 flex items-center">
                                    <Settings className="w-4 h-4 mr-2" /> Parameters
                                </h2>
                                <div className="space-y-3">
                                    <div className="grid grid-cols-2 gap-2">
                                        <InputGroup label="Power (kW)" value={specs.ratedPower} onChange={v => handleInputChange('ratedPower', v)} icon={<Zap className="w-3 h-3" />} />
                                        <InputGroup label="Speed (RPM)" value={specs.ratedSpeed} onChange={v => handleInputChange('ratedSpeed', v)} icon={<Activity className="w-3 h-3" />} />
                                        <InputGroup label="Voltage (V)" value={specs.ratedVoltage} onChange={v => handleInputChange('ratedVoltage', v)} />
                                        <InputGroup label="Current Dens. (A/mm²)" value={specs.currentDensity} onChange={v => handleInputChange('currentDensity', v)} />
                                        <InputGroup label="Outer Dia. (mm)" value={specs.outerDiameterLimit} onChange={v => handleInputChange('outerDiameterLimit', v)} icon={<Maximize className="w-3 h-3" />} />
                                        <InputGroup label="Length (mm)" value={specs.axialLengthLimit} onChange={v => handleInputChange('axialLengthLimit', v)} />
                                        <InputGroup label="Air Gap (mm)" value={specs.airGap} step={0.1} onChange={v => handleInputChange('airGap', v)} />
                                        <InputGroup label="Slot Count" value={specs.slotCount} step={3} onChange={v => handleInputChange('slotCount', v)} />
                                        <InputGroup label="Pole Count" value={specs.poleCount} step={2} onChange={v => handleInputChange('poleCount', v)} />
                                    </div>
                                </div>
                            </div>

                            {/* CSV Chart Visualizer */}
                            {path2FEACsv && (
                                <div className="bg-muted/50 rounded-lg p-4 border border-border">
                                    <h2 className="text-sm font-semibold text-foreground mb-3 flex items-center">
                                        <BarChart3 className="w-4 h-4 mr-2" /> CSV 图表可视化
                                    </h2>
                                    <div className="h-[300px]">
                                        <CsvChartVisualizer
                                            path2FEACsv={path2FEACsv}
                                            projectName={metadata?.name}
                                        />
                                    </div>
                                </div>
                            )}
                        </div>
                    </div>
                </div>

                {/* Main Content: Visualization & Results - One Column */}
                <div className="flex-1 bg-background flex flex-col h-full overflow-hidden">

                    <div className="flex-1 overflow-y-auto p-6">
                        {/* Visualization Panel - Full Width */}
                        <div className="mb-8">
                            <div className="flex items-center justify-between mb-4">
                                <h3 className="text-lg font-medium text-foreground">Geometry</h3>
                            </div>

                            <div className="grid grid-cols-1 lg:grid-cols-10 gap-6 mb-6">
                                {/* PDF Cross Section Viewer - 30% width (3 columns, 60% of original 50%) */}
                                <div className="lg:col-span-3 h-[500px] bg-card rounded-lg border border-border overflow-hidden">
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

                                {/* Parameters Display - Classified by Type - 70% width (7 columns) */}
                                <div className="lg:col-span-7 h-[500px] bg-card rounded-lg border border-border overflow-y-auto">
                                    {isLoadingMetadata ? (
                                        <div className="h-full flex items-center justify-center">
                                            <div className="text-center">
                                                <Loader2 className="w-6 h-6 animate-spin text-primary mx-auto mb-2" />
                                                <p className="text-sm text-muted-foreground">加载参数中...</p>
                                            </div>
                                        </div>
                                    ) : metadata?.parameters ? (
                                        <div className="p-2">
                                            <h4 className="text-sm font-semibold text-foreground mb-2">Parameters</h4>
                                            {(() => {
                                                const params = metadata.parameters;
                                                const fixedParams: Array<[string, any]> = [];
                                                const freeParams: Array<[string, any]> = [];
                                                const derivedParams: Array<[string, any]> = [];

                                                Object.entries(params).forEach(([key, param]: [string, any]) => {
                                                    if (param?.type === 'fixed') {
                                                        fixedParams.push([key, param]);
                                                    } else if (param?.type === 'free') {
                                                        freeParams.push([key, param]);
                                                    } else if (param?.type === 'derived') {
                                                        derivedParams.push([key, param]);
                                                    }
                                                });

                                                return (
                                                    <div className="space-y-3">
                                                        {/* Fixed Parameters */}
                                                        {fixedParams.length > 0 && (
                                                            <div>
                                                                <h5 className="text-xs font-semibold text-muted-foreground uppercase tracking-wider mb-2 flex items-center">
                                                                    <span className="w-2 h-2 rounded-full bg-blue-500 mr-2"></span>
                                                                    Fixed Parameters ({fixedParams.length})
                                                                </h5>
                                                                <div className="grid grid-cols-7 gap-1.5">
                                                                    {fixedParams.map(([key, param]) => (
                                                                        <div key={key} className="bg-muted/30 rounded p-1 border border-border/50">
                                                                            <div className="flex flex-col">
                                                                                <div className="text-xs font-mono font-semibold text-foreground truncate" title={key}>{key}</div>
                                                                                <div className="text-xs text-muted-foreground truncate" title={param.name || key}>{param.name || key}</div>
                                                                                <div className="text-xs font-mono text-foreground mt-0.5">
                                                                                    {typeof param.value === 'number' ? param.value.toFixed(4) : String(param.value)}
                                                                                </div>
                                                                                {param.unit && (
                                                                                    <div className="text-xs text-muted-foreground">{param.unit}</div>
                                                                                )}
                                                                            </div>
                                                                        </div>
                                                                    ))}
                                                                </div>
                                                            </div>
                                                        )}

                                                        {/* Free Parameters */}
                                                        {freeParams.length > 0 && (
                                                            <div>
                                                                <h5 className="text-xs font-semibold text-muted-foreground uppercase tracking-wider mb-2 flex items-center">
                                                                    <span className="w-2 h-2 rounded-full bg-green-500 mr-2"></span>
                                                                    Free Parameters ({freeParams.length})
                                                                </h5>
                                                                <div className="grid grid-cols-7 gap-1.5">
                                                                    {freeParams.map(([key, param]) => (
                                                                        <div key={key} className="bg-muted/30 rounded p-1 border border-border/50">
                                                                            <div className="flex flex-col">
                                                                                <div className="text-xs font-mono font-semibold text-foreground truncate" title={key}>{key}</div>
                                                                                <div className="text-xs text-muted-foreground truncate" title={param.name || key}>{param.name || key}</div>
                                                                                <div className="text-xs font-mono text-foreground mt-0.5">
                                                                                    {typeof param.value === 'number' ? param.value.toFixed(4) : String(param.value)}
                                                                                </div>
                                                                                {param.unit && (
                                                                                    <div className="text-xs text-muted-foreground">{param.unit}</div>
                                                                                )}
                                                                                {param.bounds && Array.isArray(param.bounds) && (
                                                                                    <div className="text-xs text-muted-foreground mt-0.5 truncate" title={`Bounds: [${param.bounds[0]?.toFixed(4)}, ${param.bounds[1]?.toFixed(4)}]`}>
                                                                                        [{param.bounds[0]?.toFixed(2)}, {param.bounds[1]?.toFixed(2)}]
                                                                                    </div>
                                                                                )}
                                                                            </div>
                                                                        </div>
                                                                    ))}
                                                                </div>
                                                            </div>
                                                        )}

                                                        {/* Derived Parameters */}
                                                        {derivedParams.length > 0 && (
                                                            <div>
                                                                <h5 className="text-xs font-semibold text-muted-foreground uppercase tracking-wider mb-2 flex items-center">
                                                                    <span className="w-2 h-2 rounded-full bg-purple-500 mr-2"></span>
                                                                    Derived Parameters ({derivedParams.length})
                                                                </h5>
                                                                <div className="grid grid-cols-7 gap-1.5">
                                                                    {derivedParams.map(([key, param]) => (
                                                                        <div key={key} className="bg-muted/30 rounded p-1 border border-border/50">
                                                                            <div className="flex flex-col">
                                                                                <div className="text-xs font-mono font-semibold text-foreground truncate" title={key}>{key}</div>
                                                                                <div className="text-xs text-muted-foreground truncate" title={param.name || key}>{param.name || key}</div>
                                                                                <div className="text-xs font-mono text-foreground mt-0.5">
                                                                                    {typeof param.value === 'number' ? param.value.toFixed(4) : String(param.value)}
                                                                                </div>
                                                                                {param.unit && (
                                                                                    <div className="text-xs text-muted-foreground">{param.unit}</div>
                                                                                )}
                                                                                {param._calc_info?.source && (
                                                                                    <div className="text-xs text-muted-foreground mt-0.5 italic truncate" title={param._calc_info.source}>
                                                                                        calc: {param._calc_info.source.substring(0, 30)}...
                                                                                    </div>
                                                                                )}
                                                                            </div>
                                                                        </div>
                                                                    ))}
                                                                </div>
                                                            </div>
                                                        )}
                                                    </div>
                                                );
                                            })()}
                                        </div>
                                    ) : (
                                        <div className="h-full flex items-center justify-center">
                                            <div className="text-center p-6">
                                                <AlertCircle className="w-8 h-8 text-muted-foreground mx-auto mb-2" />
                                                <p className="text-sm text-muted-foreground">参数未加载</p>
                                            </div>
                                        </div>
                                    )}
                                </div>
                            </div>

                            {/* Linear Geometry View - Full Width with More Space */}
                            {wilyData && (
                                <div className="w-full mb-6">
                                    <div className="flex items-center justify-between mb-4">
                                        <h4 className="text-sm font-medium text-foreground">
                                            Linear View (Auto-scaled)
                                            {coilPitchY !== null && (
                                                <span className="ml-2 text-muted-foreground font-normal">
                                                    - coil_pitch_y = {coilPitchY}
                                                </span>
                                            )}
                                        </h4>
                                        <label className="flex items-center space-x-2 cursor-pointer">
                                            <input
                                                type="checkbox"
                                                checked={showWindingInfo}
                                                onChange={(e) => setShowWindingInfo(e.target.checked)}
                                                className="w-4 h-4 rounded border-border text-primary focus:ring-primary"
                                            />
                                            <span className="text-xs text-muted-foreground">显示绕组信息</span>
                                        </label>
                                    </div>
                                    <div className="w-full bg-card rounded-lg border border-border p-4" style={{ minHeight: '700px' }}>
                                        <LinearMachineView
                                            Qs={wilyData.stator_slot_number_Qs || wilyData.Qs || 12}
                                            p={wilyData.pole_pair_number_p || wilyData.p || 2}
                                            ps={wilyData.suspension_pole_pair_number_ps || wilyData.ps || 3}
                                            coilPitchY={coilPitchY}
                                            layer_X_phases={wilyData.layer_X_phases || null}
                                            layer_Y_phases={wilyData.layer_Y_phases || null}
                                            layer_X_signs={wilyData.layer_X_signs || null}
                                            layer_Y_signs={wilyData.layer_Y_signs || null}
                                            grouping_AC={wilyData.grouping_AC || null}
                                        />
                                        {showWindingInfo && (
                                            <div className="mt-6 pt-6 border-t border-border">
                                                <h5 className="text-xs font-semibold text-foreground mb-3 uppercase tracking-wider">绕组信息 (Winding Information)</h5>
                                                <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
                                                    {/* Winding Factor */}
                                                    {wilyData.kw1 !== undefined && (
                                                        <div className="bg-muted/30 rounded p-3 border border-border/50">
                                                            <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold mb-1">绕组因子 (kw1)</div>
                                                            <div className="text-sm font-mono text-foreground">{typeof wilyData.kw1 === 'number' ? wilyData.kw1.toFixed(4) : String(wilyData.kw1)}</div>
                                                        </div>
                                                    )}

                                                    {/* Layer X Phases */}
                                                    {wilyData.layer_X_phases && Array.isArray(wilyData.layer_X_phases) && (
                                                        <div className="bg-muted/30 rounded p-3 border border-border/50">
                                                            <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold mb-2">Layer X 相序</div>
                                                            <div className="text-xs font-mono text-foreground break-all">
                                                                {wilyData.layer_X_phases.map((phase: any, idx: number) => {
                                                                    let phaseStr = '-';
                                                                    if (phase === null || phase === undefined) {
                                                                        phaseStr = '-';
                                                                    } else if (typeof phase === 'string') {
                                                                        phaseStr = phase;
                                                                    } else if (typeof phase === 'number') {
                                                                        phaseStr = String(phase);
                                                                    } else if (typeof phase === 'object') {
                                                                        // Try to extract meaningful value from object
                                                                        phaseStr = phase.toString ? phase.toString() : (phase.value || phase.phase || JSON.stringify(phase));
                                                                    } else {
                                                                        phaseStr = String(phase);
                                                                    }
                                                                    return <span key={idx} className="inline-block mr-1">{phaseStr}</span>;
                                                                })}
                                                            </div>
                                                        </div>
                                                    )}

                                                    {/* Layer X Signs */}
                                                    {wilyData.layer_X_signs && Array.isArray(wilyData.layer_X_signs) && (
                                                        <div className="bg-muted/30 rounded p-3 border border-border/50">
                                                            <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold mb-2">Layer X 符号</div>
                                                            <div className="text-xs font-mono text-foreground break-all">
                                                                {wilyData.layer_X_signs.map((sign: any, idx: number) => {
                                                                    let signStr = '-';
                                                                    if (sign === null || sign === undefined) {
                                                                        signStr = '-';
                                                                    } else if (typeof sign === 'string') {
                                                                        signStr = sign;
                                                                    } else if (typeof sign === 'number') {
                                                                        signStr = String(sign);
                                                                    } else if (typeof sign === 'object') {
                                                                        signStr = sign.toString ? sign.toString() : (sign.value || sign.sign || JSON.stringify(sign));
                                                                    } else {
                                                                        signStr = String(sign);
                                                                    }
                                                                    return <span key={idx} className="inline-block mr-1">{signStr}</span>;
                                                                })}
                                                            </div>
                                                        </div>
                                                    )}

                                                    {/* Layer Y Phases */}
                                                    {wilyData.layer_Y_phases && Array.isArray(wilyData.layer_Y_phases) && (
                                                        <div className="bg-muted/30 rounded p-3 border border-border/50">
                                                            <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold mb-2">Layer Y 相序</div>
                                                            <div className="text-xs font-mono text-foreground break-all">
                                                                {wilyData.layer_Y_phases.map((phase: any, idx: number) => {
                                                                    let phaseStr = '-';
                                                                    if (phase === null || phase === undefined) {
                                                                        phaseStr = '-';
                                                                    } else if (typeof phase === 'string') {
                                                                        phaseStr = phase;
                                                                    } else if (typeof phase === 'number') {
                                                                        phaseStr = String(phase);
                                                                    } else if (typeof phase === 'object') {
                                                                        phaseStr = phase.toString ? phase.toString() : (phase.value || phase.phase || JSON.stringify(phase));
                                                                    } else {
                                                                        phaseStr = String(phase);
                                                                    }
                                                                    return <span key={idx} className="inline-block mr-1">{phaseStr}</span>;
                                                                })}
                                                            </div>
                                                        </div>
                                                    )}

                                                    {/* Layer Y Signs */}
                                                    {wilyData.layer_Y_signs && Array.isArray(wilyData.layer_Y_signs) && (
                                                        <div className="bg-muted/30 rounded p-3 border border-border/50">
                                                            <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold mb-2">Layer Y 符号</div>
                                                            <div className="text-xs font-mono text-foreground break-all">
                                                                {wilyData.layer_Y_signs.map((sign: any, idx: number) => {
                                                                    let signStr = '-';
                                                                    if (sign === null || sign === undefined) {
                                                                        signStr = '-';
                                                                    } else if (typeof sign === 'string') {
                                                                        signStr = sign;
                                                                    } else if (typeof sign === 'number') {
                                                                        signStr = String(sign);
                                                                    } else if (typeof sign === 'object') {
                                                                        signStr = sign.toString ? sign.toString() : (sign.value || sign.sign || JSON.stringify(sign));
                                                                    } else {
                                                                        signStr = String(sign);
                                                                    }
                                                                    return <span key={idx} className="inline-block mr-1">{signStr}</span>;
                                                                })}
                                                            </div>
                                                        </div>
                                                    )}

                                                    {/* Grouping AC */}
                                                    {wilyData.grouping_AC && Array.isArray(wilyData.grouping_AC) && (
                                                        <div className="bg-muted/30 rounded p-3 border border-border/50">
                                                            <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold mb-2">Grouping AC</div>
                                                            <div className="text-xs font-mono text-foreground break-all">
                                                                {wilyData.grouping_AC.map((group: any, idx: number) => {
                                                                    const groupValue = group === 1 || group === true || group === '1' ? '1' : '0';
                                                                    return <span key={idx} className="inline-block mr-1">{groupValue}</span>;
                                                                })}
                                                            </div>
                                                        </div>
                                                    )}

                                                    {/* Additional Winding Information */}
                                                    {wilyData.SIict_kw_els && typeof wilyData.SIict_kw_els === 'object' && (
                                                        <div className="bg-muted/30 rounded p-3 border border-border/50 md:col-span-2 lg:col-span-3">
                                                            <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold mb-2">绕组因子详情 (SIict_kw_els)</div>
                                                            <div className="text-xs font-mono text-foreground space-y-1">
                                                                {Object.entries(wilyData.SIict_kw_els).map(([key, value]: [string, any]) => (
                                                                    <div key={key} className="flex justify-between">
                                                                        <span className="text-muted-foreground">{key}:</span>
                                                                        <span>{typeof value === 'number' ? value.toFixed(4) : String(value)}</span>
                                                                    </div>
                                                                ))}
                                                            </div>
                                                        </div>
                                                    )}

                                                    {wilyData.SIict_kw_cjh && typeof wilyData.SIict_kw_cjh === 'object' && (
                                                        <div className="bg-muted/30 rounded p-3 border border-border/50 md:col-span-2 lg:col-span-3">
                                                            <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold mb-2">绕组因子详情 (SIict_kw_cjh)</div>
                                                            <div className="text-xs font-mono text-foreground space-y-1">
                                                                {Object.entries(wilyData.SIict_kw_cjh).map(([key, value]: [string, any]) => (
                                                                    <div key={key} className="flex justify-between">
                                                                        <span className="text-muted-foreground">{key}:</span>
                                                                        <span>{typeof value === 'number' ? value.toFixed(4) : String(value)}</span>
                                                                    </div>
                                                                ))}
                                                            </div>
                                                        </div>
                                                    )}
                                                </div>

                                                {/* Winding Diagrams */}
                                                <div className="mt-6 pt-6 border-t border-border">
                                                    <h5 className="text-xs font-semibold text-foreground mb-3 uppercase tracking-wider">绕组相量图 (Winding Phasor Diagrams)</h5>
                                                    <WindingDiagrams
                                                        Qs={wilyData.stator_slot_number_Qs || wilyData.Qs || 12}
                                                        p={wilyData.pole_pair_number_p || wilyData.p || 2}
                                                        m={wilyData.m || 3}
                                                        layer_X_phases={wilyData.layer_X_phases || []}
                                                        layer_X_signs={wilyData.layer_X_signs || []}
                                                    />
                                                </div>
                                            </div>
                                        )}
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
                                    <div className="grid grid-cols-2 md:grid-cols-4 lg:grid-cols-5 gap-2 p-3">
                                        <div className="space-y-0.5">
                                            <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">名称</div>
                                            <div className="text-xs font-mono text-foreground truncate" title={metadata.name || 'N/A'}>{metadata.name || 'N/A'}</div>
                                        </div>
                                        <div className="space-y-0.5">
                                            <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">永磁体</div>
                                            <div className="text-xs text-foreground">
                                                {metadata.bool_PermanentMagnet ? (
                                                    <span className="text-emerald-500">是</span>
                                                ) : (
                                                    <span className="text-red-500">否</span>
                                                )}
                                            </div>
                                        </div>
                                        <div className="space-y-0.5">
                                            <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">定子槽封闭</div>
                                            <div className="text-xs text-foreground">
                                                {metadata.bool_StatorSlotClosed ? (
                                                    <span className="text-emerald-500">是</span>
                                                ) : (
                                                    <span className="text-red-500">否</span>
                                                )}
                                            </div>
                                        </div>
                                        <div className="space-y-0.5">
                                            <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">转子开槽</div>
                                            <div className="text-xs text-foreground">
                                                {metadata.bool_RotorNotched ? (
                                                    <span className="text-emerald-500">是</span>
                                                ) : (
                                                    <span className="text-red-500">否</span>
                                                )}
                                            </div>
                                        </div>
                                        <div className="space-y-0.5">
                                            <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">FEA 工具</div>
                                            <div className="text-xs font-mono text-foreground truncate" title={metadata.select_FEA_tool || 'N/A'}>{metadata.select_FEA_tool || 'N/A'}</div>
                                        </div>
                                        <div className="space-y-0.5">
                                            <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">FEA 配置</div>
                                            <div className="text-xs text-foreground truncate" title={metadata.select_fea_config_dict || 'N/A'}>{metadata.select_fea_config_dict || 'N/A'}</div>
                                        </div>
                                        <div className="space-y-0.5">
                                            <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">计算后删除结果</div>
                                            <div className="text-xs text-foreground">
                                                {metadata.bool_jmagDeleteResultsAfterCalculation ? (
                                                    <span className="text-emerald-500">是</span>
                                                ) : (
                                                    <span className="text-red-500">否</span>
                                                )}
                                            </div>
                                        </div>
                                        <div className="space-y-0.5">
                                            <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">计数器</div>
                                            <div className="text-xs font-mono text-foreground">{metadata.counter ?? 'N/A'}</div>
                                        </div>
                                        <div className="space-y-0.5">
                                            <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">参数数量</div>
                                            <div className="text-xs font-mono text-foreground">
                                                {metadata.parameters && typeof metadata.parameters === 'object'
                                                    ? Object.keys(metadata.parameters).length
                                                    : 0}
                                            </div>
                                        </div>
                                    </div>

                                    {/* Path and Project Information (lines 867-870) */}
                                    <div className="border-t border-border pt-3 mt-3">
                                        <h4 className="text-xs font-semibold text-foreground mb-2 flex items-center">
                                            <Settings className="w-3 h-3 mr-1" /> 路径和项目信息
                                        </h4>
                                        <div className="grid grid-cols-1 md:grid-cols-2 gap-2">
                                            {metadata.path2SwarmData && (
                                                <div className="space-y-1">
                                                    <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">Swarm 数据路径</div>
                                                    <div className="text-sm font-mono text-foreground break-words">{metadata.path2SwarmData}</div>
                                                </div>
                                            )}
                                            {metadata.project_name && (
                                                <div className="space-y-1">
                                                    <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">项目名称</div>
                                                    <div className="text-sm font-mono text-foreground">{metadata.project_name}</div>
                                                </div>
                                            )}
                                            {metadata.expected_project_file && (
                                                <div className="space-y-1 md:col-span-2">
                                                    <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">预期项目文件</div>
                                                    <div className="text-sm font-mono text-foreground break-words">{metadata.expected_project_file}</div>
                                                </div>
                                            )}
                                            {metadata.path2FEACsv && (
                                                <div className="space-y-1 md:col-span-2">
                                                    <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">FEA CSV 路径</div>
                                                    <div className="text-sm font-mono text-foreground break-words">{metadata.path2FEACsv}</div>
                                                </div>
                                            )}
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
                                        <div className="bg-muted/50 px-3 py-2 border-b border-border">
                                            <h4 className="text-xs font-semibold text-foreground flex items-center">
                                                <Maximize className="w-3 h-3 mr-1" /> 几何参数
                                            </h4>
                                        </div>
                                        <div className="grid grid-cols-2 md:grid-cols-4 lg:grid-cols-5 gap-2 p-3">
                                            <div className="space-y-0.5">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">叠片长度</div>
                                                <div className="text-xs font-mono text-foreground">{exData.mm_stack_length_specified?.toFixed(2) || 'N/A'} <span className="text-muted-foreground">mm</span></div>
                                            </div>
                                            <div className="space-y-0.5">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">槽面积</div>
                                                <div className="text-xs font-mono text-foreground">{exData.mm2_slot_area?.toFixed(2) || 'N/A'} <span className="text-muted-foreground">mm²</span></div>
                                            </div>
                                            <div className="space-y-0.5">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">磁体面积</div>
                                                <div className="text-xs font-mono text-foreground">{exData.mm2_magnet_area?.toFixed(2) || 'N/A'} <span className="text-muted-foreground">mm²</span></div>
                                            </div>
                                            <div className="space-y-0.5">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">磁体起始角度</div>
                                                <div className="text-xs font-mono text-foreground">{exData.Magnet_StartAngle?.toFixed(2) || 'N/A'} <span className="text-muted-foreground">°</span></div>
                                            </div>
                                            <div className="space-y-0.5">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">初始旋转角度</div>
                                                <div className="text-xs font-mono text-foreground">{exData.InitialRotationAngle?.toFixed(2) || 'N/A'} <span className="text-muted-foreground">°</span></div>
                                            </div>
                                        </div>
                                    </div>

                                    {/* 材料参数 */}
                                    <div className="bg-card rounded-lg border border-border overflow-hidden">
                                        <div className="bg-muted/50 px-3 py-2 border-b border-border">
                                            <h4 className="text-xs font-semibold text-foreground flex items-center">
                                                <Cpu className="w-3 h-3 mr-1" /> 材料参数
                                            </h4>
                                        </div>
                                        <div className="grid grid-cols-2 md:grid-cols-4 lg:grid-cols-5 gap-2 p-3">
                                            <div className="space-y-0.5">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">磁体材料</div>
                                                <div className="text-xs font-mono text-foreground truncate" title={exData.Magnet_Name || 'N/A'}>{exData.Magnet_Name || 'N/A'}</div>
                                            </div>
                                            <div className="space-y-0.5">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">钢材材料</div>
                                                <div className="text-xs font-mono text-foreground truncate" title={exData.SteelMaterial || 'N/A'}>{exData.SteelMaterial || 'N/A'}</div>
                                            </div>
                                            <div className="space-y-0.5">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">定子铁心材料</div>
                                                <div className="text-xs font-mono text-foreground truncate" title={exData.StatorCore_Material || 'N/A'}>{exData.StatorCore_Material || 'N/A'}</div>
                                            </div>
                                            <div className="space-y-0.5">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">转子铁心材料</div>
                                                <div className="text-xs font-mono text-foreground truncate" title={exData.RotorCore_Material || 'N/A'}>{exData.RotorCore_Material || 'N/A'}</div>
                                            </div>
                                            <div className="space-y-0.5">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">叠片系数</div>
                                                <div className="text-xs font-mono text-foreground">{exData.LaminationFactor?.toFixed(1) || 'N/A'} <span className="text-muted-foreground">%</span></div>
                                            </div>
                                        </div>
                                    </div>

                                    {/* 电气参数 */}
                                    <div className="bg-card rounded-lg border border-border overflow-hidden">
                                        <div className="bg-muted/50 px-3 py-2 border-b border-border">
                                            <h4 className="text-xs font-semibold text-foreground flex items-center">
                                                <Zap className="w-3 h-3 mr-1" /> 电气参数
                                            </h4>
                                        </div>
                                        <div className="grid grid-cols-2 md:grid-cols-4 lg:grid-cols-5 gap-2 p-3">
                                            <div className="space-y-0.5">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">额定功率</div>
                                                <div className="text-xs font-mono text-foreground">{(exData.RatedPower / 1000)?.toFixed(1) || 'N/A'} <span className="text-muted-foreground">kW</span></div>
                                            </div>
                                            <div className="space-y-0.5">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">额定转速</div>
                                                <div className="text-xs font-mono text-foreground">{exData.RatedSpeed?.toLocaleString() || 'N/A'} <span className="text-muted-foreground">RPM</span></div>
                                            </div>
                                            <div className="space-y-0.5">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">激励频率</div>
                                                <div className="text-xs font-mono text-foreground">{exData.ExcitationFreqSimulated?.toFixed(1) || 'N/A'} <span className="text-muted-foreground">Hz</span></div>
                                            </div>
                                            <div className="space-y-0.5">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">直流母线电压</div>
                                                <div className="text-xs font-mono text-foreground">{exData.DCBusVoltage?.toFixed(0) || 'N/A'} <span className="text-muted-foreground">V</span></div>
                                            </div>
                                            <div className="space-y-0.5">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">电流密度</div>
                                                <div className="text-xs font-mono text-foreground">{(exData.Js / 1e6)?.toFixed(1) || 'N/A'} <span className="text-muted-foreground">A/mm²</span></div>
                                            </div>
                                            <div className="space-y-0.5">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">驱动绕组电阻</div>
                                                <div className="text-xs font-mono text-foreground">{exData.DriveW_Rs?.toFixed(2) || 'N/A'} <span className="text-muted-foreground">Ω</span></div>
                                            </div>
                                            <div className="space-y-0.5">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">悬浮绕组电阻</div>
                                                <div className="text-xs font-mono text-foreground">{exData.BeariW_Rs?.toFixed(2) || 'N/A'} <span className="text-muted-foreground">Ω</span></div>
                                            </div>
                                            <div className="space-y-0.5">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">连接方式</div>
                                                <div className="text-xs text-foreground">
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
                                        <div className="bg-muted/50 px-3 py-2 border-b border-border">
                                            <h4 className="text-xs font-semibold text-foreground flex items-center">
                                                <Activity className="w-3 h-3 mr-1" /> 温度参数
                                            </h4>
                                        </div>
                                        <div className="grid grid-cols-2 md:grid-cols-4 lg:grid-cols-5 gap-2 p-3">
                                            <div className="space-y-0.5">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">工作温度</div>
                                                <div className="text-xs font-mono text-foreground">{exData.Temperature?.toFixed(0) || 'N/A'} <span className="text-muted-foreground">°C</span></div>
                                            </div>
                                            <div className="space-y-0.5">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">磁体温度</div>
                                                <div className="text-xs font-mono text-foreground">{exData.Magnet_Temperature?.toFixed(0) || 'N/A'} <span className="text-muted-foreground">°C</span></div>
                                            </div>
                                        </div>
                                    </div>

                                    {/* 绕组参数 */}
                                    <div className="bg-card rounded-lg border border-border overflow-hidden">
                                        <div className="bg-muted/50 px-3 py-2 border-b border-border">
                                            <h4 className="text-xs font-semibold text-foreground flex items-center">
                                                <Wind className="w-3 h-3 mr-1" /> 绕组参数
                                            </h4>
                                        </div>
                                        <div className="grid grid-cols-2 md:grid-cols-4 lg:grid-cols-5 gap-2 p-3">
                                            <div className="space-y-0.5">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">串联线圈匝数</div>
                                                <div className="text-xs font-mono text-foreground">{exData.no_series_coil_turns_N?.toFixed(0) || 'N/A'}</div>
                                            </div>
                                            <div className="space-y-0.5">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">驱动绕组每槽导体数</div>
                                                <div className="text-xs font-mono text-foreground">{exData.DriveW_zQ?.toFixed(0) || 'N/A'}</div>
                                            </div>
                                            <div className="space-y-0.5">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">悬浮绕组每槽导体数</div>
                                                <div className="text-xs font-mono text-foreground">{exData.BeariW_zQ?.toFixed(0) || 'N/A'}</div>
                                            </div>
                                            <div className="space-y-0.5">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">绕组填充因子</div>
                                                <div className="text-xs font-mono text-foreground">{(exData.WindingFill * 100)?.toFixed(2) || 'N/A'} <span className="text-muted-foreground">%</span></div>
                                            </div>
                                            <div className="space-y-0.5">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">转矩电流比例</div>
                                                <div className="text-xs font-mono text-foreground">{(exData.TORQUE_CURRENT_RATIO * 100)?.toFixed(1) || 'N/A'} <span className="text-muted-foreground">%</span></div>
                                            </div>
                                            <div className="space-y-0.5">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">悬浮电流比例</div>
                                                <div className="text-xs font-mono text-foreground">{(exData.SUSPENSION_CURRENT_RATIO * 100)?.toFixed(1) || 'N/A'} <span className="text-muted-foreground">%</span></div>
                                            </div>
                                            <div className="space-y-0.5">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">槽电流利用率（转矩）</div>
                                                <div className="text-xs font-mono text-foreground">{(exData.slot_current_utilizing_ratio_for_torque * 100)?.toFixed(2) || 'N/A'} <span className="text-muted-foreground">%</span></div>
                                            </div>
                                        </div>
                                    </div>

                                    {/* 电流参数 */}
                                    <div className="bg-card rounded-lg border border-border overflow-hidden">
                                        <div className="bg-muted/50 px-3 py-2 border-b border-border">
                                            <h4 className="text-xs font-semibold text-foreground flex items-center">
                                                <BarChart3 className="w-3 h-3 mr-1" /> 电流参数
                                            </h4>
                                        </div>
                                        <div className="grid grid-cols-2 md:grid-cols-4 lg:grid-cols-5 gap-2 p-3">
                                            <div className="space-y-0.5">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">槽内电流</div>
                                                <div className="text-xs font-mono text-foreground">{exData.CurrentAmp_in_the_slot?.toFixed(2) || 'N/A'} <span className="text-muted-foreground">A</span></div>
                                            </div>
                                            <div className="space-y-0.5">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">每导体电流</div>
                                                <div className="text-xs font-mono text-foreground">{exData.CurrentAmp_per_conductor?.toFixed(2) || 'N/A'} <span className="text-muted-foreground">A</span></div>
                                            </div>
                                            <div className="space-y-0.5">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">每相电流</div>
                                                <div className="text-xs font-mono text-foreground">{exData.CurrentAmp_per_phase?.toFixed(2) || 'N/A'} <span className="text-muted-foreground">A</span></div>
                                            </div>
                                            <div className="space-y-0.5">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">驱动绕组电流</div>
                                                <div className="text-xs font-mono text-foreground">{exData.DriveW_CurrentAmp?.toFixed(2) || 'N/A'} <span className="text-muted-foreground">A</span></div>
                                            </div>
                                            <div className="space-y-0.5">
                                                <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">悬浮绕组电流</div>
                                                <div className="text-xs font-mono text-foreground">{exData.BeariW_CurrentAmp?.toFixed(2) || 'N/A'} <span className="text-muted-foreground">A</span></div>
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

                        {/* Spec Performance Dictionary Section */}
                        {specPerformanceDict && (
                            <div className="mb-8">
                                <h3 className="text-lg font-medium text-foreground mb-4 flex items-center">
                                    <BarChart3 className="w-5 h-5 mr-2" /> Performance Metrics (spec_performance_dict)
                                </h3>
                                <div className="space-y-6">
                                    {/* Optimization Parameters */}
                                    {specPerformanceDict.x_denorm_dict && (
                                        <div className="bg-card rounded-lg border border-border overflow-hidden">
                                            <div className="bg-muted/50 px-3 py-2 border-b border-border">
                                                <h4 className="text-xs font-semibold text-foreground flex items-center">
                                                    <Settings className="w-3 h-3 mr-1" /> Optimization Parameters (x_denorm_dict)
                                                </h4>
                                            </div>
                                            <div className="grid grid-cols-2 md:grid-cols-4 lg:grid-cols-5 gap-2 p-3">
                                                {Object.entries(specPerformanceDict.x_denorm_dict).map(([key, value]) => (
                                                    <div key={key} className="space-y-0.5">
                                                        <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold truncate" title={key.replace(/_/g, ' ')}>
                                                            {key.replace(/_/g, ' ')}
                                                        </div>
                                                        <div className="text-xs font-mono text-foreground">
                                                            {typeof value === 'number' ? value.toFixed(4) : String(value)}
                                                        </div>
                                                    </div>
                                                ))}
                                            </div>
                                        </div>
                                    )}

                                    {/* Fitness Values */}
                                    <div className="bg-card rounded-lg border border-border overflow-hidden">
                                        <div className="bg-muted/50 px-3 py-2 border-b border-border">
                                            <h4 className="text-xs font-semibold text-foreground flex items-center">
                                                <Activity className="w-3 h-3 mr-1" /> Fitness Values
                                            </h4>
                                        </div>
                                        <div className="grid grid-cols-2 md:grid-cols-4 lg:grid-cols-5 gap-2 p-3">
                                            {specPerformanceDict.f1 !== undefined && (
                                                <div className="space-y-0.5">
                                                    <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">F1 (Cost)</div>
                                                    <div className="text-xs font-mono text-foreground">{specPerformanceDict.f1.toFixed(4)}</div>
                                                </div>
                                            )}
                                            {specPerformanceDict.f2 !== undefined && (
                                                <div className="space-y-0.5">
                                                    <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">F2 (Efficiency)</div>
                                                    <div className="text-xs font-mono text-foreground">{specPerformanceDict.f2.toFixed(4)}</div>
                                                </div>
                                            )}
                                            {specPerformanceDict.f3 !== undefined && (
                                                <div className="space-y-0.5">
                                                    <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">F3 (Ripple)</div>
                                                    <div className="text-xs font-mono text-foreground">{specPerformanceDict.f3.toFixed(4)}</div>
                                                </div>
                                            )}
                                            {specPerformanceDict.Cost !== undefined && (
                                                <div className="space-y-0.5">
                                                    <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">Total Cost</div>
                                                    <div className="text-xs font-mono text-foreground">{specPerformanceDict.Cost.toFixed(2)}</div>
                                                </div>
                                            )}
                                        </div>
                                    </div>

                                    {/* Torque and Force Metrics */}
                                    <div className="bg-card rounded-lg border border-border overflow-hidden">
                                        <div className="bg-muted/50 px-3 py-2 border-b border-border">
                                            <h4 className="text-xs font-semibold text-foreground flex items-center">
                                                <RotateCw className="w-3 h-3 mr-1" /> Torque & Force Metrics
                                            </h4>
                                        </div>
                                        <div className="grid grid-cols-2 md:grid-cols-4 lg:grid-cols-5 gap-2 p-3">
                                            {specPerformanceDict.torque_average !== undefined && (
                                                <div className="space-y-0.5">
                                                    <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">Average Torque</div>
                                                    <div className="text-xs font-mono text-foreground">{specPerformanceDict.torque_average.toFixed(4)} <span className="text-muted-foreground">Nm</span></div>
                                                </div>
                                            )}
                                            {specPerformanceDict.TRV !== undefined && (
                                                <div className="space-y-0.5">
                                                    <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">TRV</div>
                                                    <div className="text-xs font-mono text-foreground">{specPerformanceDict.TRV.toFixed(2)}</div>
                                                </div>
                                            )}
                                            {specPerformanceDict.ss_avg_force_magnitude !== undefined && (
                                                <div className="space-y-0.5">
                                                    <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">Avg Force Magnitude</div>
                                                    <div className="text-xs font-mono text-foreground">{specPerformanceDict.ss_avg_force_magnitude.toFixed(4)} <span className="text-muted-foreground">N</span></div>
                                                </div>
                                            )}
                                            {specPerformanceDict.FRW !== undefined && (
                                                <div className="space-y-0.5">
                                                    <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">FRW</div>
                                                    <div className="text-xs font-mono text-foreground">{specPerformanceDict.FRW.toFixed(4)}</div>
                                                </div>
                                            )}
                                            {specPerformanceDict.normalized_torque_ripple !== undefined && (
                                                <div className="space-y-0.5">
                                                    <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">Normalized Torque Ripple</div>
                                                    <div className="text-xs font-mono text-foreground">{(specPerformanceDict.normalized_torque_ripple * 100).toFixed(2)} <span className="text-muted-foreground">%</span></div>
                                                </div>
                                            )}
                                            {specPerformanceDict.normalized_force_error_magnitude !== undefined && (
                                                <div className="space-y-0.5">
                                                    <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">Normalized Force Error</div>
                                                    <div className="text-xs font-mono text-foreground">{(specPerformanceDict.normalized_force_error_magnitude * 100).toFixed(2)} <span className="text-muted-foreground">%</span></div>
                                                </div>
                                            )}
                                            {specPerformanceDict.force_error_angle !== undefined && (
                                                <div className="space-y-0.5">
                                                    <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">Force Error Angle</div>
                                                    <div className="text-xs font-mono text-foreground">{specPerformanceDict.force_error_angle.toFixed(4)} <span className="text-muted-foreground">°</span></div>
                                                </div>
                                            )}
                                        </div>
                                    </div>

                                    {/* Losses */}
                                    <div className="bg-card rounded-lg border border-border overflow-hidden">
                                        <div className="bg-muted/50 px-3 py-2 border-b border-border">
                                            <h4 className="text-xs font-semibold text-foreground flex items-center">
                                                <Zap className="w-3 h-3 mr-1" /> Losses
                                            </h4>
                                        </div>
                                        <div className="grid grid-cols-1 lg:grid-cols-2 gap-3 p-3">
                                            {/* Losses Table */}
                                            <div className="grid grid-cols-2 md:grid-cols-3 lg:grid-cols-3 gap-2">
                                                {specPerformanceDict.rated_total_loss !== undefined && (
                                                    <div className="space-y-0.5 bg-muted/30 rounded p-2 border border-border/50">
                                                        <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">Total Loss</div>
                                                        <div className="text-xs font-mono text-foreground">{specPerformanceDict.rated_total_loss.toFixed(2)} <span className="text-muted-foreground">W</span></div>
                                                    </div>
                                                )}
                                                {specPerformanceDict.rated_stator_copper_loss_along_stack !== undefined && (
                                                    <div className="space-y-0.5 bg-muted/30 rounded p-2 border border-border/50">
                                                        <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold truncate" title="Stator Copper Loss (Stack)">Stator Copper Loss (Stack)</div>
                                                        <div className="text-xs font-mono text-foreground">{specPerformanceDict.rated_stator_copper_loss_along_stack.toFixed(2)} <span className="text-muted-foreground">W</span></div>
                                                    </div>
                                                )}
                                                {specPerformanceDict.stator_copper_loss_in_end_turn !== undefined && (
                                                    <div className="space-y-0.5 bg-muted/30 rounded p-2 border border-border/50">
                                                        <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold truncate" title="Stator Copper Loss (End Turn)">Stator Copper Loss (End Turn)</div>
                                                        <div className="text-xs font-mono text-foreground">{specPerformanceDict.stator_copper_loss_in_end_turn.toFixed(2)} <span className="text-muted-foreground">W</span></div>
                                                    </div>
                                                )}
                                                {specPerformanceDict.rated_iron_loss !== undefined && (
                                                    <div className="space-y-0.5 bg-muted/30 rounded p-2 border border-border/50">
                                                        <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">Iron Loss</div>
                                                        <div className="text-xs font-mono text-foreground">{specPerformanceDict.rated_iron_loss.toFixed(2)} <span className="text-muted-foreground">W</span></div>
                                                    </div>
                                                )}
                                                {specPerformanceDict.rated_magnet_Joule_loss !== undefined && (
                                                    <div className="space-y-0.5 bg-muted/30 rounded p-2 border border-border/50">
                                                        <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">Magnet Joule Loss</div>
                                                        <div className="text-xs font-mono text-foreground">{specPerformanceDict.rated_magnet_Joule_loss.toFixed(2)} <span className="text-muted-foreground">W</span></div>
                                                    </div>
                                                )}
                                                {specPerformanceDict.rated_windage_loss !== undefined && (
                                                    <div className="space-y-0.5 bg-muted/30 rounded p-2 border border-border/50">
                                                        <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">Windage Loss</div>
                                                        <div className="text-xs font-mono text-foreground">{specPerformanceDict.rated_windage_loss.toFixed(2)} <span className="text-muted-foreground">W</span></div>
                                                    </div>
                                                )}
                                            </div>

                                            {/* Donut Chart */}
                                            <div className="h-[300px]">
                                                {(() => {
                                                    const lossData = [];
                                                    if (specPerformanceDict.rated_stator_copper_loss_along_stack !== undefined && specPerformanceDict.rated_stator_copper_loss_along_stack > 0) {
                                                        lossData.push({
                                                            name: 'Stator Copper (Stack)',
                                                            value: specPerformanceDict.rated_stator_copper_loss_along_stack
                                                        });
                                                    }
                                                    if (specPerformanceDict.stator_copper_loss_in_end_turn !== undefined && specPerformanceDict.stator_copper_loss_in_end_turn > 0) {
                                                        lossData.push({
                                                            name: 'Stator Copper (End Turn)',
                                                            value: specPerformanceDict.stator_copper_loss_in_end_turn
                                                        });
                                                    }
                                                    if (specPerformanceDict.rated_iron_loss !== undefined && specPerformanceDict.rated_iron_loss > 0) {
                                                        lossData.push({
                                                            name: 'Iron Loss',
                                                            value: specPerformanceDict.rated_iron_loss
                                                        });
                                                    }
                                                    if (specPerformanceDict.rated_magnet_Joule_loss !== undefined && specPerformanceDict.rated_magnet_Joule_loss > 0) {
                                                        lossData.push({
                                                            name: 'Magnet Joule Loss',
                                                            value: specPerformanceDict.rated_magnet_Joule_loss
                                                        });
                                                    }
                                                    if (specPerformanceDict.rated_windage_loss !== undefined && specPerformanceDict.rated_windage_loss > 0) {
                                                        lossData.push({
                                                            name: 'Windage Loss',
                                                            value: specPerformanceDict.rated_windage_loss
                                                        });
                                                    }

                                                    return lossData.length > 0 ? (
                                                        <DonutChart data={lossData} title="Loss Distribution" />
                                                    ) : (
                                                        <div className="h-full flex items-center justify-center text-muted-foreground text-sm">
                                                            No loss data available
                                                        </div>
                                                    );
                                                })()}
                                            </div>
                                        </div>
                                    </div>

                                    {/* Cost Breakdown */}
                                    {(specPerformanceDict.Cost_Fe !== undefined || specPerformanceDict.Cost_Cu !== undefined || specPerformanceDict.Cost_PM !== undefined) && (
                                        <div className="bg-card rounded-lg border border-border overflow-hidden">
                                            <div className="bg-muted/50 px-3 py-2 border-b border-border">
                                                <h4 className="text-xs font-semibold text-foreground flex items-center">
                                                    <BarChart3 className="w-3 h-3 mr-1" /> Cost Breakdown
                                                </h4>
                                            </div>
                                            <div className="grid grid-cols-2 md:grid-cols-4 lg:grid-cols-5 gap-2 p-3">
                                                {specPerformanceDict.Cost_Fe !== undefined && (
                                                    <div className="space-y-0.5">
                                                        <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">Iron Cost</div>
                                                        <div className="text-xs font-mono text-foreground">{specPerformanceDict.Cost_Fe.toFixed(2)}</div>
                                                    </div>
                                                )}
                                                {specPerformanceDict.Cost_Cu !== undefined && (
                                                    <div className="space-y-0.5">
                                                        <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">Copper Cost</div>
                                                        <div className="text-xs font-mono text-foreground">{specPerformanceDict.Cost_Cu.toFixed(2)}</div>
                                                    </div>
                                                )}
                                                {specPerformanceDict.Cost_PM !== undefined && (
                                                    <div className="space-y-0.5">
                                                        <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">Magnet Cost</div>
                                                        <div className="text-xs font-mono text-foreground">{specPerformanceDict.Cost_PM.toFixed(2)}</div>
                                                    </div>
                                                )}
                                            </div>
                                        </div>
                                    )}

                                    {/* Other Performance Metrics */}
                                    <div className="bg-card rounded-lg border border-border overflow-hidden">
                                        <div className="bg-muted/50 px-3 py-2 border-b border-border">
                                            <h4 className="text-xs font-semibold text-foreground flex items-center">
                                                <Cpu className="w-3 h-3 mr-1" /> Other Metrics
                                            </h4>
                                        </div>
                                        <div className="grid grid-cols-2 md:grid-cols-4 lg:grid-cols-5 gap-2 p-3">
                                            {specPerformanceDict.power_factor !== undefined && (
                                                <div className="space-y-0.5">
                                                    <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">Power Factor</div>
                                                    <div className="text-xs font-mono text-foreground">{specPerformanceDict.power_factor.toFixed(4)}</div>
                                                </div>
                                            )}
                                            {specPerformanceDict.rated_ratio !== undefined && (
                                                <div className="space-y-0.5">
                                                    <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">Rated Ratio</div>
                                                    <div className="text-xs font-mono text-foreground">{specPerformanceDict.rated_ratio.toFixed(4)}</div>
                                                </div>
                                            )}
                                            {specPerformanceDict.rated_stack_length_mm !== undefined && (
                                                <div className="space-y-0.5">
                                                    <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">Rated Stack Length</div>
                                                    <div className="text-xs font-mono text-foreground">{specPerformanceDict.rated_stack_length_mm.toFixed(2)} <span className="text-muted-foreground">mm</span></div>
                                                </div>
                                            )}
                                            {specPerformanceDict.rotor_weight !== undefined && (
                                                <div className="space-y-1">
                                                    <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">Rotor Weight</div>
                                                    <div className="text-sm font-mono text-foreground">{specPerformanceDict.rotor_weight.toFixed(4)} <span className="text-muted-foreground">kg</span></div>
                                                </div>
                                            )}
                                            {specPerformanceDict.coil_flux_linkage_peak2peak_value !== undefined && (
                                                <div className="space-y-1">
                                                    <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">Coil Flux Linkage (Pk-Pk)</div>
                                                    <div className="text-sm font-mono text-foreground">{specPerformanceDict.coil_flux_linkage_peak2peak_value.toFixed(6)} <span className="text-muted-foreground">Wb</span></div>
                                                </div>
                                            )}
                                            {specPerformanceDict.mm2_slot_area !== undefined && (
                                                <div className="space-y-1">
                                                    <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">Slot Area</div>
                                                    <div className="text-sm font-mono text-foreground">{specPerformanceDict.mm2_slot_area.toFixed(2)} <span className="text-muted-foreground">mm²</span></div>
                                                </div>
                                            )}
                                            {specPerformanceDict.project_name && (
                                                <div className="space-y-1">
                                                    <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">Project Name</div>
                                                    <div className="text-sm font-mono text-foreground">{specPerformanceDict.project_name}</div>
                                                </div>
                                            )}
                                            {specPerformanceDict.individual_name && (
                                                <div className="space-y-1">
                                                    <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">Individual Name</div>
                                                    <div className="text-sm font-mono text-foreground">{specPerformanceDict.individual_name}</div>
                                                </div>
                                            )}
                                            {specPerformanceDict.number_current_generation !== undefined && (
                                                <div className="space-y-1">
                                                    <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">Generation</div>
                                                    <div className="text-sm font-mono text-foreground">{specPerformanceDict.number_current_generation}</div>
                                                </div>
                                            )}
                                            {specPerformanceDict.individual_index !== undefined && (
                                                <div className="space-y-1">
                                                    <div className="text-xs text-muted-foreground uppercase tracking-wider font-semibold">Individual Index</div>
                                                    <div className="text-sm font-mono text-foreground">{specPerformanceDict.individual_index}</div>
                                                </div>
                                            )}
                                        </div>
                                    </div>
                                </div>
                            </div>
                        )}

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
