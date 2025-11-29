"use client";

import React, { useState, useEffect } from 'react';
import { Activity, Settings, Zap, Cpu, Maximize, RotateCw, Wind, BarChart3, Loader2, Sparkles, AlertCircle, Download } from 'lucide-react';
import { DesignSpecs, OptimizationResult } from '../../types';
import { calculateMachineDesign, generateEfficiencyCurve } from '../../services/physicsEngine';
import { suggestSpecsFromDescription, analyzeDesignResult } from '../../services/geminiService';
import MotorVisualizer from '../../components/MotorVisualizer';
import LinearMachineView from '../../components/LinearMachineView';
import { EfficiencyChart } from '../../components/Charts';
import { InputGroup, KpiCard, ResultRow } from '../../components/DesignHelpers';
import axios from 'axios';
import CsvVisualizer from '../../components/CsvVisualizer';

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

    // View mode: 'circular' or 'linear'
    const [viewMode, setViewMode] = useState<'circular' | 'linear'>('circular');

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
                    </div>
                </aside>

                {/* Center: Visualization & Results */}
                <div className="flex-1 bg-background flex flex-col h-full overflow-hidden">

                    <div className="flex-1 overflow-y-auto p-6">
                        <div className="grid grid-cols-1 xl:grid-cols-2 gap-6 mb-6">

                            {/* Visualization Panel */}
                            <div className="flex flex-col space-y-4">
                                <div className="flex items-center justify-between">
                                    <h3 className="text-lg font-medium text-foreground">Geometry</h3>
                                    <div className="flex items-center space-x-2">
                                        <button
                                            onClick={() => setViewMode('circular')}
                                            className={`text-xs px-3 py-1 rounded transition-colors ${viewMode === 'circular'
                                                    ? 'bg-primary text-primary-foreground'
                                                    : 'bg-muted text-muted-foreground hover:bg-muted/80'
                                                }`}
                                        >
                                            Circular
                                        </button>
                                        <button
                                            onClick={() => setViewMode('linear')}
                                            className={`text-xs px-3 py-1 rounded transition-colors ${viewMode === 'linear'
                                                    ? 'bg-primary text-primary-foreground'
                                                    : 'bg-muted text-muted-foreground hover:bg-muted/80'
                                                }`}
                                        >
                                            Linear
                                        </button>
                                    </div>
                                </div>
                                {result && (
                                    viewMode === 'circular'
                                        ? <MotorVisualizer geometry={result.geometry} />
                                        : <LinearMachineView geometry={result.geometry} />
                                )}

                                {/* Stats Summary under Viz */}
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

                            {/* Performance Panel */}
                            <div className="flex flex-col space-y-4">
                                <div className="flex items-center justify-between">
                                    <h3 className="text-lg font-medium text-foreground">Performance</h3>
                                    <button
                                        onClick={handleAiAnalyze}
                                        disabled={isAnalyzing}
                                        className="text-xs flex items-center text-primary hover:text-primary/80 transition-colors"
                                    >
                                        {isAnalyzing ? "Analyzing..." : "AI Analysis"} <Sparkles className="w-3 h-3 ml-1" />
                                    </button>
                                </div>

                                {/* KPI Cards */}
                                <div className="grid grid-cols-2 gap-4">
                                    <KpiCard label="Efficiency" value={result?.performance.efficiency.toFixed(1)} unit="%" icon={<Wind className="w-4 h-4 text-emerald-500" />} />
                                    <KpiCard label="Rated Torque" value={result?.performance.torque.toFixed(2)} unit="Nm" icon={<RotateCw className="w-4 h-4 text-amber-500" />} />
                                    <KpiCard label="Suspension Force" value={result?.performance.suspensionForce.toFixed(1)} unit="N" icon={<Activity className="w-4 h-4 text-purple-500" />} />
                                    <KpiCard label="Cost Est." value={result?.performance.materialCost.toFixed(2)} unit="$" icon={<BarChart3 className="w-4 h-4 text-blue-500" />} />
                                </div>

                                {/* Charts */}
                                <EfficiencyChart data={efficiencyData} />

                                {/* AI Analysis Box */}
                                {aiAnalysis && (
                                    <div className="bg-primary/10/20 border border-primary/30 p-4 rounded-lg mt-4 animate-in fade-in slide-in-from-bottom-2 duration-500">
                                        <h4 className="text-primary text-xs font-bold uppercase mb-2 flex items-center">
                                            <Cpu className="w-3 h-3 mr-1" /> Engineer's Note
                                        </h4>
                                        <p className="text-sm text-primary-foreground leading-relaxed">
                                            {aiAnalysis}
                                        </p>
                                    </div>
                                )}
                            </div>
                        </div>

                        {/* CSV Visualizer Section */}
                        <div className="mb-8">
                            <h3 className="text-lg font-medium text-foreground mb-4">Simulation Results (CSV)</h3>
                            <div className="h-[500px]">
                                <CsvVisualizer projectName="SuperCoolPMSM" />
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
