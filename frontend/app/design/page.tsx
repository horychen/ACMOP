"use client"

import { useState, useEffect } from "react"
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card"
import { Button } from "@/components/ui/button"
import { Input } from "@/components/ui/input"
import { Label } from "@/components/ui/label"
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select"
import { Switch } from "@/components/ui/switch"
import { ProjectSelector } from "@/components/ProjectSelector"
import { useProject } from "@/context/ProjectContext"
import { Loader2, Settings, Play, Copy, Check, ChevronDown, ChevronUp } from "lucide-react"
import axios from "axios"
import { useTheme } from "@/context/ThemeContext"
import MotorVisualizer from "@/components/MotorVisualizer"
import InteractiveDrawingEditor from "@/components/InteractiveDrawingEditor"

interface DesignConfig {
    select_spec: string
    select_fea_config_dict: string
    project_loc: string
    bool_show_GUI: boolean
}

export default function DesignViewerPage() {
    const { currentProject } = useProject()
    const { theme } = useTheme()
    const [availableSpecs, setAvailableSpecs] = useState<string[]>([])
    const [availableFEAConfigs, setAvailableFEAConfigs] = useState<string[]>([])
    const [loading, setLoading] = useState(false)
    const [loadingOptions, setLoadingOptions] = useState(true)
    const [copied, setCopied] = useState(false)
    const [copiedSpec, setCopiedSpec] = useState(false)
    const [copiedFea, setCopiedFea] = useState(false)
    const [specDetails, setSpecDetails] = useState<any>(null)
    const [feaConfigDetails, setFeaConfigDetails] = useState<any>(null)
    const [loadingSpecDetails, setLoadingSpecDetails] = useState(false)
    const [loadingFeaDetails, setLoadingFeaDetails] = useState(false)
    const [expandedSpec, setExpandedSpec] = useState(false)
    const [expandedFea, setExpandedFea] = useState(false)
    const [geometricParams, setGeometricParams] = useState<any>(null)
    const [loadingGeometricParams, setLoadingGeometricParams] = useState(false)
    const [selectedComponent, setSelectedComponent] = useState<string>("rotorCore")

    // Helper function to render spec details in a structured way
    const renderSpecDetails = (spec: any) => {
        if (!spec) return null

        return (
            <div className="space-y-6">
                {/* Machine Topology */}
                <div>
                    <h4 className="text-sm font-semibold mb-3 text-foreground">Machine Topology</h4>
                    <div className="grid grid-cols-2 gap-3 text-sm">
                        <div className="flex justify-between">
                            <span className="text-muted-foreground">Machine Type:</span>
                            <span className="font-mono font-medium">{spec.machine_type || "N/A"}</span>
                        </div>
                        <div className="flex justify-between">
                            <span className="text-muted-foreground">Slots (Qs):</span>
                            <span className="font-mono font-medium">{spec.Qs || "N/A"}</span>
                        </div>
                        <div className="flex justify-between">
                            <span className="text-muted-foreground">Poles (p):</span>
                            <span className="font-mono font-medium">{spec.p || "N/A"}</span>
                        </div>
                        <div className="flex justify-between">
                            <span className="text-muted-foreground">Pole Pairs (ps):</span>
                            <span className="font-mono font-medium">{spec.ps || "N/A"}</span>
                        </div>
                        <div className="flex justify-between">
                            <span className="text-muted-foreground">Coil Pitch (y):</span>
                            <span className="font-mono font-medium">{spec.coil_pitch_y || "N/A"}</span>
                        </div>
                        <div className="flex justify-between">
                            <span className="text-muted-foreground">Phases (m):</span>
                            <span className="font-mono font-medium">{spec.m || "N/A"}</span>
                        </div>
                        <div className="flex justify-between">
                            <span className="text-muted-foreground">DPNV/Separate:</span>
                            <span className="font-mono font-medium">{spec.DPNV_or_SEPA ? "DPNV" : "Separate"}</span>
                        </div>
                        <div className="flex justify-between">
                            <span className="text-muted-foreground">Rated Power:</span>
                            <span className="font-mono font-medium">{spec.mec_power ? `${spec.mec_power} W` : "N/A"}</span>
                        </div>
                    </div>
                </div>

                {/* Circuit Excitation (EX-user) */}
                {spec["EX-user"] && (
                    <div>
                        <h4 className="text-sm font-semibold mb-3 text-foreground">Circuit Excitation</h4>
                        <div className="grid grid-cols-2 gap-3 text-sm">
                            <div className="flex justify-between">
                                <span className="text-muted-foreground">Frequency (Hz):</span>
                                <span className="font-mono font-medium">{spec["EX-user"].ExcitationFreqSimulated || "N/A"}</span>
                            </div>
                            <div className="flex justify-between">
                                <span className="text-muted-foreground">Voltage (V):</span>
                                <span className="font-mono font-medium">{spec["EX-user"].VoltageRating || "N/A"}</span>
                            </div>
                            <div className="flex justify-between">
                                <span className="text-muted-foreground">Stack Length (mm):</span>
                                <span className="font-mono font-medium">{spec["EX-user"].mm_stack_length || "N/A"}</span>
                            </div>
                            <div className="flex justify-between">
                                <span className="text-muted-foreground">Steel:</span>
                                <span className="font-mono font-medium">{spec["EX-user"].Steel || "N/A"}</span>
                            </div>
                            <div className="flex justify-between">
                                <span className="text-muted-foreground">Current Density:</span>
                                <span className="font-mono font-medium">{spec["EX-user"].Js ? `${(spec["EX-user"].Js / 1e6).toFixed(1)} A/mm²` : "N/A"}</span>
                            </div>
                            <div className="flex justify-between">
                                <span className="text-muted-foreground">Temperature (°C):</span>
                                <span className="font-mono font-medium">{spec["EX-user"].Temperature || "N/A"}</span>
                            </div>
                            <div className="flex justify-between">
                                <span className="text-muted-foreground">Winding Fill:</span>
                                <span className="font-mono font-medium">{spec["EX-user"].WindingFill ? `${(spec["EX-user"].WindingFill * 100).toFixed(1)}%` : "N/A"}</span>
                            </div>
                            <div className="flex justify-between">
                                <span className="text-muted-foreground">AWG:</span>
                                <span className="font-mono font-medium">{spec["EX-user"].AWG || "N/A"}</span>
                            </div>
                        </div>
                    </div>
                )}

                {/* Geometric Parameters (GP-user) */}
                {spec["GP-user"] && (
                    <div>
                        <h4 className="text-sm font-semibold mb-3 text-foreground">Geometric Parameters</h4>
                        <div className="space-y-2">
                            {Object.entries(spec["GP-user"]).map(([key, param]: [string, any]) => {
                                const paramType = param.type || "unknown"
                                const getTypeColor = (type: string) => {
                                    if (type === "fixed") {
                                        return theme === 'dark' 
                                            ? "bg-blue-950/30 border-blue-700/50 text-blue-200" 
                                            : "bg-blue-50 border-blue-200 text-blue-900"
                                    } else if (type === "free") {
                                        return theme === 'dark' 
                                            ? "bg-green-950/30 border-green-700/50 text-green-200" 
                                            : "bg-green-50 border-green-200 text-green-900"
                                    } else if (type === "derived") {
                                        return theme === 'dark' 
                                            ? "bg-purple-950/30 border-purple-700/50 text-purple-200" 
                                            : "bg-purple-50 border-purple-200 text-purple-900"
                                    }
                                    return theme === 'dark' 
                                        ? "bg-slate-800/30 border-slate-700 text-slate-300" 
                                        : "bg-slate-50 border-slate-300 text-slate-900"
                                }
                                
                                const getTypeBadgeColor = (type: string) => {
                                    if (type === "fixed") {
                                        return theme === 'dark' 
                                            ? "bg-blue-700 text-blue-100" 
                                            : "bg-blue-500 text-white"
                                    } else if (type === "free") {
                                        return theme === 'dark' 
                                            ? "bg-green-700 text-green-100" 
                                            : "bg-green-500 text-white"
                                    } else if (type === "derived") {
                                        return theme === 'dark' 
                                            ? "bg-purple-700 text-purple-100" 
                                            : "bg-purple-500 text-white"
                                    }
                                    return theme === 'dark' 
                                        ? "bg-slate-700 text-slate-100" 
                                        : "bg-slate-500 text-white"
                                }
                                
                                return (
                                    <div key={key} className={`flex items-center justify-between p-2.5 rounded border text-sm ${getTypeColor(paramType)}`}>
                                        <div className="flex items-center gap-2 flex-1">
                                            <span className="font-medium">{key}:</span>
                                            <span className={`px-2 py-0.5 rounded text-xs font-semibold ${getTypeBadgeColor(paramType)}`}>
                                                {paramType}
                                            </span>
                                        </div>
                                        <div className="flex items-center gap-3">
                                            {param.value !== null && param.value !== undefined ? (
                                                <span className="font-mono font-medium">
                                                    {typeof param.value === 'number' ? param.value.toFixed(3).replace(/\.?0+$/, '') : param.value}
                                                    {key.includes('mm_') && ' mm'}
                                                </span>
                                            ) : (
                                                <span className="italic opacity-70">derived</span>
                                            )}
                                            {param.bounds && Array.isArray(param.bounds) && param.bounds[0] !== null && param.bounds[1] !== null && (
                                                <span className="text-xs opacity-70">
                                                    [{param.bounds[0]}, {param.bounds[1]}]
                                                </span>
                                            )}
                                        </div>
                                    </div>
                                )
                            })}
                        </div>
                        {/* Legend */}
                        <div className="mt-3 pt-3 border-t flex flex-wrap gap-3 text-xs">
                            <div className="flex items-center gap-1.5">
                                <div className={`w-3 h-3 rounded ${theme === 'dark' ? 'bg-blue-700' : 'bg-blue-500'}`}></div>
                                <span className="text-muted-foreground">Fixed</span>
                            </div>
                            <div className="flex items-center gap-1.5">
                                <div className={`w-3 h-3 rounded ${theme === 'dark' ? 'bg-green-700' : 'bg-green-500'}`}></div>
                                <span className="text-muted-foreground">Free</span>
                            </div>
                            <div className="flex items-center gap-1.5">
                                <div className={`w-3 h-3 rounded ${theme === 'dark' ? 'bg-purple-700' : 'bg-purple-500'}`}></div>
                                <span className="text-muted-foreground">Derived</span>
                            </div>
                        </div>
                    </div>
                )}

                {/* Initial Guesses */}
                <div>
                    <h4 className="text-sm font-semibold mb-3 text-foreground">Initial Design Guesses</h4>
                    <div className="grid grid-cols-2 gap-3 text-sm">
                        {spec.guess_air_gap_flux_density_Bg && (
                            <div className="flex justify-between">
                                <span className="text-muted-foreground">Air Gap Flux (T):</span>
                                <span className="font-mono font-medium">{spec.guess_air_gap_flux_density_Bg}</span>
                            </div>
                        )}
                        {spec.guess_stator_tooth_flux_density_Bst && (
                            <div className="flex justify-between">
                                <span className="text-muted-foreground">Tooth Flux (T):</span>
                                <span className="font-mono font-medium">{spec.guess_stator_tooth_flux_density_Bst}</span>
                            </div>
                        )}
                        {spec.guess_efficiency && (
                            <div className="flex justify-between">
                                <span className="text-muted-foreground">Efficiency:</span>
                                <span className="font-mono font-medium">{(spec.guess_efficiency * 100).toFixed(1)}%</span>
                            </div>
                        )}
                        {spec.guess_power_factor && (
                            <div className="flex justify-between">
                                <span className="text-muted-foreground">Power Factor:</span>
                                <span className="font-mono font-medium">{spec.guess_power_factor}</span>
                            </div>
                        )}
                    </div>
                </div>
            </div>
        )
    }

    // Helper function to render FEA config details in a structured way
    const renderFeaConfigDetails = (config: any) => {
        if (!config) return null

        return (
            <div className="space-y-6">
                {/* Circuit Configuration */}
                {(config["circuit.TORQUE_CURRENT_RATIO"] !== undefined || config["circuit.SUSPENSION_CURRENT_RATIO"] !== undefined) && (
                    <div>
                        <h4 className="text-sm font-semibold mb-3 text-foreground">Circuit Configuration</h4>
                        <div className="grid grid-cols-2 gap-3 text-sm">
                            {config["circuit.TORQUE_CURRENT_RATIO"] !== undefined && (
                                <div className="flex justify-between">
                                    <span className="text-muted-foreground">Torque Current Ratio:</span>
                                    <span className="font-mono font-medium">{config["circuit.TORQUE_CURRENT_RATIO"]}</span>
                                </div>
                            )}
                            {config["circuit.SUSPENSION_CURRENT_RATIO"] !== undefined && (
                                <div className="flex justify-between">
                                    <span className="text-muted-foreground">Suspension Current Ratio:</span>
                                    <span className="font-mono font-medium">{config["circuit.SUSPENSION_CURRENT_RATIO"]}</span>
                                </div>
                            )}
                        </div>
                    </div>
                )}

                {/* FEA Simulation Settings */}
                {Object.keys(config).some(k => k.startsWith("designer.")) && (
                    <div>
                        <h4 className="text-sm font-semibold mb-3 text-foreground">FEA Simulation Settings</h4>
                        <div className="space-y-2">
                            {/* General Settings */}
                            {Object.entries(config)
                                .filter(([key]) => key.startsWith("designer.") && 
                                    !key.includes("number_of_steps") && 
                                    !key.includes("number_cycles") && 
                                    !key.includes("StepPerCycle") &&
                                    !key.includes("TranRef"))
                                .map(([key, value]) => {
                                    const displayKey = key.replace("designer.", "").replace(/_/g, " ")
                                    return (
                                        <div key={key} className="flex items-center justify-between p-2 rounded border bg-card text-sm">
                                            <span className="text-muted-foreground capitalize">{displayKey}:</span>
                                            <span className="font-mono font-medium">
                                                {typeof value === 'boolean' ? (value ? 'Yes' : 'No') : String(value)}
                                            </span>
                                        </div>
                                    )
                                })}
                            
                            {/* Time Stepping Settings */}
                            {(Object.keys(config).some(k => k.includes("number_of_steps") || k.includes("number_cycles") || k.includes("StepPerCycle"))) && (
                                <div className="mt-3 pt-3 border-t">
                                    <h5 className="text-xs font-semibold mb-2 text-muted-foreground uppercase">Time Stepping</h5>
                                    {Object.entries(config)
                                        .filter(([key]) => key.startsWith("designer.") && 
                                            (key.includes("number_of_steps") || 
                                             key.includes("number_cycles") || 
                                             key.includes("StepPerCycle") ||
                                             key.includes("TranRef")))
                                        .map(([key, value]) => {
                                            const displayKey = key.replace("designer.", "").replace(/_/g, " ")
                                            return (
                                                <div key={key} className="flex items-center justify-between p-2 rounded border bg-card text-sm">
                                                    <span className="text-muted-foreground capitalize">{displayKey}:</span>
                                                    <span className="font-mono font-medium">{String(value)}</span>
                                                </div>
                                            )
                                        })}
                                </div>
                            )}
                        </div>
                    </div>
                )}

                {/* Optimization Settings */}
                {Object.keys(config).some(k => k.startsWith("moo.")) && (
                    <div>
                        <h4 className="text-sm font-semibold mb-3 text-foreground">Optimization Settings</h4>
                        <div className="space-y-2">
                            {Object.entries(config)
                                .filter(([key]) => key.startsWith("moo."))
                                .map(([key, value]) => {
                                    const displayKey = key.replace("moo.", "")
                                    return (
                                        <div key={key} className="flex items-center justify-between p-2 rounded border bg-card text-sm">
                                            <span className="text-muted-foreground">{displayKey}:</span>
                                            <span className="font-mono font-medium">
                                                {value === null ? 'null' : typeof value === 'boolean' ? (value ? 'Yes' : 'No') : String(value)}
                                            </span>
                                        </div>
                                    )
                                })}
                        </div>
                    </div>
                )}

                {/* Other Settings */}
                {Object.keys(config).some(k => !k.startsWith("designer.") && !k.startsWith("moo.") && !k.startsWith("circuit.") && !k.startsWith("femm.")) && (
                    <div>
                        <h4 className="text-sm font-semibold mb-3 text-foreground">Other Settings</h4>
                        <div className="space-y-2">
                            {Object.entries(config)
                                .filter(([key]) => !key.startsWith("designer.") && !key.startsWith("moo.") && !key.startsWith("circuit.") && !key.startsWith("femm."))
                                .map(([key, value]) => (
                                    <div key={key} className="flex items-center justify-between p-2 rounded border bg-card text-sm">
                                        <span className="text-muted-foreground">{key}:</span>
                                        <span className="font-mono font-medium">
                                            {value === null ? 'null' : typeof value === 'boolean' ? (value ? 'Yes' : 'No') : String(value)}
                                        </span>
                                    </div>
                                ))}
                        </div>
                    </div>
                )}
            </div>
        )
    }
    const [config, setConfig] = useState<DesignConfig>({
        select_spec: "",
        select_fea_config_dict: "",
        project_loc: "../_default/",
        bool_show_GUI: true
    })

    // Load available options
    useEffect(() => {
        const loadOptions = async () => {
            setLoadingOptions(true)
            try {
                const [specsRes, feaRes] = await Promise.all([
                    axios.get('http://localhost:8000/api/design/specs'),
                    axios.get('http://localhost:8000/api/design/fea-configs')
                ])
                setAvailableSpecs(specsRes.data.specs || [])
                setAvailableFEAConfigs(feaRes.data.fea_configs || [])
                
                // Set defaults if available
                if (specsRes.data.specs && specsRes.data.specs.length > 0 && !config.select_spec) {
                    setConfig(prev => ({ ...prev, select_spec: specsRes.data.specs[0] }))
                }
                if (feaRes.data.fea_configs && feaRes.data.fea_configs.length > 0 && !config.select_fea_config_dict) {
                    setConfig(prev => ({ ...prev, select_fea_config_dict: feaRes.data.fea_configs[0] }))
                }
            } catch (error) {
                console.error("Failed to load options:", error)
                // Show user-friendly error message
                if (axios.isAxiosError(error)) {
                    console.error("API Error:", error.response?.data || error.message)
                }
            } finally {
                setLoadingOptions(false)
            }
        }
        loadOptions()
        // eslint-disable-next-line react-hooks/exhaustive-deps
    }, [])

    // Load detailed spec when selected
    useEffect(() => {
        const loadSpecDetails = async () => {
            if (!config.select_spec) {
                setSpecDetails(null)
                setGeometricParams(null)
                return
            }
            setLoadingSpecDetails(true)
            try {
                const [specResponse, geoResponse] = await Promise.all([
                    axios.get(`http://localhost:8000/api/design/specs/${encodeURIComponent(config.select_spec)}`),
                    axios.get(`http://localhost:8000/api/design/visualization/geometric-params?spec_name=${encodeURIComponent(config.select_spec)}`)
                ])
                setSpecDetails(specResponse.data)
                setGeometricParams(geoResponse.data)
            } catch (error) {
                console.error("Failed to load spec details:", error)
                setSpecDetails(null)
                setGeometricParams(null)
            } finally {
                setLoadingSpecDetails(false)
            }
        }
        loadSpecDetails()
    }, [config.select_spec])

    // Load detailed FEA config when selected
    useEffect(() => {
        const loadFeaConfigDetails = async () => {
            if (!config.select_fea_config_dict) {
                setFeaConfigDetails(null)
                return
            }
            setLoadingFeaDetails(true)
            try {
                const response = await axios.get(`http://localhost:8000/api/design/fea-configs/${encodeURIComponent(config.select_fea_config_dict)}`)
                setFeaConfigDetails(response.data)
            } catch (error) {
                console.error("Failed to load FEA config details:", error)
                setFeaConfigDetails(null)
            } finally {
                setLoadingFeaDetails(false)
            }
        }
        loadFeaConfigDetails()
    }, [config.select_fea_config_dict])

    const handleInitialize = async () => {
        if (!config.select_spec || !config.select_fea_config_dict) {
            alert("Please select both Specification and FEA Configuration")
            return
        }

        setLoading(true)
        try {
            const response = await axios.post('http://localhost:8000/api/design/initialize', config)
            console.log("Design initialized:", response.data)
            alert("Design initialized successfully!")
            // TODO: Handle success - maybe show a success message or redirect
        } catch (error) {
            console.error("Failed to initialize design:", error)
            if (axios.isAxiosError(error)) {
                const errorMsg = error.response?.data?.detail || error.message
                alert(`Failed to initialize design: ${errorMsg}`)
            } else {
                alert("Failed to initialize design. Please check the console for details.")
            }
        } finally {
            setLoading(false)
        }
    }

    return (
        <div className="space-y-6">
            <div className="flex items-center justify-between">
                <h1 className="text-3xl font-bold tracking-tight">Design Viewer</h1>
                <ProjectSelector />
            </div>

            <div className="grid gap-4 md:grid-cols-2">
                {/* Configuration Form */}
                <Card className="col-span-2 md:col-span-1">
                    <CardHeader>
                        <CardTitle className="flex items-center gap-2">
                            <Settings className="w-5 h-5" />
                            Design Configuration
                        </CardTitle>
                        <CardDescription>
                            Configure the design process parameters
                        </CardDescription>
                    </CardHeader>
                    <CardContent className="space-y-4">
                        {loadingOptions ? (
                            <div className="flex items-center justify-center py-8">
                                <Loader2 className="w-6 h-6 animate-spin text-muted-foreground" />
                            </div>
                        ) : (
                            <>
                                {/* Select Specification */}
                                <div className="space-y-2">
                                    <Label htmlFor="select_spec">Machine Specification</Label>
                                    <Select
                                        value={config.select_spec}
                                        onValueChange={(value) => setConfig(prev => ({ ...prev, select_spec: value }))}
                                    >
                                        <SelectTrigger id="select_spec">
                                            <SelectValue placeholder="Select specification" />
                                        </SelectTrigger>
                                        <SelectContent>
                                            {availableSpecs.map((spec) => (
                                                <SelectItem key={spec} value={spec}>
                                                    {spec}
                                                </SelectItem>
                                            ))}
                                        </SelectContent>
                                    </Select>
                                </div>

                                {/* Select FEA Configuration */}
                                <div className="space-y-2">
                                    <Label htmlFor="select_fea_config">FEA Configuration</Label>
                                    <Select
                                        value={config.select_fea_config_dict}
                                        onValueChange={(value) => setConfig(prev => ({ ...prev, select_fea_config_dict: value }))}
                                    >
                                        <SelectTrigger id="select_fea_config">
                                            <SelectValue placeholder="Select FEA configuration" />
                                        </SelectTrigger>
                                        <SelectContent>
                                            {availableFEAConfigs.map((fea) => (
                                                <SelectItem key={fea} value={fea}>
                                                    {fea}
                                                </SelectItem>
                                            ))}
                                        </SelectContent>
                                    </Select>
                                </div>

                                {/* Project Location */}
                                <div className="space-y-2">
                                    <Label htmlFor="project_loc">Project Location</Label>
                                    <Input
                                        id="project_loc"
                                        value={config.project_loc}
                                        onChange={(e) => setConfig(prev => ({ ...prev, project_loc: e.target.value }))}
                                        placeholder="../_default/"
                                    />
                                </div>

                                {/* Show GUI Toggle */}
                                <div className="flex items-center justify-between space-x-2">
                                    <Label htmlFor="bool_show_GUI" className="flex-1">
                                        Show GUI
                                    </Label>
                                    <Switch
                                        id="bool_show_GUI"
                                        checked={config.bool_show_GUI}
                                        onCheckedChange={(checked) => setConfig(prev => ({ ...prev, bool_show_GUI: checked }))}
                                    />
                                </div>

                                {/* Initialize Button */}
                                <Button
                                    onClick={handleInitialize}
                                    disabled={loading || !config.select_spec || !config.select_fea_config_dict}
                                    className="w-full"
                                >
                                    {loading ? (
                                        <>
                                            <Loader2 className="w-4 h-4 mr-2 animate-spin" />
                                            Initializing...
                                        </>
                                    ) : (
                                        <>
                                            <Play className="w-4 h-4 mr-2" />
                                            Initialize Design
                                        </>
                                    )}
                                </Button>
                            </>
                        )}
                    </CardContent>
                </Card>

                {/* Preview/Status Card */}
                <Card className="col-span-2 md:col-span-1">
                    <CardHeader>
                        <CardTitle>Configuration Summary</CardTitle>
                        <CardDescription>Review your design configuration</CardDescription>
                    </CardHeader>
                    <CardContent>
                        <div className="space-y-3">
                            <div>
                                <span className="text-sm font-medium text-muted-foreground">Specification:</span>
                                <p className="text-sm mt-1">{config.select_spec || "Not selected"}</p>
                            </div>
                            <div>
                                <span className="text-sm font-medium text-muted-foreground">FEA Config:</span>
                                <p className="text-sm mt-1">{config.select_fea_config_dict || "Not selected"}</p>
                            </div>
                            <div>
                                <span className="text-sm font-medium text-muted-foreground">Project Location:</span>
                                <p className="text-sm mt-1 font-mono">{config.project_loc}</p>
                            </div>
                            <div>
                                <span className="text-sm font-medium text-muted-foreground">Show GUI:</span>
                                <p className="text-sm mt-1">{config.bool_show_GUI ? "Yes" : "No"}</p>
                            </div>
                        </div>
                    </CardContent>
                </Card>
            </div>

            {/* Machine Visualization */}
            {geometricParams && (
                <Card>
                    <CardHeader>
                        <CardTitle>Machine Geometry Visualization</CardTitle>
                        <CardDescription>Visual representation of the machine based on geometric parameters</CardDescription>
                    </CardHeader>
                    <CardContent>
                        <div className="grid gap-6 md:grid-cols-2">
                            {/* Visualization */}
                            <div className="space-y-4">
                                <h4 className="text-sm font-semibold">Cross-Section View</h4>
                                {geometricParams.geometric_parameters && (() => {
                                    const gp = geometricParams.geometric_parameters
                                    // Convert GP parameters to MachineGeometry format
                                    const geometry = {
                                        statorOuterRadius: gp.mm_r_so?.value || 6.5,
                                        statorInnerRadius: gp.mm_r_si?.value || (gp.mm_r_so?.value || 6.5) * (gp.split_ratio?.value || 0.6),
                                        rotorOuterRadius: gp.mm_r_ro?.value || 3.0,
                                        rotorInnerRadius: (gp.mm_r_ro?.value || 3.0) * 0.3, // Estimate shaft radius
                                        slotDepth: (gp.mm_d_st?.value || 0) + (gp.mm_d_sts?.value || 0),
                                        toothWidth: gp.mm_w_st?.value || 1.0,
                                        magnetThickness: gp.mm_d_pm?.value || 1.0,
                                        airGap: gp.mm_d_mech_air_gap?.value || 0.3,
                                        slots: geometricParams.slots || 12,
                                        poles: geometricParams.poles || 4
                                    }
                                    return <MotorVisualizer geometry={geometry} />
                                })()}
                            </div>
                            
                            {/* Geometric Parameters Used for Drawing */}
                            <div className="space-y-4">
                                <h4 className="text-sm font-semibold">Geometric Parameters for Drawing</h4>
                                <div className="space-y-2 max-h-[400px] overflow-y-auto">
                                    {geometricParams.geometric_parameters && Object.entries(geometricParams.geometric_parameters).map(([key, param]: [string, any]) => {
                                        const paramType = param.type || "unknown"
                                        const getTypeColor = (type: string) => {
                                            if (type === "fixed") {
                                                return theme === 'dark' 
                                                    ? "bg-blue-950/30 border-blue-700/50 text-blue-200" 
                                                    : "bg-blue-50 border-blue-200 text-blue-900"
                                            } else if (type === "free") {
                                                return theme === 'dark' 
                                                    ? "bg-green-950/30 border-green-700/50 text-green-200" 
                                                    : "bg-green-50 border-green-200 text-green-900"
                                            } else if (type === "derived") {
                                                return theme === 'dark' 
                                                    ? "bg-purple-950/30 border-purple-700/50 text-purple-200" 
                                                    : "bg-purple-50 border-purple-200 text-purple-900"
                                            }
                                            return theme === 'dark' 
                                                ? "bg-slate-800/30 border-slate-700 text-slate-300" 
                                                : "bg-slate-50 border-slate-300 text-slate-900"
                                        }
                                        
                                        const getTypeBadgeColor = (type: string) => {
                                            if (type === "fixed") {
                                                return theme === 'dark' 
                                                    ? "bg-blue-700 text-blue-100" 
                                                    : "bg-blue-500 text-white"
                                            } else if (type === "free") {
                                                return theme === 'dark' 
                                                    ? "bg-green-700 text-green-100" 
                                                    : "bg-green-500 text-white"
                                            } else if (type === "derived") {
                                                return theme === 'dark' 
                                                    ? "bg-purple-700 text-purple-100" 
                                                    : "bg-purple-500 text-white"
                                            }
                                            return theme === 'dark' 
                                                ? "bg-slate-700 text-slate-100" 
                                                : "bg-slate-500 text-white"
                                        }
                                        
                                        return (
                                            <div key={key} className={`flex items-center justify-between p-2.5 rounded border text-sm ${getTypeColor(paramType)}`}>
                                                <div className="flex items-center gap-2 flex-1">
                                                    <span className="font-medium">{key}:</span>
                                                    <span className={`px-2 py-0.5 rounded text-xs font-semibold ${getTypeBadgeColor(paramType)}`}>
                                                        {paramType}
                                                    </span>
                                                </div>
                                                <div className="flex items-center gap-3">
                                                    {param.value !== null && param.value !== undefined ? (
                                                        <span className="font-mono font-medium">
                                                            {typeof param.value === 'number' ? param.value.toFixed(3).replace(/\.?0+$/, '') : param.value}
                                                            {key.includes('mm_') && ' mm'}
                                                        </span>
                                                    ) : (
                                                        <span className="italic opacity-70">derived</span>
                                                    )}
                                                    {param.bounds && Array.isArray(param.bounds) && param.bounds[0] !== null && param.bounds[1] !== null && (
                                                        <span className="text-xs opacity-70">
                                                            [{param.bounds[0]}, {param.bounds[1]}]
                                                        </span>
                                                    )}
                                                </div>
                                            </div>
                                        )
                                    })}
                                </div>
                            </div>
                        </div>
                    </CardContent>
                </Card>
            )}

            {/* Interactive Drawing Editor */}
            {specDetails && specDetails.GeometricComponents && (
                <InteractiveDrawingEditor
                    availableComponents={specDetails.GeometricComponents}
                    selectedComponent={selectedComponent}
                    onComponentChange={setSelectedComponent}
                    parameters={{
                        r_ri: geometricParams?.geometric_parameters?.mm_r_ri?.value || 40,
                        d_ri: geometricParams?.geometric_parameters?.mm_d_ri?.value || 8,
                        d_rp: geometricParams?.geometric_parameters?.mm_d_rp?.value || 5,
                        alpha_rp: Math.PI / (specDetails.ps || 2),
                        alpha_rm: (geometricParams?.geometric_parameters?.deg_alpha_rm?.value || 60) * Math.PI / 180,
                        alpha_rs: (geometricParams?.geometric_parameters?.deg_alpha_rs?.value || 10) * Math.PI / 180,
                        r_or: (geometricParams?.geometric_parameters?.mm_r_ro?.value || 40) + 
                              (geometricParams?.geometric_parameters?.mm_d_ri?.value || 8) + 
                              (geometricParams?.geometric_parameters?.mm_d_pm?.value || 2),
                        d_sleeve: geometricParams?.geometric_parameters?.mm_d_sleeve?.value || 1,
                        p: specDetails.ps || 2,
                        r_si: geometricParams?.geometric_parameters?.mm_r_si?.value || 42,
                        r_so: geometricParams?.geometric_parameters?.mm_r_so?.value || 50,
                        Qs: specDetails.Qs || 12,
                        slotDepth: (geometricParams?.geometric_parameters?.mm_d_st?.value || 5) + 
                                  (geometricParams?.geometric_parameters?.mm_d_sts?.value || 0),
                        d_pm: geometricParams?.geometric_parameters?.mm_d_pm?.value || 2
                    }}
                />
            )}

            {/* Detailed Specifications and FEA Config */}
            <div className="grid gap-4 md:grid-cols-2">
                {/* Machine Specification Details */}
                {config.select_spec && (
                    <Card>
                        <CardHeader>
                            <div className="flex items-center justify-between">
                                <div>
                                    <CardTitle>Machine Specification Details</CardTitle>
                                    <CardDescription>{config.select_spec}</CardDescription>
                                </div>
                                <Button
                                    variant="ghost"
                                    size="sm"
                                    onClick={() => setExpandedSpec(!expandedSpec)}
                                >
                                    {expandedSpec ? (
                                        <ChevronUp className="w-4 h-4" />
                                    ) : (
                                        <ChevronDown className="w-4 h-4" />
                                    )}
                                </Button>
                            </div>
                        </CardHeader>
                        {expandedSpec && (
                            <CardContent>
                                {loadingSpecDetails ? (
                                    <div className="flex items-center justify-center py-8">
                                        <Loader2 className="w-6 h-6 animate-spin text-muted-foreground" />
                                    </div>
                                ) : specDetails ? (
                                    <div className="space-y-4">
                                        <div className="relative">
                                            <div className="max-h-[600px] overflow-y-auto pr-2">
                                                {renderSpecDetails(specDetails)}
                                            </div>
                                            <div className="mt-4 pt-4 border-t">
                                                <div className="flex items-center justify-between">
                                                    <span className="text-xs text-muted-foreground">Raw JSON:</span>
                                                    <Button
                                                        variant="ghost"
                                                        size="sm"
                                                        onClick={() => {
                                                            navigator.clipboard.writeText(JSON.stringify(specDetails, null, 2))
                                                            setCopiedSpec(true)
                                                            setTimeout(() => setCopiedSpec(false), 2000)
                                                        }}
                                                    >
                                                        {copiedSpec ? (
                                                            <>
                                                                <Check className="w-3 h-3 mr-1" />
                                                                <span className="text-xs">Copied</span>
                                                            </>
                                                        ) : (
                                                            <>
                                                                <Copy className="w-3 h-3 mr-1" />
                                                                <span className="text-xs">Copy JSON</span>
                                                            </>
                                                        )}
                                                    </Button>
                                                </div>
                                            </div>
                                        </div>
                                    </div>
                                ) : (
                                    <p className="text-sm text-muted-foreground">No details available</p>
                                )}
                            </CardContent>
                        )}
                    </Card>
                )}

                {/* FEA Configuration Details */}
                {config.select_fea_config_dict && (
                    <Card>
                        <CardHeader>
                            <div className="flex items-center justify-between">
                                <div>
                                    <CardTitle>FEA Configuration Details</CardTitle>
                                    <CardDescription>{config.select_fea_config_dict}</CardDescription>
                                </div>
                                <Button
                                    variant="ghost"
                                    size="sm"
                                    onClick={() => setExpandedFea(!expandedFea)}
                                >
                                    {expandedFea ? (
                                        <ChevronUp className="w-4 h-4" />
                                    ) : (
                                        <ChevronDown className="w-4 h-4" />
                                    )}
                                </Button>
                            </div>
                        </CardHeader>
                        {expandedFea && (
                            <CardContent>
                                {loadingFeaDetails ? (
                                    <div className="flex items-center justify-center py-8">
                                        <Loader2 className="w-6 h-6 animate-spin text-muted-foreground" />
                                    </div>
                                ) : feaConfigDetails ? (
                                    <div className="space-y-4">
                                        <div className="relative">
                                            <div className="max-h-[600px] overflow-y-auto pr-2">
                                                {renderFeaConfigDetails(feaConfigDetails)}
                                            </div>
                                            <div className="mt-4 pt-4 border-t">
                                                <div className="flex items-center justify-between">
                                                    <span className="text-xs text-muted-foreground">Raw JSON:</span>
                                                    <Button
                                                        variant="ghost"
                                                        size="sm"
                                                        onClick={() => {
                                                            navigator.clipboard.writeText(JSON.stringify(feaConfigDetails, null, 2))
                                                            setCopiedFea(true)
                                                            setTimeout(() => setCopiedFea(false), 2000)
                                                        }}
                                                    >
                                                        {copiedFea ? (
                                                            <>
                                                                <Check className="w-3 h-3 mr-1" />
                                                                <span className="text-xs">Copied</span>
                                                            </>
                                                        ) : (
                                                            <>
                                                                <Copy className="w-3 h-3 mr-1" />
                                                                <span className="text-xs">Copy JSON</span>
                                                            </>
                                                        )}
                                                    </Button>
                                                </div>
                                            </div>
                                        </div>
                                    </div>
                                ) : (
                                    <p className="text-sm text-muted-foreground">No details available</p>
                                )}
                            </CardContent>
                        )}
                    </Card>
                )}
            </div>

            {/* Terminal Command Card */}
            {config.select_spec && config.select_fea_config_dict && (
                <Card>
                    <CardHeader>
                        <CardTitle>Suggested Terminal Command</CardTitle>
                        <CardDescription>
                            Copy this Python code to initialize AC_Machine_Optiomization_Wrapper with your configuration
                        </CardDescription>
                    </CardHeader>
                    <CardContent>
                        <div className="space-y-4">
                            <div className="relative">
                                <pre className={`p-4 rounded-lg overflow-x-auto text-sm font-mono border ${
                                    theme === 'dark' 
                                        ? 'bg-slate-950 text-slate-100 border-slate-700' 
                                        : 'bg-slate-50 text-slate-900 border-slate-300'
                                }`}>
                                    <code>{`from acmop import AC_Machine_Optiomization_Wrapper

mop = AC_Machine_Optiomization_Wrapper(
    select_spec="${config.select_spec}",
    select_fea_config_dict="${config.select_fea_config_dict}",
    project_loc=r'${config.project_loc}',
    bool_show_GUI=${config.bool_show_GUI}
)`}</code>
                                </pre>
                                <Button
                                    variant="outline"
                                    size="sm"
                                    className="absolute top-2 right-2"
                                    onClick={() => {
                                        const code = `from acmop import AC_Machine_Optiomization_Wrapper

mop = AC_Machine_Optiomization_Wrapper(
    select_spec="${config.select_spec}",
    select_fea_config_dict="${config.select_fea_config_dict}",
    project_loc=r'${config.project_loc}',
    bool_show_GUI=${config.bool_show_GUI}
)`
                                        navigator.clipboard.writeText(code)
                                        setCopied(true)
                                        setTimeout(() => setCopied(false), 2000)
                                    }}
                                >
                                    {copied ? (
                                        <>
                                            <Check className="w-4 h-4 mr-2" />
                                            Copied!
                                        </>
                                    ) : (
                                        <>
                                            <Copy className="w-4 h-4 mr-2" />
                                            Copy
                                        </>
                                    )}
                                </Button>
                            </div>
                            <div className="text-xs text-muted-foreground">
                                <p className="mb-2">To run this code:</p>
                                <ol className="list-decimal list-inside space-y-1 ml-2">
                                    <li>Navigate to the codes4 directory: <code className="bg-muted px-1 py-0.5 rounded">cd backend/codes4</code></li>
                                    <li>Activate your conda environment: <code className="bg-muted px-1 py-0.5 rounded">conda activate acmop</code></li>
                                    <li>Run the code in a Python script or interactive session</li>
                                </ol>
                            </div>
                        </div>
                    </CardContent>
                </Card>
            )}
        </div>
    )
}
export const dynamic = 'force-dynamic'
