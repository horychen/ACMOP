export enum MachineType {
    PERMANENT_MAGNET = 'Permanent Magnet',
    INDUCTION = 'Induction',
    SWITCHED_RELUCTANCE = 'Switched Reluctance'
}

export enum OptimizationObjective {
    MAX_EFFICIENCY = 'Maximize Efficiency',
    MAX_TORQUE_DENSITY = 'Maximize Torque Density',
    MIN_COST = 'Minimize Cost',
    MIN_VOLUME = 'Minimize Volume'
}

export interface DesignSpecs {
    ratedPower: number; // kW
    ratedSpeed: number; // RPM
    ratedVoltage: number; // V
    outerDiameterLimit: number; // mm
    axialLengthLimit: number; // mm
    airGap: number; // mm
    slotCount: number;
    poleCount: number;
    currentDensity: number; // A/mm^2
}

export interface MachineGeometry {
    statorOuterRadius: number;
    statorInnerRadius: number;
    rotorOuterRadius: number;
    rotorInnerRadius: number; // Shaft
    slotDepth: number;
    toothWidth: number;
    magnetThickness: number;
    airGap: number;
    slots: number;
    poles: number;
}

export interface PerformanceMetrics {
    efficiency: number; // %
    torque: number; // Nm
    suspensionForce: number; // N (Bearingless specific)
    copperLoss: number; // W
    ironLoss: number; // W
    powerFactor: number;
    torqueRipple: number; // %
    materialCost: number; // USD
}

export interface OptimizationResult {
    id: string;
    timestamp: number;
    geometry: MachineGeometry;
    performance: PerformanceMetrics;
    specs: DesignSpecs;
    aiAnalysis?: string;
}

export interface ChartPoint {
    x: number;
    y: number;
    label?: string;
}
