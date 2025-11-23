import { DesignSpecs, MachineGeometry, PerformanceMetrics, OptimizationResult } from '../types';

/**
 * Mocks the complex finite element analysis/optimization logic from the original Python project.
 * It provides deterministic, physically plausible outputs based on inputs to drive the visualization.
 */
export const calculateMachineDesign = (specs: DesignSpecs): OptimizationResult => {
    // 1. Basic Geometry Derivation
    // In a real scenario, this would be an iterative solver.

    const outerRadius = specs.outerDiameterLimit / 2;
    const statorYokeThickness = outerRadius * 0.15; // approximate rule of thumb
    const slotDepth = outerRadius * 0.25;

    const statorInnerRadius = outerRadius - statorYokeThickness - slotDepth;
    const rotorOuterRadius = statorInnerRadius - specs.airGap;
    const magnetThickness = rotorOuterRadius * 0.1;
    const rotorInnerRadius = rotorOuterRadius * 0.3; // Shaft

    const geometry: MachineGeometry = {
        statorOuterRadius: outerRadius,
        statorInnerRadius,
        rotorOuterRadius,
        rotorInnerRadius,
        slotDepth,
        toothWidth: (2 * Math.PI * statorInnerRadius) / specs.slotCount * 0.5, // 50% tooth/slot ratio
        magnetThickness,
        airGap: specs.airGap,
        slots: specs.slotCount,
        poles: specs.poleCount,
    };

    // 2. Physics Estimation
    const angularVelocity = (specs.ratedSpeed * 2 * Math.PI) / 60; // rad/s

    // Torque T = P / w
    // Rated Power is in kW
    const ratedTorque = (specs.ratedPower * 1000) / (angularVelocity || 1);

    // Estimate losses
    const copperLoss = 0.05 * specs.ratedPower * 1000; // 5% estimate
    const ironLoss = 0.02 * specs.ratedPower * 1000 * (specs.ratedSpeed / 3000); // Speed dependent

    const totalLoss = copperLoss + ironLoss;
    const efficiency = (specs.ratedPower * 1000) / ((specs.ratedPower * 1000) + totalLoss) * 100;

    // Bearingless specific: Suspension force capacity roughly proportional to airgap flux and surface area
    const suspensionForce = (ratedTorque / rotorOuterRadius) * 0.5; // Heuristic for bearingless capacity

    const performance: PerformanceMetrics = {
        efficiency: Math.min(efficiency, 99.9),
        torque: ratedTorque,
        suspensionForce,
        copperLoss,
        ironLoss,
        powerFactor: 0.85 + (specs.ratedSpeed / 100000), // Pseudo-calc
        torqueRipple: 5 + (20 / specs.slotCount), // Less slots, more ripple
        materialCost: outerRadius * specs.axialLengthLimit * 0.05 // Heuristic cost
    };

    return {
        id: crypto.randomUUID(),
        timestamp: Date.now(),
        specs,
        geometry,
        performance
    };
};

export const generateEfficiencyCurve = (specs: DesignSpecs): { speed: number, efficiency: number }[] => {
    const data = [];
    for (let s = 0; s <= specs.ratedSpeed * 1.2; s += specs.ratedSpeed / 10) {
        // Simple efficiency curve shape
        const normSpeed = s / specs.ratedSpeed;
        let eff = 0;
        if (s > 0) {
            eff = 95 * (1 - Math.exp(-3 * normSpeed)) * (1 - 0.1 * Math.pow(normSpeed - 1, 2));
        }
        data.push({
            speed: Math.round(s),
            efficiency: Math.max(0, Math.min(99, eff))
        });
    }
    return data;
};
