
export interface Point {
    [key: number]: number; // Array like [x, y]
}

export interface GeometricComponent {
    name: string;
    color: string;
    alpha_notch?: number | null;
    list_region?: any[][];
    [key: string]: any; // For dynamic points like P1, P2, etc.
}

export interface GeometricComponentsObjects {
    rotorCore: GeometricComponent;
    shaft: GeometricComponent;
    rotorMagnet: GeometricComponent;
    sleeve: GeometricComponent;
    statorCore: GeometricComponent;
    coils: GeometricComponent;
}

export interface WindingLayout {
    layer_X_phases: string[];
    layer_X_signs: string[];
    layer_Y_phases: string[];
    layer_Y_signs: string[];
    ox_distribution_three_phase: string[];
    coil_pitch_y?: number;
    kd1?: number;
    kp1?: number;
    [key: string]: any;
}

export interface ExUser {
    ExcitationFreqSimulated: number;
    VoltageRating: number;
    mm_stack_length: number;
    Steel: string;
    Coil: string;
    Js: number;
    Temperature: number;
    WindingFill: number;
    wily: WindingLayout;
    DriveW_CurrentAmp: number;
    DriveW_Freq: number;
    DriveW_zQ?: number;
    BeariW_CurrentAmp: number;
    BeariW_Freq: number;
    BeariW_zQ?: number;
    the_speed: number;
    Omega: number;
    [key: string]: any;
}

export interface GpUser {
    [key: string]: {
        value: number | null;
        type: string;
    } | number | null;
}

export interface GPParameter {
    type: 'fixed' | 'free' | 'derived';
    description: string;
    value: number;
    bounds: [number, number] | null;
}

export interface GP {
    [key: string]: GPParameter;
}

export interface FEA_Performance {
    project_name: string;
    individual_name: string;
    f1: number;
    f2: number;
    f3: number;
    TRV: number;
    FRW: number;
    torque_average: number;
    ss_avg_force_magnitude: number;
    rotor_weight: number;
    normalized_torque_ripple: number;
    normalized_force_error_magnitude: number;
    force_error_angle: number;
    coil_flux_linkage_peak2peak_value: number;
    mm2_slot_area: number;
    Cost: number;
    Cost_Fe: number;
    Cost_Cu: number;
    Cost_PM: number;
    power_factor: number;
    rated_ratio: number;
    rated_stack_length_mm: number;
    rated_total_loss: number;
    rated_stator_copper_loss_along_stack: number;
    rated_rotor_copper_loss_along_stack: number;
    rated_magnet_Joule_loss: number;
    stator_copper_loss_in_end_turn: number;
    rotor_copper_loss_in_end_turn: number;
    rated_iron_loss: number;
    rated_windage_loss: number;
    [key: string]: any;
}

export interface DesignData {
    machine_type: string;
    m: number;
    Qs: number;
    p: number;
    ps: number;
    mec_power: number;
    guess_efficiency: number;
    GeometricComponentsObjects: GeometricComponentsObjects;
    "EX-user": ExUser;
    "GP-user": GpUser;
    GP?: GP;
    "FEA_Evaluated_Performance--1-Initial"?: FEA_Performance;
    [key: string]: any;
}

export function parseDesignData(json: any): DesignData {
    return json as DesignData;
}

export function getPoints(component: GeometricComponent): [number, number][] {
    if (!component) return [];
    const points: [number, number][] = [];
    const keys = Object.keys(component).filter(k => /^P\d+/.test(k) || /^P[A-Za-z]+/.test(k));

    const pointKeys = keys.sort((a, b) => {
        const numA = parseFloat(a.replace('P', '').replace('_', '.'));
        const numB = parseFloat(b.replace('P', '').replace('_', '.'));
        if (!isNaN(numA) && !isNaN(numB)) return numA - numB;
        return a.localeCompare(b);
    });

    pointKeys.forEach(key => {
        const val = component[key];
        if (Array.isArray(val) && val.length === 2 && val[0] !== null && val[1] !== null) {
            points.push(val as [number, number]);
        }
    });

    return points;
}

/**
 * Parse GP (Geometric Parameters) from the JSON pickle structure
 * Converts the complex py/object structure to a simple GP interface
 */
export function parseGPData(gpRaw: any): GP | undefined {
    if (!gpRaw || !gpRaw['py/reduce']) return undefined;

    const tuples = gpRaw['py/reduce']?.[4]?.['py/tuple'];
    if (!tuples || !Array.isArray(tuples)) return undefined;

    const gp: GP = {};

    tuples.forEach((item: any) => {
        if (!item['py/tuple'] || item['py/tuple'].length < 2) return;

        const key = item['py/tuple'][0];
        const paramObj = item['py/tuple'][1];

        if (!paramObj?.['py/state']?.['py/tuple']) return;

        const state = paramObj['py/state']['py/tuple'];
        const type = state[0] as 'fixed' | 'free' | 'derived';
        const description = state[1];
        const value = state[2];
        const boundsRaw = state[3];

        let bounds: [number, number] | null = null;
        if (Array.isArray(boundsRaw) && boundsRaw.length === 2 &&
            boundsRaw[0] !== null && boundsRaw[1] !== null) {
            bounds = [boundsRaw[0], boundsRaw[1]];
        }

        gp[key] = {
            type,
            description,
            value,
            bounds
        };
    });

    return gp;
}

