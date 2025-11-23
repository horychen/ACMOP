
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
