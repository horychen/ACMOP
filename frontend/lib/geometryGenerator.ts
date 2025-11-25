import { GP, DesignData, GeometricComponentsObjects, GeometricComponent } from './DesignData';

/**
 * Generate geometric components from GP parameters
 */
export function generateGeometryFromGP(data: DesignData): GeometricComponentsObjects | null {
    if (!data.GP) return null;

    const gp = data.GP;
    const Qs = data.Qs || 12;
    const p = data.p || 4;
    const ps = data.ps || p * 2;

    // Extract parameters with defaults
    const r_so = gp.mm_r_so?.value ?? 125;
    const r_si = gp.mm_r_si?.value ?? (r_so * 0.5);
    const r_ro = gp.mm_r_ro?.value ?? (r_si * 0.95);
    const r_ri = gp.mm_r_ri?.value ?? (r_ro * 0.4);
    const d_pm = gp.mm_d_pm?.value ?? 4;
    const d_mech_air_gap = gp.mm_d_mech_air_gap?.value ?? 0.5;
    const d_sleeve = gp.mm_d_sleeve?.value ?? 1.0;
    const d_st = gp.mm_d_st?.value ?? ((r_si - r_ro - d_mech_air_gap) * 0.6);
    const d_sts = gp.mm_d_sts?.value ?? 3;
    const d_sy = gp.mm_d_sy?.value ?? ((r_so - r_si - d_st - d_sts) * 0.5);
    const w_st = gp.mm_w_st?.value ?? ((2 * Math.PI * r_si) / Qs * 0.4);
    const deg_alpha_st = gp.deg_alpha_st?.value ?? (360 / Qs * 0.7);
    const deg_alpha_rm = gp.deg_alpha_rm?.value ?? (360 / ps * 0.85);
    const deg_alpha_rs = gp.deg_alpha_rs?.value ?? deg_alpha_rm;
    const d_ri = gp.mm_d_ri?.value ?? (r_ro - r_ri - d_pm);
    const d_rp = gp.mm_d_rp?.value ?? 3;

    // Convert angles to radians
    const alpha_st = (deg_alpha_st * Math.PI) / 180;
    const alpha_rm = (deg_alpha_rm * Math.PI) / 180;
    const alpha_rs = (deg_alpha_rs * Math.PI) / 180;
    const slotAngle = (2 * Math.PI) / Qs;
    const poleAngle = (2 * Math.PI) / ps;

    // Generate stator core points
    const statorCore = generateStatorCore(r_so, r_si, r_si + d_st + d_sts, d_sy, Qs, alpha_st, slotAngle);
    
    // Generate rotor core points
    const rotorCore = generateRotorCore(r_ro, r_ri, d_ri, ps, alpha_rm, poleAngle);
    
    // Generate magnets
    const rotorMagnet = generateMagnets(r_ro, d_pm, ps, alpha_rs, poleAngle);
    
    // Generate shaft
    const shaft = generateShaft(r_ri);
    
    // Generate sleeve
    const sleeve = generateSleeve(r_ro + d_mech_air_gap, d_sleeve);
    
    // Generate coils (simplified representation)
    const coils = generateCoils(r_si, d_st, Qs, slotAngle);

    return {
        statorCore,
        rotorCore,
        rotorMagnet,
        shaft,
        sleeve,
        coils
    };
}

function generateStatorCore(
    r_so: number, r_si: number, r_st: number, d_sy: number,
    Qs: number, alpha_st: number, slotAngle: number
): GeometricComponent {
    const points: [number, number][] = [];
    const numPoints = 64; // More points for smoother circle
    
    // Outer circle (stator outer radius)
    for (let i = 0; i <= numPoints; i++) {
        const angle = (i * 2 * Math.PI) / numPoints;
        points.push([r_so * Math.cos(angle), r_so * Math.sin(angle)]);
    }
    
    // Inner circle (stator inner radius) - reverse order
    for (let i = numPoints; i >= 0; i--) {
        const angle = (i * 2 * Math.PI) / numPoints;
        points.push([r_st * Math.cos(angle), r_st * Math.sin(angle)]);
    }
    
    const region = points.map((p, i) => {
        if (i === 0) {
            return { move_to: p };
        } else {
            return { line_to: p };
        }
    });
    
    return {
        name: 'StatorCore',
        color: '#BAFD01',
        list_region: [region],
        ...Object.fromEntries(points.map((p, i) => [`P${i + 1}`, p]))
    };
}

function generateRotorCore(
    r_ro: number, r_ri: number, d_ri: number,
    ps: number, alpha_rm: number, poleAngle: number
): GeometricComponent {
    const points: [number, number][] = [];
    const r_iron = r_ro - d_ri;
    const numPoints = 64;
    
    // Outer circle (rotor outer radius, but excluding magnet area)
    for (let i = 0; i <= numPoints; i++) {
        const angle = (i * 2 * Math.PI) / numPoints;
        points.push([r_iron * Math.cos(angle), r_iron * Math.sin(angle)]);
    }
    
    // Inner circle (rotor inner radius) - reverse order
    for (let i = numPoints; i >= 0; i--) {
        const angle = (i * 2 * Math.PI) / numPoints;
        points.push([r_ri * Math.cos(angle), r_ri * Math.sin(angle)]);
    }
    
    const region = points.map((p, i) => {
        if (i === 0) {
            return { move_to: p };
        } else {
            return { line_to: p };
        }
    });
    
    return {
        name: 'RotorCore',
        color: '#CBD5E1',
        list_region: [region],
        ...Object.fromEntries(points.map((p, i) => [`P${i + 1}`, p]))
    };
}

function generateMagnets(
    r_ro: number, d_pm: number, ps: number, alpha_rs: number, poleAngle: number
): GeometricComponent {
    const regions: any[][] = [];
    const r_inner = r_ro - d_pm;
    const numArcPoints = 20;
    
    for (let i = 0; i < ps; i++) {
        const startAngle = i * poleAngle;
        const endAngle = startAngle + alpha_rs;
        const points: [number, number][] = [];
        
        // Outer arc
        for (let j = 0; j <= numArcPoints; j++) {
            const a = startAngle + (j / numArcPoints) * alpha_rs;
            points.push([r_ro * Math.cos(a), r_ro * Math.sin(a)]);
        }
        // Inner arc (reverse)
        for (let j = numArcPoints; j >= 0; j--) {
            const a = startAngle + (j / numArcPoints) * alpha_rs;
            points.push([r_inner * Math.cos(a), r_inner * Math.sin(a)]);
        }
        
        const region = points.map((p, idx) => {
            if (idx === 0) {
                return { move_to: p };
            } else {
                return { line_to: p };
            }
        });
        
        regions.push(region);
    }
    
    return {
        name: 'RotorMagnet',
        color: '#EF4444',
        list_region: regions,
        ...Object.fromEntries(regions.flat().map((r, i) => [`P${i + 1}`, r]))
    };
}

function generateShaft(r_ri: number): GeometricComponent {
    const numPoints = 32;
    const points: [number, number][] = [];
    
    for (let i = 0; i <= numPoints; i++) {
        const angle = (i * 2 * Math.PI) / numPoints;
        points.push([r_ri * Math.cos(angle), r_ri * Math.sin(angle)]);
    }
    
    return {
        name: 'Shaft',
        color: '#94A3B8',
        list_region: [points.map((p, i) => ({
            move_to: i === 0 ? p : undefined,
            line_to: i > 0 ? p : undefined
        }))],
        ...Object.fromEntries(points.map((p, i) => [`P${i + 1}`, p]))
    };
}

function generateSleeve(r_outer: number, d_sleeve: number): GeometricComponent {
    const r_inner = r_outer - d_sleeve;
    const numPoints = 32;
    const points: [number, number][] = [];
    
    // Outer circle
    for (let i = 0; i <= numPoints; i++) {
        const angle = (i * 2 * Math.PI) / numPoints;
        points.push([r_outer * Math.cos(angle), r_outer * Math.sin(angle)]);
    }
    // Inner circle (reverse)
    for (let i = numPoints; i >= 0; i--) {
        const angle = (i * 2 * Math.PI) / numPoints;
        points.push([r_inner * Math.cos(angle), r_inner * Math.sin(angle)]);
    }
    
    return {
        name: 'Sleeve',
        color: '#64748B',
        list_region: [points.map((p, i) => ({
            move_to: i === 0 ? p : undefined,
            line_to: i > 0 ? p : undefined
        }))],
        ...Object.fromEntries(points.map((p, i) => [`P${i + 1}`, p]))
    };
}

function generateCoils(r_si: number, d_st: number, Qs: number, slotAngle: number): GeometricComponent {
    const regions: any[][] = [];
    const coilRadius = d_st * 0.3;
    
    for (let i = 0; i < Qs; i++) {
        const angle = i * slotAngle + slotAngle / 2;
        const r = r_si + d_st * 0.5;
        const centerX = r * Math.cos(angle);
        const centerY = r * Math.sin(angle);
        
        const numPoints = 16;
        const points: [number, number][] = [];
        for (let j = 0; j <= numPoints; j++) {
            const a = (j * 2 * Math.PI) / numPoints;
            points.push([centerX + coilRadius * Math.cos(a), centerY + coilRadius * Math.sin(a)]);
        }
        
        regions.push(points.map((p, idx) => ({
            move_to: idx === 0 ? p : undefined,
            line_to: idx > 0 ? p : undefined
        })));
    }
    
    return {
        name: 'Coils',
        color: '#B87333',
        list_region: regions,
        ...Object.fromEntries(regions.flat().map((r, i) => [`P${i + 1}`, r]))
    };
}

