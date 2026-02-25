import React, { useState, useMemo } from 'react';

const MotorWindingSimulator = () => {
    const [activeGauge, setActiveGauge] = useState("AWG 31");

    // 核心几何参数 (单位: mm)
    const motor = {
        od: 13.0,
        rotorOd: 8.0,
        airGap: 0.15,
        statorId: 8.3,
        toothDepth: 2.1,
        toothWidth: 1.2,
        yoke: 0.25,
        slots: 12,
        liner: 0.1
    };

    const wires = {
        "AWG 30": { d_bare: 0.254, d_od: 0.290, color: '#f87171' },
        "AWG 31": { d_bare: 0.226, d_od: 0.260, color: '#60a5fa' },
        "AWG 32": { d_bare: 0.203, d_od: 0.235, color: '#4ade80' }
    };

    const currentGauge = wires[activeGauge];

    // 1. 槽面积几何计算 (mm^2)
    const slotArea = useMemo(() => {
        const rIn = motor.statorId / 2;
        const rOut = rIn + motor.toothDepth;
        const totalAnnulus = Math.PI * (Math.pow(rOut, 2) - Math.pow(rIn, 2));
        const totalTeethArea = motor.slots * motor.toothWidth * motor.toothDepth;
        return (totalAnnulus - totalTeethArea) / motor.slots;
    }, []);

    // 2. 坐标变换函数 (以 12 点钟方向为 0 度)
    const toCartesian = (r, theta, S) => ({
        x: r * Math.sin(theta) * S,
        y: -r * Math.cos(theta) * S
    });

    // 3. 绕线排布算法与电机学分析
    const windingData = useMemo(() => {
        const w = currentGauge;
        const rIn = motor.statorId / 2;
        const slotAngle = (2 * Math.PI) / motor.slots;
        const results = [];
        // 两个相邻齿，齿1在左，齿2在右，中间形成槽
        const toothAngles = [-slotAngle / 2, slotAngle / 2];

        toothAngles.forEach((tAngle, idx) => {
            const side = idx === 0 ? 1 : -1;
            const coilColor = idx === 0 ? '#3b82f6' : '#ef4444';

            for (let layer = 0; layer < 6; layer++) {
                // 计算层位横向偏移
                const hOffset = (motor.toothWidth / 2 + motor.liner + w.d_od / 2 + layer * w.d_od * 0.88);
                for (let row = 0; row < 15; row++) {
                    const vOffset = motor.liner + w.d_od / 2 + row * w.d_od;
                    if (vOffset > motor.toothDepth - motor.liner) break;
                    const r = rIn + vOffset;
                    if (hOffset >= r) continue;

                    const thetaOffset = Math.asin(hOffset / r) * side;
                    const finalAngle = tAngle + thetaOffset;

                    // 槽中心碰撞检测 (允许极小重叠或预留间隙)
                    if (side === 1 && finalAngle > -0.002) continue;
                    if (side === -1 && finalAngle < 0.002) continue;

                    results.push({
                        ...toCartesian(r, finalAngle, 1),
                        r: w.d_od / 2,
                        color: coilColor
                    });
                }
            }
        });

        const totalWires = results.length;
        const bareCuArea = totalWires * (Math.PI * Math.pow(w.d_bare / 2, 2));
        const targetJ = 14; // A/mm^2
        const iRated = (Math.PI * Math.pow(w.d_bare / 2, 2)) * targetJ;

        // 线负荷 A (安培/厘米)
        const lineLoadA_cm = ((totalWires * motor.slots * iRated) / (Math.PI * motor.statorId)) * 10;
        const aj = lineLoadA_cm * targetJ;

        return {
            wires: results,
            totalWires,
            cuFill: (bareCuArea / slotArea) * 100,
            totalAmps: totalWires * iRated, // 单槽总安匝
            lineLoadA: lineLoadA_cm,
            aj: aj,
            iRated
        };
    }, [activeGauge, slotArea]);

    // 绘图比例: 1mm = 200px (确保在电脑端能看清微米级间隙)
    const S = 200;

    return (
        <div className="flex flex-col lg:flex-row w-full h-screen bg-slate-950 text-slate-200 overflow-hidden font-sans">

            {/* 左侧面板: 核心控制与电机数据 */}
            <div className="w-full lg:w-[380px] p-6 bg-slate-900/50 border-b lg:border-b-0 lg:border-r border-slate-800 flex flex-col overflow-y-auto">
                <div className="mb-6">
                    <h2 className="text-2xl font-black bg-gradient-to-r from-blue-400 to-cyan-400 bg-clip-text text-transparent">
                        电机下线设计仿真
                    </h2>
                    <p className="text-slate-500 text-[10px] font-mono mt-1">12S10P MICRO MOTOR DESIGN SYSTEM</p>
                </div>

                <div className="flex gap-2 mb-6">
                    {Object.keys(wires).map(g => (
                        <button
                            key={g}
                            onClick={() => setActiveGauge(g)}
                            className={`flex-1 py-3 rounded-xl border-2 text-xs font-bold transition-all ${activeGauge === g ? 'bg-blue-600 border-blue-400 text-white shadow-xl scale-105' : 'bg-slate-800 border-slate-700 opacity-50 hover:opacity-100'
                                }`}
                        >
                            {g}
                        </button>
                    ))}
                </div>

                <div className="space-y-4">
                    <StatBox label="单槽总导线数 (Z_slot)" value={windingData.totalWires} unit="根" sub="双层绕组总和" />
                    <StatBox label="安匝激励 (NI)" value={windingData.totalAmps.toFixed(1)} unit="A·t" color="text-cyan-400" sub="单槽总磁动势" />
                    <StatBox label="铜槽满率 (η_cu)" value={windingData.cuFill.toFixed(1)} unit="%" color="text-emerald-400" />

                    <div className="h-[1px] bg-slate-800 my-4" />

                    <StatBox label="线负荷 (A)" value={windingData.lineLoadA.toFixed(1)} unit="A/cm" color="text-yellow-400" />
                    <StatBox label="热负荷 (AJ)" value={windingData.aj.toFixed(0)} unit="A²/(cm·mm²)" color="text-red-500" />
                </div>

                <div className="mt-auto pt-6 border-t border-slate-800">
                    <div className="p-4 bg-red-950/20 border border-red-900/40 rounded-2xl text-[10px] text-red-300 leading-relaxed">
                        <strong>工程警告:</strong> AJ 值超过 1000 表明产热极高。目前 0.25mm 的轭部将导致磁饱和，并阻碍齿部热量向外壳传导。建议校核瞬态温升。
                    </div>
                </div>
            </div>

            {/* 右侧主绘图区: 自适应 SVG */}
            <div className="flex-1 bg-black relative flex flex-col">
                <div className="absolute top-6 left-8 z-20 flex gap-6 text-[10px] uppercase font-bold text-slate-500">
                    <div className="flex items-center gap-2"><div className="w-3 h-3 rounded-full bg-[#3b82f6]"></div> 齿1绕组边</div>
                    <div className="flex items-center gap-2"><div className="w-3 h-3 rounded-full bg-[#ef4444]"></div> 齿2绕组边</div>
                </div>

                <div className="flex-1 w-full h-full flex items-center justify-center p-4">
                    {/* viewBox 优化: 宽度 800, 高度 700. 覆盖 y轴从 -1300 到 -600 的槽部区域 */}
                    <svg
                        viewBox="-400 -1350 800 700"
                        className="w-full h-full max-h-[90vh] transition-all"
                        preserveAspectRatio="xMidYMid meet"
                    >
                        {/* 1. 定子轭部 (Yoke) */}
                        <path
                            d={`M -300 -6.5 A 6.5 6.5 0 0 1 300 -6.5`}
                            fill="none" stroke="#1e293b" strokeWidth={motor.yoke * S}
                            transform={`scale(${S / S})`}
                        />
                        <text x="0" y="-1320" textAnchor="middle" fill="#ef4444" fontSize="14" fontWeight="black">YOKE 0.25mm</text>

                        {/* 2. 气隙与转子边界 (显示下方) */}
                        <line x1="-350" y1="-415" x2="350" y2="-415" stroke="#334155" strokeWidth="1" strokeDasharray="5,5" />
                        <text x="0" y="-780" textAnchor="middle" fill="#475569" fontSize="12" fontWeight="bold">AIR GAP 0.15mm</text>

                        {/* 3. 相邻齿绘制 (形成物理槽位) */}
                        {[-1, 1].map(side => {
                            const angle = side * (Math.PI / 12);
                            const rIn = motor.statorId / 2;
                            const rOut = rIn + motor.toothDepth;
                            const hW = motor.toothWidth / 2;
                            const p1 = toCartesian(rIn, angle + Math.asin(hW / rIn), S);
                            const p2 = toCartesian(rIn, angle - Math.asin(hW / rIn), S);
                            const p3 = toCartesian(rOut, angle - Math.asin(hW / rOut), S);
                            const p4 = toCartesian(rOut, angle + Math.asin(hW / rOut), S);
                            return (
                                <path
                                    key={side}
                                    d={`M ${p1.x} ${p1.y} L ${p2.x} ${p2.y} L ${p3.x} ${p3.y} L ${p4.x} ${p4.y} Z`}
                                    fill="#1e293b" stroke="#475569" strokeWidth="1.5"
                                />
                            );
                        })}

                        {/* 4. 槽绝缘线 (Liner) */}
                        <path
                            d={`M ${-0.038 * S} ${-4.3 * S} L ${-0.038 * S} ${-6.2 * S} L ${0.038 * S} ${-6.2 * S} L ${0.038 * S} ${-4.3 * S}`}
                            fill="none" stroke="#ca8a04" strokeWidth="2" strokeDasharray="8,4" opacity="0.4"
                        />

                        {/* 5. 渲染槽内所有导线 */}
                        {windingData.wires.map((w, i) => (
                            <circle
                                key={i} cx={w.x * S} cy={w.y * S} r={w.r * S}
                                fill={w.color} stroke="black" strokeWidth="0.5"
                            />
                        ))}

                        <text x="0" y="-1150" textAnchor="middle" fill="#ca8a04" fontSize="30" fontWeight="black" opacity="0.1">SLOT AREA</text>
                    </svg>
                </div>

                {/* 底部信息栏 */}
                <div className="p-6 bg-slate-900/80 backdrop-blur-lg border-t border-slate-800 flex justify-between items-center">
                    <div className="text-[10px] text-slate-500 font-mono">
                        [GEOMETRY] SLOT WIDTH: 0.97-2.07mm | DEPTH: 2.1mm
                    </div>
                    <div className="text-right">
                        <p className="text-[10px] text-slate-500 font-bold uppercase">仿真比例</p>
                        <p className="text-xs text-blue-400 font-bold">1mm : 200px</p>
                    </div>
                </div>
            </div>
        </div>
    );
};

const StatBox = ({ label, value, unit, color = "text-white", sub }) => (
    <div className="bg-slate-950/50 p-4 rounded-2xl border border-slate-800 hover:border-slate-700 transition-colors">
        <p className="text-[10px] text-slate-500 uppercase font-black tracking-widest">{label}</p>
        <div className="flex items-baseline gap-1 my-1">
            <span className={`text-2xl font-mono font-black ${color}`}>{value}</span>
            <span className="text-[10px] text-slate-600 font-bold uppercase">{unit}</span>
        </div>
        {sub && <p className="text-[9px] text-slate-600 italic">{sub}</p>}
    </div>
);

export default MotorWindingSimulator;