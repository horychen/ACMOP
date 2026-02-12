"use client";

import React, { useState } from "react";
import { CheckCircle2, Circle, Lock } from "lucide-react";
import { cn } from "@/lib/utils";
import { Button } from "@/components/ui/button";
import MachineVisualizer from "./MachineVisualizer";

const STEPS = [
    { id: "materials", title: "材料选择", description: "Material Selection" },
    { id: "geometry", title: "几何尺寸", description: "Geometry & Space Exploration" },
    { id: "winding", title: "绕组与电磁", description: "Winding & Excitation" },
    { id: "targets", title: "优化与结果", description: "Targets & FEA Waves" },
];

export default function DesignStepper() {
    const [currentStep, setCurrentStep] = useState(0);
    const [confirmedSteps, setConfirmedSteps] = useState<string[]>([]);
    const [isCurrentStepValid, setIsCurrentStepValid] = useState(true);

    const handleConfirm = () => {
        if (!isCurrentStepValid) return;

        const stepId = STEPS[currentStep].id;
        if (!confirmedSteps.includes(stepId)) {
            setConfirmedSteps([...confirmedSteps, stepId]);
        }
        if (currentStep < STEPS.length - 1) {
            setCurrentStep(currentStep + 1);
        }
    };

    const handleReset = () => {
        if (confirm("确定要重置所有设计进度吗？")) {
            setCurrentStep(0);
            setConfirmedSteps([]);
        }
    };

    return (
        <div className="space-y-6">
            {/* Header with Steps and Navigation */}
            <div className="bg-white border rounded-xl p-4 shadow-sm sticky top-0 z-20">
                <div className="flex flex-col md:flex-row items-center gap-6">
                    {/* Step Bar */}
                    <div className="relative flex-1 flex justify-between w-full">
                        <div className="absolute top-5 left-0 w-full h-0.5 bg-slate-100 -z-10" />
                        {STEPS.map((step, index) => {
                            const isConfirmed = confirmedSteps.includes(step.id);
                            const isActive = currentStep === index;
                            const isLocked = index > currentStep && !isConfirmed;

                            return (
                                <div
                                    key={step.id}
                                    className={cn(
                                        "flex flex-col items-center gap-1 bg-white z-10 transition-all",
                                        isActive ? "px-4 scale-105" : "px-2",
                                        isLocked ? "text-slate-300" : "text-slate-600"
                                    )}
                                >
                                    <div
                                        className={cn(
                                            "w-9 h-9 rounded-full flex items-center justify-center border-2 transition-all shadow-sm",
                                            isConfirmed
                                                ? "bg-green-500 border-green-500 text-white"
                                                : isActive
                                                    ? "bg-primary border-primary text-primary-foreground font-bold ring-4 ring-primary/10"
                                                    : "bg-white border-slate-200"
                                        )}
                                    >
                                        {isConfirmed ? (
                                            <CheckCircle2 className="w-5 h-5" />
                                        ) : isLocked ? (
                                            <Lock className="w-4 h-4" />
                                        ) : (
                                            <span className="text-sm">{index + 1}</span>
                                        )}
                                    </div>
                                    <p className={cn("text-[11px] font-bold uppercase tracking-tighter", isActive ? "text-primary" : "text-slate-400")}>
                                        {step.title}
                                    </p>
                                </div>
                            );
                        })}
                    </div>

                    {/* Navigation Buttons in Header */}
                    <div className="flex items-center gap-2 border-l pl-6">
                        <Button
                            variant="ghost"
                            size="sm"
                            onClick={() => currentStep > 0 && setCurrentStep(currentStep - 1)}
                            disabled={currentStep === 0}
                            className="h-9 px-4"
                        >
                            上一步
                        </Button>
                        <Button
                            variant="outline"
                            size="sm"
                            onClick={handleReset}
                            className="h-9 text-slate-500 hover:text-red-500"
                        >
                            重置
                        </Button>
                        <Button
                            onClick={handleConfirm}
                            disabled={!isCurrentStepValid}
                            size="sm"
                            className={cn(
                                "h-9 min-w-[120px] font-bold shadow-md transition-all",
                                isCurrentStepValid ? "bg-green-600 hover:bg-green-700 text-white" : "bg-slate-200 text-slate-400"
                            )}
                        >
                            {currentStep === STEPS.length - 1 ? "完成设计" : `下一步: ${STEPS[currentStep + 1].title}`}
                        </Button>
                    </div>
                </div>
            </div>

            {/* Content Area */}
            <div className="bg-slate-50/30 rounded-xl min-h-[600px] border shadow-inner">
                <MachineVisualizer
                    currentStepId={STEPS[currentStep].id}
                    onValidationChange={setIsCurrentStepValid}
                />
            </div>
        </div>
    );
}
