"use client";

import React from 'react';
import MultiCsvVisualizer from '../../components/MultiCsvVisualizer';

export default function CsvVisualizerPage() {
    // Hardcoded project name for now, could be dynamic later
    const projectName = "SuperCoolPMSM";

    return (
        <div className="min-h-screen flex flex-col font-sans bg-background text-foreground">
            <header className="bg-card border-b border-border p-4 flex items-center justify-between sticky top-0 z-50">
                <div className="flex items-center space-x-3">
                    <h1 className="text-xl font-bold tracking-tight">CSV <span className="text-primary font-mono text-sm">Visualizer</span></h1>
                </div>
            </header>

            <main className="flex-1 p-6 overflow-y-auto">
                <MultiCsvVisualizer projectName={projectName} />
            </main>
        </div>
    );
}
