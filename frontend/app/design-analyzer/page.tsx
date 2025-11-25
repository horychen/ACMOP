'use client';

import React, { useState } from 'react';
import { DesignData } from '@/lib/DesignData';
import FileUploader from '@/components/design-analyzer/FileUploader';
import AnalyzerDashboard from '@/components/design-analyzer/AnalyzerDashboard';
import { Button } from '@/components/ui/button';
import { ArrowLeft } from 'lucide-react';

export default function DesignAnalyzerPage() {
    const [designData, setDesignData] = useState<DesignData | null>(null);
    const [fileName, setFileName] = useState<string>("");

    const handleDataLoaded = (data: DesignData, name: string) => {
        setDesignData(data);
        setFileName(name);
    };

    const handleReset = () => {
        setDesignData(null);
        setFileName("");
    };

    return (
        <div className="container mx-auto p-4 h-screen flex flex-col">
            <div className="flex items-center justify-between mb-6">
                <div className="flex items-center gap-4">
                    {designData && (
                        <Button variant="ghost" size="icon" onClick={handleReset}>
                            <ArrowLeft className="h-4 w-4" />
                        </Button>
                    )}
                    <h1 className="text-3xl font-bold tracking-tight">Design Analyzer</h1>
                </div>
            </div>

            <div className="flex-1 overflow-hidden">
                {!designData ? (
                    <FileUploader onDataLoaded={handleDataLoaded} />
                ) : (
                    <AnalyzerDashboard initialData={designData} fileName={fileName} />
                )}
            </div>
        </div>
    );
}
