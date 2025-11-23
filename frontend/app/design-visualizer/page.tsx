
import React from 'react';
import fs from 'fs';
import path from 'path';
import DesignVisualizerClient from '@/components/design-visualizer/DesignVisualizerClient';
import { DesignData } from '@/lib/DesignData';

async function getDesignData(): Promise<DesignData | null> {
    // Construct path to backend/DesignVisualizationPickle.json
    // Assuming the app is running in the root or we can navigate up
    // process.cwd() usually points to the project root (frontend) in Next.js
    const filePath = path.resolve(process.cwd(), '../backend/DesignVisualizationPickle.json');

    try {
        if (!fs.existsSync(filePath)) {
            return null;
        }
        const fileContent = fs.readFileSync(filePath, 'utf-8');
        const data = JSON.parse(fileContent);
        return data as DesignData;
    } catch (error) {
        console.error("Error reading design data:", error);
        return null;
    }
}

export default async function DesignVisualizerPage() {
    const data = await getDesignData();

    if (!data) {
        return (
            <div className="flex flex-col items-center justify-center h-screen space-y-4">
                <h1 className="text-2xl font-bold text-red-600">Design Data Not Found</h1>
                <p className="text-muted-foreground">
                    Could not find <code>backend/DesignVisualizationPickle.json</code>.
                </p>
                <p>
                    Please run <code>acmop.py</code> manually to generate the visualization data.
                </p>
            </div>
        );
    }

    return <DesignVisualizerClient data={data} />;
}
