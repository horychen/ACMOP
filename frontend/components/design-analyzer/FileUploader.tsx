'use client';

import React, { useCallback, useState } from 'react';
import { useDropzone } from 'react-dropzone';
import { Upload, FileJson, AlertCircle } from 'lucide-react';
import { Card, CardContent } from '@/components/ui/card';
import { Alert, AlertDescription, AlertTitle } from '@/components/ui/alert';
import { DesignData, parseGPData } from '@/lib/DesignData';

interface FileUploaderProps {
    onDataLoaded: (data: DesignData, fileName: string) => void;
}

export default function FileUploader({ onDataLoaded }: FileUploaderProps) {
    const [error, setError] = useState<string | null>(null);

    const onDrop = useCallback((acceptedFiles: File[]) => {
        setError(null);
        const file = acceptedFiles[0];
        if (!file) return;

        const reader = new FileReader();
        reader.onload = (event) => {
            try {
                const text = event.target?.result as string;
                let json = JSON.parse(text);

                // Handle nested structure: if machine_type is not at top level,
                // try to extract from the first key that contains machine_type
                if (!json.machine_type) {
                    // Check if JSON has a single top-level key with nested data
                    const keys = Object.keys(json);
                    if (keys.length === 1 && typeof json[keys[0]] === 'object' && json[keys[0]] !== null) {
                        // Extract the nested object
                        json = json[keys[0]];
                    } else {
                        // Try to find any key that contains machine_type
                        for (const key of keys) {
                            if (json[key] && typeof json[key] === 'object' && json[key].machine_type) {
                                json = json[key];
                                break;
                            }
                        }
                    }
                }

                // Basic validation
                if (!json.machine_type) {
                    throw new Error("Invalid design file: Missing machine_type");
                }

                // Parse GP data if present
                if (json.GP) {
                    json.GP = parseGPData(json.GP);
                }

                onDataLoaded(json as DesignData, file.name);
            } catch (err) {
                console.error("Error parsing JSON:", err);
                setError(err instanceof Error ? err.message : "Failed to parse JSON file");
            }
        };
        reader.readAsText(file);
    }, [onDataLoaded]);

    const { getRootProps, getInputProps, isDragActive } = useDropzone({
        onDrop,
        accept: {
            'application/json': ['.json']
        },
        multiple: false
    });

    return (
        <div className="w-full max-w-2xl mx-auto mt-10">
            <Card className={`border-2 border-dashed transition-colors ${isDragActive ? 'border-primary bg-primary/5' : 'border-muted-foreground/25'}`}>
                <CardContent className="p-0">
                    <div
                        {...getRootProps()}
                        className="flex flex-col items-center justify-center h-64 cursor-pointer p-6 text-center"
                    >
                        <input {...getInputProps()} />
                        <div className="p-4 rounded-full bg-muted mb-4">
                            {isDragActive ? (
                                <Upload className="h-10 w-10 text-primary animate-bounce" />
                            ) : (
                                <FileJson className="h-10 w-10 text-muted-foreground" />
                            )}
                        </div>
                        <h3 className="text-lg font-semibold mb-1">
                            {isDragActive ? "Drop the file here" : "Drag & drop your design JSON"}
                        </h3>
                        <p className="text-sm text-muted-foreground max-w-xs">
                            Upload a PMSM design file (e.g., PMSM-Q12p4ps5y1-50kW.json) to analyze geometry and performance.
                        </p>
                    </div>
                </CardContent>
            </Card>

            {error && (
                <Alert variant="destructive" className="mt-4">
                    <AlertCircle className="h-4 w-4" />
                    <AlertTitle>Error</AlertTitle>
                    <AlertDescription>{error}</AlertDescription>
                </Alert>
            )}
        </div>
    );
}
