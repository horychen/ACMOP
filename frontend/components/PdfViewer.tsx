"use client";

import React, { useState, useEffect } from 'react';
import { Loader2, AlertCircle, ZoomIn, ZoomOut, RotateCcw } from 'lucide-react';

interface PdfViewerProps {
    pdfUrl: string;
    className?: string;
}

export default function PdfViewer({ pdfUrl, className = "" }: PdfViewerProps) {
    const [loading, setLoading] = useState(true);
    const [error, setError] = useState<string | null>(null);
    const [scale, setScale] = useState(1.0);
    const iframeRef = React.useRef<HTMLIFrameElement>(null);
    const timeoutRef = React.useRef<NodeJS.Timeout | null>(null);

    useEffect(() => {
        setLoading(true);
        setError(null);
        
        // Set timeout for loading (10 seconds)
        timeoutRef.current = setTimeout(() => {
            setError("PDF 加载超时，请检查文件是否存在");
            setLoading(false);
        }, 10000);

        return () => {
            if (timeoutRef.current) {
                clearTimeout(timeoutRef.current);
            }
        };
    }, [pdfUrl]);

    const handleZoomIn = () => {
        setScale(prev => Math.min(prev + 0.1, 3.0));
    };

    const handleZoomOut = () => {
        setScale(prev => Math.max(prev - 0.1, 0.5));
    };

    const handleReset = () => {
        setScale(1.0);
    };

    return (
        <div className={`flex flex-col h-full ${className}`}>
            {/* Controls */}
            <div className="flex items-center justify-between p-2 bg-card border-b border-border">
                <div className="flex items-center space-x-2">
                    <button
                        onClick={handleZoomOut}
                        className="p-1.5 rounded hover:bg-muted transition-colors"
                        title="缩小"
                    >
                        <ZoomOut className="w-4 h-4" />
                    </button>
                    <span className="text-xs text-muted-foreground min-w-[60px] text-center">
                        {Math.round(scale * 100)}%
                    </span>
                    <button
                        onClick={handleZoomIn}
                        className="p-1.5 rounded hover:bg-muted transition-colors"
                        title="放大"
                    >
                        <ZoomIn className="w-4 h-4" />
                    </button>
                    <button
                        onClick={handleReset}
                        className="p-1.5 rounded hover:bg-muted transition-colors ml-2"
                        title="重置"
                    >
                        <RotateCcw className="w-4 h-4" />
                    </button>
                </div>
            </div>

            {/* PDF Viewer */}
            <div className="flex-1 relative overflow-auto bg-muted/30">
                {error ? (
                    <div className="absolute inset-0 flex items-center justify-center">
                        <div className="text-center p-6">
                            <AlertCircle className="w-8 h-8 text-destructive mx-auto mb-2" />
                            <p className="text-sm text-muted-foreground">{error}</p>
                        </div>
                    </div>
                ) : (
                    <>
                        {loading && (
                            <div className="absolute inset-0 flex items-center justify-center">
                                <Loader2 className="w-6 h-6 animate-spin text-primary" />
                            </div>
                        )}
                        <div 
                            className="w-full h-full flex items-center justify-center p-4 overflow-auto"
                            style={{ transform: `scale(${scale})`, transformOrigin: 'top center' }}
                        >
                            <iframe
                                ref={iframeRef}
                                src={`${pdfUrl}#toolbar=0&navpanes=0&scrollbar=1`}
                                className="w-full h-full border-0"
                                onLoad={() => {
                                    if (timeoutRef.current) {
                                        clearTimeout(timeoutRef.current);
                                    }
                                    setLoading(false);
                                }}
                                onError={() => {
                                    if (timeoutRef.current) {
                                        clearTimeout(timeoutRef.current);
                                    }
                                    setError("无法加载 PDF 文件");
                                    setLoading(false);
                                }}
                                title="Machine Geometry PDF"
                                style={{ minHeight: '600px' }}
                            />
                        </div>
                    </>
                )}
            </div>
        </div>
    );
}

