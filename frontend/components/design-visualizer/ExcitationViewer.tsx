'use client';

import React from 'react';
import { ExUser } from '@/lib/DesignData';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import { Activity, Zap, Gauge } from 'lucide-react';

interface ExcitationViewerProps {
    data: ExUser;
}

export default function ExcitationViewer({ data }: ExcitationViewerProps) {
    const renderCard = (title: string, value: string | number, unit: string, icon: React.ReactNode) => (
        <Card>
            <CardHeader className="flex flex-row items-center justify-between space-y-0 pb-2">
                <CardTitle className="text-sm font-medium">
                    {title}
                </CardTitle>
                {icon}
            </CardHeader>
            <CardContent>
                <div className="text-2xl font-bold">{value} <span className="text-sm font-normal text-muted-foreground">{unit}</span></div>
            </CardContent>
        </Card>
    );

    return (
        <div className="grid gap-4 md:grid-cols-2 lg:grid-cols-4">
            {renderCard("Drive Frequency", data.DriveW_Freq, "Hz", <Activity className="h-4 w-4 text-muted-foreground" />)}
            {renderCard("Drive Current", data.DriveW_CurrentAmp?.toFixed(2), "A", <Zap className="h-4 w-4 text-muted-foreground" />)}
            {renderCard("Bearing Frequency", data.BeariW_Freq, "Hz", <Activity className="h-4 w-4 text-muted-foreground" />)}
            {renderCard("Bearing Current", data.BeariW_CurrentAmp?.toFixed(2), "A", <Zap className="h-4 w-4 text-muted-foreground" />)}
            {renderCard("Voltage Rating", data.VoltageRating, "V", <Zap className="h-4 w-4 text-muted-foreground" />)}
            {renderCard("Speed", data.the_speed, "RPM", <Gauge className="h-4 w-4 text-muted-foreground" />)}
            {renderCard("Current Density", (data.Js / 1e6).toFixed(2), "A/mm²", <Zap className="h-4 w-4 text-muted-foreground" />)}
            {renderCard("Fill Factor", (data.WindingFill * 100).toFixed(1), "%", <Activity className="h-4 w-4 text-muted-foreground" />)}
        </div>
    );
}
