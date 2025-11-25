"use client"

import { useState } from "react";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import DeveloperMode from "@/components/modes/DeveloperMode";
import VisualizationMode from "@/components/modes/VisualizationMode";

export default function Home() {
  const [activeMode, setActiveMode] = useState<"developer" | "visualization">("developer");

  return (
    <div className="container mx-auto py-8 space-y-6">
      <div className="flex items-center justify-between">
        <div>
          <h1 className="text-4xl font-bold tracking-tight">ACMOP</h1>
          <p className="text-muted-foreground mt-2">
            AC Machine Optimization Platform
          </p>
        </div>
      </div>

      <Tabs value={activeMode} onValueChange={(v) => setActiveMode(v as "developer" | "visualization")} className="w-full">
        <TabsList className="grid w-full max-w-md grid-cols-2">
          <TabsTrigger value="developer">开发者模式</TabsTrigger>
          <TabsTrigger value="visualization">可视化模式</TabsTrigger>
        </TabsList>
        
        <TabsContent value="developer" className="mt-6">
          <DeveloperMode />
        </TabsContent>
        
        <TabsContent value="visualization" className="mt-6">
          <VisualizationMode />
        </TabsContent>
      </Tabs>
    </div>
  );
}

export const dynamic = 'force-dynamic'
