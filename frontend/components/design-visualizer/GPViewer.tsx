"use client";

import React from 'react';
import { GP, GPParameter } from '@/lib/DesignData';
import { Input } from '@/components/ui/input';
import { Label } from '@/components/ui/label';
import { ScrollArea } from '@/components/ui/scroll-area';
import { Badge } from '@/components/ui/badge';

interface GPViewerProps {
  data: GP;
  onParameterChange?: (key: string, value: number) => void;
}

export default function GPViewer({ data, onParameterChange }: GPViewerProps) {
  const handleValueChange = (key: string, value: string) => {
    const numValue = parseFloat(value);
    if (!isNaN(numValue) && onParameterChange) {
      onParameterChange(key, numValue);
    }
  };

  const getTypeColor = (type: string) => {
    switch (type) {
      case 'fixed':
        return 'bg-blue-100 text-blue-800 dark:bg-blue-900 dark:text-blue-200';
      case 'free':
        return 'bg-green-100 text-green-800 dark:bg-green-900 dark:text-green-200';
      case 'derived':
        return 'bg-gray-100 text-gray-800 dark:bg-gray-900 dark:text-gray-200';
      default:
        return 'bg-gray-100 text-gray-800 dark:bg-gray-900 dark:text-gray-200';
    }
  };

  return (
    <ScrollArea className="h-full">
      <div className="space-y-3 pr-4">
        {Object.entries(data).map(([key, param]: [string, GPParameter]) => (
          <div key={key} className="space-y-1.5">
            <div className="flex items-center justify-between">
              <Label htmlFor={key} className="text-xs font-medium">
                {key}
              </Label>
              <Badge variant="outline" className={`text-xs ${getTypeColor(param.type)}`}>
                {param.type}
              </Badge>
            </div>
            <div className="space-y-1">
              <Input
                id={key}
                type="number"
                value={param.value ?? ''}
                onChange={(e) => handleValueChange(key, e.target.value)}
                disabled={param.type === 'derived' || !onParameterChange}
                className="h-8 text-xs"
                step="any"
              />
              {param.bounds && (
                <div className="text-xs text-muted-foreground px-1">
                  [{param.bounds[0].toFixed(4)}, {param.bounds[1].toFixed(4)}]
                </div>
              )}
              {param.description && param.description !== key && (
                <div className="text-xs text-muted-foreground px-1">
                  {param.description}
                </div>
              )}
            </div>
          </div>
        ))}
      </div>
    </ScrollArea>
  );
}

