"use client"

import * as React from "react"
import {
    Select,
    SelectContent,
    SelectGroup,
    SelectItem,
    SelectLabel,
    SelectTrigger,
    SelectValue,
} from "@/components/ui/select"
import { useProject } from "@/context/ProjectContext"

export function ProjectSelector({ onSelect }: { onSelect?: (project: string) => void }) {
    const { projects, currentProject, setCurrentProject } = useProject()

    const handleSelect = (value: string) => {
        setCurrentProject(value)
        if (onSelect) onSelect(value)
    }

    return (
        <Select value={currentProject} onValueChange={handleSelect}>
            <SelectTrigger className="w-[280px]">
                <SelectValue placeholder="Select a project" />
            </SelectTrigger>
            <SelectContent>
                <SelectGroup>
                    <SelectLabel>Projects</SelectLabel>
                    {projects.map((p) => (
                        <SelectItem key={p} value={p}>
                            {p}
                        </SelectItem>
                    ))}
                </SelectGroup>
            </SelectContent>
        </Select>
    )
}
