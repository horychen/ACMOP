"use client"

import { ProjectProvider } from "@/context/ProjectContext"
import { ThemeProvider } from "@/context/ThemeContext"

export function ClientProviders({ children }: { children: React.ReactNode }) {
    return (
        <ThemeProvider>
            <ProjectProvider>
                {children}
            </ProjectProvider>
        </ThemeProvider>
    )
}
