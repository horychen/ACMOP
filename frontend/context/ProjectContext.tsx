"use client"

import React, { createContext, useContext, useState, useEffect } from "react"
import axios from "axios"

interface ProjectContextType {
    projects: string[]
    currentProject: string
    setCurrentProject: (project: string) => void
    loading: boolean
}

const ProjectContext = createContext<ProjectContextType | undefined>(undefined)

export function ProjectProvider({ children }: { children: React.ReactNode }) {
    const [projects, setProjects] = useState<string[]>([])
    const [currentProject, setCurrentProject] = useState<string>("")
    const [loading, setLoading] = useState(true)

    useEffect(() => {
        axios.get("http://localhost:8000/api/projects/")
            .then((res) => {
                setProjects(res.data)
                if (res.data.length > 0) {
                    // Check if there's a stored project in localStorage
                    const stored = localStorage.getItem("currentProject")
                    if (stored && res.data.includes(stored)) {
                        setCurrentProject(stored)
                    } else {
                        setCurrentProject(res.data[0])
                    }
                }
            })
            .catch((err) => console.error("Failed to fetch projects", err))
            .finally(() => setLoading(false))
    }, [])

    const handleSetProject = (project: string) => {
        setCurrentProject(project)
        localStorage.setItem("currentProject", project)
    }

    return (
        <ProjectContext.Provider value={{ projects, currentProject, setCurrentProject: handleSetProject, loading }}>
            {children}
        </ProjectContext.Provider>
    )
}

export function useProject() {
    const context = useContext(ProjectContext)
    if (context === undefined) {
        throw new Error("useProject must be used within a ProjectProvider")
    }
    return context
}
