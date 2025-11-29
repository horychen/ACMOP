"use client"

import Link from "next/link"
import { usePathname } from "next/navigation"
import { cn } from "@/lib/utils"
import { Button } from "@/components/ui/button"
import { ScrollArea } from "@/components/ui/scroll-area"
import { LayoutDashboard, Settings, BarChart3, Box, Sparkles, Moon, Sun, FileSearch, Code2 } from "lucide-react"
import { useTheme } from "@/context/ThemeContext"
import { useState, useEffect } from "react"

interface SidebarProps extends React.HTMLAttributes<HTMLDivElement> { }

export function Sidebar({ className }: SidebarProps) {
  const pathname = usePathname()
  const { theme, toggleTheme } = useTheme()
  const [mounted, setMounted] = useState(false)

  useEffect(() => {
    setMounted(true)
  }, [])

  return (
    <div className={cn("pb-12 w-64 border-r min-h-screen bg-card", className)}>
      <div className="space-y-4 py-4">
        <div className="px-3 py-2">
          <h2 className="mb-2 px-4 text-lg font-semibold tracking-tight">
            ACMOP
          </h2>
          <div className="space-y-1">
            <Button variant={pathname === "/" ? "secondary" : "ghost"} className="w-full justify-start" asChild>
              <Link href="/">
                <LayoutDashboard className="mr-2 h-4 w-4" />
                Dashboard
              </Link>
            </Button>
            <Button variant={pathname === "/design" ? "secondary" : "ghost"} className="w-full justify-start" asChild>
              <Link href="/design">
                <Box className="mr-2 h-4 w-4" />
                Design Viewer
              </Link>
            </Button>
            <Button variant={pathname === "/optimization" ? "secondary" : "ghost"} className="w-full justify-start" asChild>
              <Link href="/optimization">
                <BarChart3 className="mr-2 h-4 w-4" />
                Optimization
              </Link>
            </Button>
            <Button variant={pathname === "/visualizer" ? "secondary" : "ghost"} className="w-full justify-start" asChild>
              <Link href="/visualizer">
                <Sparkles className="mr-2 h-4 w-4" />
                Visualizer
              </Link>
            </Button>
            <Button variant={pathname === "/design-visualizer" ? "secondary" : "ghost"} className="w-full justify-start" asChild>
              <Link href="/design-visualizer">
                <Box className="mr-2 h-4 w-4" />
                Design Visualizer (New)
              </Link>
            </Button>
            <Button variant={pathname === "/design-analyzer" ? "secondary" : "ghost"} className="w-full justify-start" asChild>
              <Link href="/design-analyzer">
                <FileSearch className="mr-2 h-4 w-4" />
                Design Analyzer
              </Link>
            </Button>
            <Button variant={pathname === "/csv-visualizer" ? "secondary" : "ghost"} className="w-full justify-start" asChild>
              <Link href="/csv-visualizer">
                <BarChart3 className="mr-2 h-4 w-4" />
                CSV Visualizer
              </Link>
            </Button>
            <Button variant={pathname === "/developer" ? "secondary" : "ghost"} className="w-full justify-start" asChild>
              <Link href="/developer">
                <Code2 className="mr-2 h-4 w-4" />
                开发者页面
              </Link>
            </Button>
          </div>
        </div>

        {mounted && (
          <div className="px-3 py-2 border-t">
            <Button
              variant="ghost"
              className="w-full justify-start"
              onClick={toggleTheme}
            >
              {theme === "dark" ? <Sun className="mr-2 h-4 w-4" /> : <Moon className="mr-2 h-4 w-4" />}
              {theme === "dark" ? "Light Mode" : "Dark Mode"}
            </Button>
          </div>
        )}
      </div>
    </div>
  )
}
