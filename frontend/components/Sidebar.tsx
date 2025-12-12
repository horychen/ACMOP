"use client"

import Link from "next/link"
import { usePathname } from "next/navigation"
import { cn } from "@/lib/utils"
import { Button } from "@/components/ui/button"
import { ScrollArea } from "@/components/ui/scroll-area"
import { LayoutDashboard, Settings, BarChart3, Box, Sparkles, Moon, Sun, FileSearch, Code2, ChevronLeft, ChevronRight } from "lucide-react"
import { useTheme } from "@/context/ThemeContext"
import { useState, useEffect } from "react"

interface SidebarProps extends React.HTMLAttributes<HTMLDivElement> { }

const SIDEBAR_COLLAPSED_KEY = "sidebar-collapsed"

export function Sidebar({ className }: SidebarProps) {
  const pathname = usePathname()
  const { theme, toggleTheme } = useTheme()
  const [mounted, setMounted] = useState(false)
  const [isCollapsed, setIsCollapsed] = useState(false)

  useEffect(() => {
    setMounted(true)
    // 从 localStorage 读取折叠状态
    const savedState = localStorage.getItem(SIDEBAR_COLLAPSED_KEY)
    if (savedState !== null) {
      setIsCollapsed(savedState === "true")
    }
  }, [])

  const toggleCollapse = () => {
    const newState = !isCollapsed
    setIsCollapsed(newState)
    localStorage.setItem(SIDEBAR_COLLAPSED_KEY, String(newState))
  }

  return (
    <div className={cn(
      "pb-12 border-r min-h-screen bg-card transition-all duration-300 relative",
      isCollapsed ? "w-16" : "w-64",
      className
    )}>
      {/* 折叠/展开按钮 */}
      <Button
        variant="ghost"
        size="icon"
        className="absolute -right-3 top-4 z-10 h-6 w-6 rounded-full border bg-background shadow-sm hover:bg-accent"
        onClick={toggleCollapse}
        title={isCollapsed ? "展开侧边栏" : "折叠侧边栏"}
      >
        {isCollapsed ? (
          <ChevronRight className="h-4 w-4" />
        ) : (
          <ChevronLeft className="h-4 w-4" />
        )}
      </Button>

      <div className="space-y-4 py-4">
        <div className="px-3 py-2">
          {!isCollapsed && (
            <h2 className="mb-2 px-4 text-lg font-semibold tracking-tight">
              ACMOP
            </h2>
          )}
          <div className="space-y-1">
            <Button 
              variant={pathname === "/" ? "secondary" : "ghost"} 
              className={cn("w-full", isCollapsed ? "justify-center px-0" : "justify-start")} 
              asChild
              title={isCollapsed ? "Dashboard" : undefined}
            >
              <Link href="/">
                <LayoutDashboard className={cn("h-4 w-4", !isCollapsed && "mr-2")} />
                {!isCollapsed && "Dashboard"}
              </Link>
            </Button>
            <Button 
              variant={pathname === "/design" ? "secondary" : "ghost"} 
              className={cn("w-full", isCollapsed ? "justify-center px-0" : "justify-start")} 
              asChild
              title={isCollapsed ? "Design Viewer" : undefined}
            >
              <Link href="/design">
                <Box className={cn("h-4 w-4", !isCollapsed && "mr-2")} />
                {!isCollapsed && "Design Viewer"}
              </Link>
            </Button>
            <Button 
              variant={pathname === "/optimization" ? "secondary" : "ghost"} 
              className={cn("w-full", isCollapsed ? "justify-center px-0" : "justify-start")} 
              asChild
              title={isCollapsed ? "Optimization" : undefined}
            >
              <Link href="/optimization">
                <BarChart3 className={cn("h-4 w-4", !isCollapsed && "mr-2")} />
                {!isCollapsed && "Optimization"}
              </Link>
            </Button>
            <Button 
              variant={pathname === "/visualizer" ? "secondary" : "ghost"} 
              className={cn("w-full", isCollapsed ? "justify-center px-0" : "justify-start")} 
              asChild
              title={isCollapsed ? "Visualizer" : undefined}
            >
              <Link href="/visualizer">
                <Sparkles className={cn("h-4 w-4", !isCollapsed && "mr-2")} />
                {!isCollapsed && "Visualizer"}
              </Link>
            </Button>
            <Button 
              variant={pathname === "/design-visualizer" ? "secondary" : "ghost"} 
              className={cn("w-full", isCollapsed ? "justify-center px-0" : "justify-start")} 
              asChild
              title={isCollapsed ? "Design Visualizer (New)" : undefined}
            >
              <Link href="/design-visualizer">
                <Box className={cn("h-4 w-4", !isCollapsed && "mr-2")} />
                {!isCollapsed && "Design Visualizer (New)"}
              </Link>
            </Button>
            <Button 
              variant={pathname === "/design-analyzer" ? "secondary" : "ghost"} 
              className={cn("w-full", isCollapsed ? "justify-center px-0" : "justify-start")} 
              asChild
              title={isCollapsed ? "Design Analyzer" : undefined}
            >
              <Link href="/design-analyzer">
                <FileSearch className={cn("h-4 w-4", !isCollapsed && "mr-2")} />
                {!isCollapsed && "Design Analyzer"}
              </Link>
            </Button>
            <Button 
              variant={pathname === "/csv-visualizer" ? "secondary" : "ghost"} 
              className={cn("w-full", isCollapsed ? "justify-center px-0" : "justify-start")} 
              asChild
              title={isCollapsed ? "CSV Visualizer" : undefined}
            >
              <Link href="/csv-visualizer">
                <BarChart3 className={cn("h-4 w-4", !isCollapsed && "mr-2")} />
                {!isCollapsed && "CSV Visualizer"}
              </Link>
            </Button>
            <Button 
              variant={pathname === "/developer" ? "secondary" : "ghost"} 
              className={cn("w-full", isCollapsed ? "justify-center px-0" : "justify-start")} 
              asChild
              title={isCollapsed ? "开发者页面" : undefined}
            >
              <Link href="/developer">
                <Code2 className={cn("h-4 w-4", !isCollapsed && "mr-2")} />
                {!isCollapsed && "开发者页面"}
              </Link>
            </Button>
          </div>
        </div>

        {mounted && (
          <div className="px-3 py-2 border-t">
            <Button
              variant="ghost"
              className={cn("w-full", isCollapsed ? "justify-center px-0" : "justify-start")}
              onClick={toggleTheme}
              title={isCollapsed ? (theme === "dark" ? "Light Mode" : "Dark Mode") : undefined}
            >
              {theme === "dark" ? (
                <Sun className={cn("h-4 w-4", !isCollapsed && "mr-2")} />
              ) : (
                <Moon className={cn("h-4 w-4", !isCollapsed && "mr-2")} />
              )}
              {!isCollapsed && (theme === "dark" ? "Light Mode" : "Dark Mode")}
            </Button>
          </div>
        )}
      </div>
    </div>
  )
}
