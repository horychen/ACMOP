"use client"

import Link from "next/link"
import { usePathname } from "next/navigation"
import { cn } from "@/lib/utils"
import { Button } from "@/components/ui/button"
import { ScrollArea } from "@/components/ui/scroll-area"
import { LayoutDashboard, Settings, BarChart3, Box } from "lucide-react"

interface SidebarProps extends React.HTMLAttributes<HTMLDivElement> {}

export function Sidebar({ className }: SidebarProps) {
  const pathname = usePathname()

  return (
    <div className={cn("pb-12 w-64 border-r min-h-screen", className)}>
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
          </div>
        </div>
      </div>
    </div>
  )
}
