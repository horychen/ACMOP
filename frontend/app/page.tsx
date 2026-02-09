"use client"

import { useEffect } from "react";
import { useRouter } from "next/navigation";

export default function Home() {
  const router = useRouter();

  useEffect(() => {
    router.replace("/optimization");
  }, [router]);

  return (
    <div className="container mx-auto py-8 flex items-center justify-center min-h-[200px]">
      <p className="text-muted-foreground">正在跳转到优化页面…</p>
    </div>
  );
}

export const dynamic = "force-dynamic";
