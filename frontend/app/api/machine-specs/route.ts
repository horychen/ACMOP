import { NextResponse } from "next/server";

const BACKEND_URL = process.env.BACKEND_URL || "http://localhost:8000";

/**
 * 代理到 FastAPI 的 /api/machine-specs，返回来自 user_minitureMachine.py 的规格（含 _debug）。
 * 前端请求同源 /api/machine-specs 即可拿到与后端一致的数据，避免直连 8000 时的缓存或未命中。
 */
export async function GET() {
  try {
    const res = await fetch(`${BACKEND_URL}/api/machine-specs`, {
      method: "GET",
      headers: { "Content-Type": "application/json" },
      cache: "no-store",
    });

    if (!res.ok) {
      const text = await res.text();
      try {
        const err = JSON.parse(text);
        return NextResponse.json(
          { error: err.detail ?? "Failed to fetch machine specs" },
          { status: res.status }
        );
      } catch {
        return NextResponse.json(
          { error: text || "Failed to fetch machine specs" },
          { status: res.status }
        );
      }
    }

    const data = await res.json();
    return NextResponse.json(data, {
      headers: { "Cache-Control": "no-store" },
    });
  } catch (e: unknown) {
    const message = e instanceof Error ? e.message : "Request failed";
    return NextResponse.json(
      { error: message },
      { status: 500 }
    );
  }
}
