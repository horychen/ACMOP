import { NextResponse } from "next/server";

const BACKEND_URL = process.env.BACKEND_URL || "http://localhost:8000";

/**
 * 获取用户配置（acmop.config.json），用于显示前端与后端所采用的虚拟环境等。
 * 默认后端为 conda 环境 "acmop"。
 */
export async function GET() {
  try {
    const response = await fetch(`${BACKEND_URL}/api/acmopv2/config`, {
      cache: "no-store",
    });
    if (!response.ok) {
      return NextResponse.json(
        { error: "获取配置失败" },
        { status: response.status }
      );
    }
    const data = await response.json();
    return NextResponse.json(data);
  } catch (error: any) {
    return NextResponse.json(
      { error: error.message || "请求失败" },
      { status: 500 }
    );
  }
}
