import { NextRequest, NextResponse } from "next/server";

const BACKEND_URL = process.env.BACKEND_URL || "http://localhost:8000";

export async function GET(request: NextRequest) {
  try {
    const searchParams = request.nextUrl.searchParams;
    const folderName = searchParams.get("folderName");
    const path2SwarmData = searchParams.get("path2SwarmData");
    const path2MachineDesignerFull = searchParams.get("path2MachineDesignerFull");

    const params = new URLSearchParams();
    if (folderName) params.append("folderName", folderName);
    if (path2SwarmData) params.append("path2SwarmData", path2SwarmData);
    if (path2MachineDesignerFull) params.append("path2MachineDesignerFull", path2MachineDesignerFull);

    const response = await fetch(`${BACKEND_URL}/api/acmopv2/pareto-front?${params.toString()}`, {
      method: "GET",
      headers: {
        "Content-Type": "application/json",
      },
    });

    if (!response.ok) {
      const error = await response.json();
      return NextResponse.json(
        { error: error.detail || "获取Pareto前沿数据失败" },
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

