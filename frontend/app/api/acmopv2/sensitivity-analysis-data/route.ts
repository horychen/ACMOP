import { NextRequest, NextResponse } from "next/server";
import fs from "fs";
import path from "path";

export async function GET(request: NextRequest) {
  try {
    const searchParams = request.nextUrl.searchParams;
    const folderPath = searchParams.get("folderPath");

    if (!folderPath) {
      return NextResponse.json(
        { error: "folderPath 参数是必需的" },
        { status: 400 }
      );
    }

    // 构建 SwarmData.json 的完整路径
    // folderPath 应该是相对于 backend_v2/output 目录的路径，例如: "_default/TIA_prototype_sensitivity_analysis"
    const swarmDataPath = path.join(process.cwd(), "../backend_v2/output", folderPath, "SwarmData.json");

    // 检查文件是否存在
    if (!fs.existsSync(swarmDataPath)) {
      return NextResponse.json(
        { error: `SwarmData.json 文件未找到: ${swarmDataPath}` },
        { status: 404 }
      );
    }

    // 读取文件内容
    const fileContent = fs.readFileSync(swarmDataPath, "utf-8");
    let data;

    try {
      data = JSON.parse(fileContent);
    } catch (parseError) {
      // 如果文件以逗号开头（某些格式），尝试处理
      const trimmedContent = fileContent.trim();
      if (trimmedContent.startsWith(",")) {
        data = JSON.parse("{" + trimmedContent.substring(1) + "}");
      } else {
        throw parseError;
      }
    }

    return NextResponse.json(data);
  } catch (error: any) {
    console.error("读取敏感性分析数据时出错:", error);
    return NextResponse.json(
      { error: error.message || "读取数据失败" },
      { status: 500 }
    );
  }
}
