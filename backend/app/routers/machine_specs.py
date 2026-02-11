import json
import os
import sys
import time
import types
from fastapi import APIRouter, HTTPException
from fastapi.responses import JSONResponse
from dataclasses import asdict

router = APIRouter()
# Project root: backend/app/routers -> backend/app -> backend -> project root
_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
_USER_MACHINE_PATH = os.path.join(_ROOT, "backend", "codes4", "user_minitureMachine.py")
_WORKSPACE_MACHINE_PATH = os.path.join(r"c:\Users\lenovo\Codes\ACMOP", "backend", "codes4", "user_minitureMachine.py")
_LOG_PATH = os.path.join(_ROOT, ".cursor", "debug.log")
_DEBUG_LOG_FALLBACK = r"c:\Users\lenovo\Codes\ACMOP\.cursor\debug.log"

def _agent_log(payload: dict) -> None:
    codes4_log = os.path.join(os.path.dirname(_USER_MACHINE_PATH), "debug_specs_log.txt")
    for path in (codes4_log, _LOG_PATH, _DEBUG_LOG_FALLBACK):
        try:
            d = os.path.dirname(path)
            if d:
                os.makedirs(d, exist_ok=True)
            with open(path, "a", encoding="utf-8") as f:
                f.write(json.dumps(payload) + "\n")
            break
        except Exception:
            continue

def _specs_to_api_response(specs):
    """Build API response from MotorSpecs; geometry/winding use asdict so all fields match user_minitureMachine.py."""
    g = specs.geometry
    w = specs.winding
    stack_opts = getattr(g, "stack_length_options", [6.0, 10.0, 16.0])
    stack_length = float(stack_opts[len(stack_opts) // 2] if stack_opts else 10.0)
    geom = asdict(g)
    geom["stack_length"] = stack_length
    wind = asdict(w)
    return {
        "geometry": geom,
        "winding": wind,
        "materials": asdict(specs.materials),
        "targets": asdict(specs.targets),
    }

def _get_machine_file_path():
    """Use workspace path if it exists, else computed path, so we always read the file the user edits."""
    if os.path.isfile(_WORKSPACE_MACHINE_PATH):
        return _WORKSPACE_MACHINE_PATH
    return _USER_MACHINE_PATH

def _load_module_from_file():
    """
    每次请求从磁盘读取 user_minitureMachine.py 源码并执行，避免 .pyc 缓存导致修改后仍读到旧默认值。
    """
    path = _get_machine_file_path()
    codes4_dir = os.path.dirname(path)
    if codes4_dir not in sys.path:
        sys.path.insert(0, codes4_dir)
    module = types.ModuleType("user_minitureMachine")
    module.__file__ = path
    with open(path, "r", encoding="utf-8") as f:
        source = f.read()
    code = compile(source, path, "exec")
    exec(code, module.__dict__)
    return module

@router.get("/machine-specs-version")
async def get_machine_specs_version():
    """仅当前工作区后端有此接口；用于确认 8000 端口跑的是本仓库代码。"""
    return JSONResponse(
        content={"_debug": "workspace", "version": "workspace_machine_specs"},
        headers={"Cache-Control": "no-store"},
    )


@router.get("/machine-specs")
async def get_machine_specs():
    """
    每次请求时从工作区 backend/codes4/user_minitureMachine.py 文件路径加载模块，返回最新 MotorSpecs。
    """
    try:
        # #region agent log
        _path_used = _get_machine_file_path()
        _agent_log({"hypothesisId": "H1-H4", "location": "machine_specs.py:get_machine_specs", "message": "backend before load", "data": {"path_used": _path_used, "path_exists": os.path.isfile(_path_used), "workspace_path": _WORKSPACE_MACHINE_PATH}, "timestamp": time.time() * 1000})
        # #endregion
        um_module = _load_module_from_file()
        specs = um_module.MotorSpecs()
        g = specs.geometry
        out = _specs_to_api_response(specs)
        # #region agent log
        _agent_log({"hypothesisId": "H1-H4", "location": "machine_specs.py:get_machine_specs", "message": "backend after load", "data": {"path_used": _get_machine_file_path(), "g_tooth_width": g.tooth_width, "g_tooth_depth": g.tooth_depth, "out_tooth_width": out["geometry"]["tooth_width"], "out_tooth_depth": out["geometry"]["tooth_depth"]}, "timestamp": time.time() * 1000})
        path_used = _get_machine_file_path()
        out["_debug"] = {
            "source": "workspace_machine_specs",
            "tooth_width_from_py": g.tooth_width,
            "path_used": path_used,
        }
        # #endregion
        return JSONResponse(content=out, headers={"Cache-Control": "no-store"})
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
