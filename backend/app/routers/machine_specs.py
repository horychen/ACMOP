import json
import os
import sys
import time
import types
from fastapi import APIRouter, HTTPException
from fastapi.responses import JSONResponse
from dataclasses import asdict, is_dataclass, fields
from fastapi.encoders import jsonable_encoder

router = APIRouter()

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
_USER_MACHINE_PATH = os.path.join(_ROOT, "backend", "codes4", "machine_geometry.py")
_WORKSPACE_MACHINE_PATH = os.path.join(r"c:\Users\lenovo\Codes\ACMOP", "backend", "codes4", "machine_geometry.py")
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

def _robust_dict(obj):
    """Deeply convert objects to dicts, manually iterating fields to avoid stale metadata issues."""
    if is_dataclass(obj):
        res = {}
        for f in fields(obj):
            try:
                # Use getattr with a sentinel to avoid AttributeError on stale fields
                val = getattr(obj, f.name, None)
                if val is not None:
                    res[f.name] = _robust_dict(val)
            except Exception:
                continue
        return res
    
    # Check by class name to catch Parameter objects from dynamic modules or utilities
    if type(obj).__name__ == "Parameter":
        return {
            "name": getattr(obj, "name", ""),
            "type": getattr(obj, "type", ""),
            "value": getattr(obj, "value", 0.0),
            "bounds": getattr(obj, "bounds", None),
            "unit": getattr(obj, "unit", "mm")
        }
    
    if isinstance(obj, dict):
        return {k: _robust_dict(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_robust_dict(v) for v in obj]
    return obj

def _specs_to_api_response(specs):
    """Build API response from MotorSpecs; using robust conversion for dynamic classes."""
    
    geometry_dict = _robust_dict(specs.geometry)
    
    # Flatten tooth_specs into geometry if it exists
    if hasattr(specs.geometry, "tooth_specs"):
        ts_dict = _robust_dict(specs.geometry.tooth_specs)
        for k, v in ts_dict.items():
            if k != "shape": # Don't overwrite the main tooth_shape if there's a conflict
                geometry_dict[k] = v

    # Manually build the structure to ensure everything is serializable
    data = {
        "geometry": geometry_dict,
        "winding": _robust_dict(specs.winding),
        "materials": _robust_dict(specs.materials),
        "targets": _robust_dict(specs.targets),
        "validations": {
            "geometry": specs.validate_inputs("geometry"),
            "materials": specs.validate_inputs("materials"),
            "winding": specs.validate_inputs("winding"),
            "targets": specs.validate_inputs("targets"),
        }
    }
    return jsonable_encoder(data)

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
    module = types.ModuleType("machine_geometry")
    module.__file__ = path
    with open(path, "r", encoding="utf-8") as f:
        source = f.read()
    _agent_log({"location": "dynamic_load", "path": path, "source_start": source[:100]})
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
        _agent_log({
            "location": "specs_init", 
            "module": specs.__class__.__module__,
            "file": getattr(sys.modules.get(specs.__class__.__module__), "__file__", "N/A"),
            "sys_path": sys.path[:5]
        })
        g = specs.geometry
        out = _specs_to_api_response(specs)
        # #region agent log
        _agent_log({
            "hypothesisId": "H1-H4", 
            "location": "machine_specs.py:get_machine_specs", 
            "message": "backend after load", 
            "data": {
                "path_used": _get_machine_file_path(), 
                "g_w_tooth": g.w_tooth.value, 
                "g_d_tooth": g.d_tooth.value, 
                "out_w_tooth": out["geometry"]["w_tooth"]["value"], 
                "out_d_tooth": out["geometry"]["d_tooth"]["value"]
            }, 
            "timestamp": time.time() * 1000
        })
        path_used = _get_machine_file_path()
        out["_debug"] = {
            "source": "workspace_machine_specs",
            "w_tooth_from_py": g.w_tooth.value,
            "path_used": path_used,
        }
        # #endregion
        return JSONResponse(content=out, headers={"Cache-Control": "no-store"})
    except Exception as e:
        import traceback
        _agent_log({"location": "API_ERROR", "error": str(e), "traceback": traceback.format_exc()})
        raise HTTPException(status_code=500, detail=str(e))
