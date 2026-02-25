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
_WORKSPACE_MACHINE_PATH = os.path.join(_ROOT, "backend", "codes4", "machine_geometry.py")
_LOG_PATH = os.path.join(_ROOT, ".cursor", "debug.log")
_DEBUG_LOG_FALLBACK = os.path.join(_ROOT, "debug_specs_log.txt")

def _agent_log(payload: dict) -> None:
    for path in (_LOG_PATH,):
        try:
            d = os.path.dirname(path)
            if d:
                os.makedirs(d, exist_ok=True)
            with open(path, "a", encoding="utf-8") as f:
                f.write(json.dumps(payload) + "\n")
                f.flush()
                os.fsync(f.fileno())
        except Exception:
            pass

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
    
    # Handle regular objects by their __dict__ if they aren't caught above
    if hasattr(obj, "__dict__") and not isinstance(obj, type):
        return {k: _robust_dict(v) for k, v in obj.__dict__.items() if not k.startswith('_')}
        
    return obj

def _specs_to_api_response(my_machine):
    """Build API response from Machine object; using robust conversion for dynamic classes."""
    
    geometry_dict = _robust_dict(my_machine.geometry)

    data = {
        "geometry": geometry_dict,
        "winding": _robust_dict(my_machine.winding),
        "materials": _robust_dict(my_machine.materials),
        "targets": _robust_dict(my_machine.target),
        "validations": {
            "geometry": [],
            "materials": [],
            "winding": [],
            "targets": [],
        }
    }
    return jsonable_encoder(data)

@router.get("/machine-specs-version")
async def get_machine_specs_version():
    """仅当前工作区后端有此接口；用于确认 8000 端口跑的是本仓库代码。"""
    return JSONResponse(
        content={"_debug": "workspace", "version": "workspace_machine_specs"},
        headers={"Cache-Control": "no-store"},
    )


@router.get("/machine-specs")
async def get_machine_specs():
    """采用新的 user_input Machine 模型替换遗留的 MotorSpecs"""
    try:
        from app.routers.debug import get_default_user_input
        from codes4.machine import Machine
        
        user_input = get_default_user_input()
        my_machine = Machine(user_input)
        my_machine.sync()
        
        out = _specs_to_api_response(my_machine)
        out["_debug"] = {
            "source": "workspace_machine_specs",
        }
        return JSONResponse(content=out, headers={"Cache-Control": "no-store"})
    except Exception as e:
        import traceback
        _agent_log({"location": "API_ERROR", "error": str(e), "traceback": traceback.format_exc()})
        raise HTTPException(status_code=500, detail=str(e))
