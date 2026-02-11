"""
加载项目根目录下的用户配置文件 acmop.config.json。
用于前后端统一指定虚拟环境（默认 conda 环境名 acmop）及后端地址。
"""

import json
import os

# 项目根目录：backend/app -> backend -> 项目根
_APP_DIR = os.path.dirname(os.path.abspath(__file__))
_BACKEND_DIR = os.path.dirname(_APP_DIR)
_PROJECT_ROOT = os.path.dirname(_BACKEND_DIR)

CONFIG_FILENAME = "acmop.config.json"
_DEFAULT_CONFIG = {
    "backend": {
        "virtualEnv": {"type": "conda", "name": "acmop"},
        "url": "http://localhost:8000",
    },
    "frontend": {
        "backendUrl": "http://localhost:8000",
    },
}


def get_project_root() -> str:
    """返回项目根目录绝对路径。"""
    return _PROJECT_ROOT


def load_acmop_config() -> dict:
    """
    加载 acmop.config.json。若文件不存在或解析失败，返回默认配置。
    默认后端虚拟环境为 conda，名称为 acmop。
    """
    path = os.path.join(_PROJECT_ROOT, CONFIG_FILENAME)
    if not os.path.isfile(path):
        return _DEFAULT_CONFIG.copy()
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        # 深合并默认值，避免缺键
        def merge(base: dict, override: dict) -> dict:
            out = base.copy()
            for k, v in override.items():
                if k in out and isinstance(out[k], dict) and isinstance(v, dict):
                    out[k] = merge(out[k], v)
                else:
                    out[k] = v
            return out
        return merge(_DEFAULT_CONFIG, data)
    except (json.JSONDecodeError, OSError):
        return _DEFAULT_CONFIG.copy()


def get_backend_virtual_env() -> dict:
    """返回 backend.virtualEnv 配置，默认 { \"type\": \"conda\", \"name\": \"acmop\" }。"""
    config = load_acmop_config()
    return config.get("backend", {}).get("virtualEnv", _DEFAULT_CONFIG["backend"]["virtualEnv"])


def get_frontend_backend_url() -> str:
    """返回 frontend.backendUrl，供前端或代理使用。"""
    config = load_acmop_config()
    return config.get("frontend", {}).get("backendUrl", _DEFAULT_CONFIG["frontend"]["backendUrl"])
