"""
ACMOP v2 API - 用于前端可视化的后端API
这个模块提供了两个主要功能：
1. 开发者模式：根据输入参数生成JSON配置文件
2. 可视化模式：提供敏感性分析和优化结果数据
"""

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import Dict, List, Any, Optional
import json
import os
from datetime import datetime

from app.config_loader import load_acmop_config

router = APIRouter(prefix="/api/acmopv2", tags=["ACMOP v2"])


# ==================== 用户配置 API ====================

@router.get("/config")
async def get_user_config() -> Dict[str, Any]:
    """
    获取用户配置文件 acmop.config.json 内容。
    用于显式指定前端与后端所采用的虚拟环境，默认后端为 conda 环境 \"acmop\"。
    """
    return load_acmop_config()


# ==================== 辅助函数：路径处理 ====================

def get_codes4_path():
    """获取codes4目录路径"""
    current_dir = os.path.dirname(os.path.abspath(__file__))  # backend/app/routers
    backend_dir = os.path.dirname(os.path.dirname(current_dir))  # backend
    codes4_dir = os.path.join(backend_dir, "codes4")
    return codes4_dir

def get_default_dir():
    """获取_default目录路径"""
    current_dir = os.path.dirname(os.path.abspath(__file__))  # backend/app/routers
    backend_dir = os.path.dirname(os.path.dirname(current_dir))  # backend
    default_dir = os.path.join(backend_dir, "_default")
    return default_dir


# ==================== 数据模型 ====================

class DesignParameters(BaseModel):
    """设计参数模型"""
    power: Optional[float] = None  # 额定功率 (kW)
    speed: Optional[float] = None  # 额定转速 (RPM)
    voltage: Optional[float] = None  # 额定电压 (V)
    slots: Optional[int] = None  # 槽数
    poles: Optional[int] = None  # 极数
    polePairs: Optional[int] = None  # 极对数
    # 可以在这里添加更多参数
    
    class Config:
        extra = "allow"  # 允许额外字段


class ParameterConfig(BaseModel):
    """单个参数配置模型"""
    name: str
    type: str  # 'fixed', 'free', 'derived'
    value: Optional[Any] = None
    bounds: Optional[List[float]] = None
    calc: Optional[str] = None  # Python lambda code as string
    unit: Optional[str] = "mm"
    comment: Optional[str] = None
    sensitivityAnalysis: Optional[bool] = False

class ConfigRequest(BaseModel):
    """配置请求模型"""
    parameters: List[ParameterConfig]


class AnalysisRequest(BaseModel):
    """分析请求模型"""
    config: Dict[str, Any]


# ==================== 开发者模式 API ====================

@router.get("/all-parameters")
async def get_all_parameters() -> Dict[str, Any]:
    """
    获取所有参数的默认值和元数据（包括fixed、free、derived）
    
    返回完整的参数列表，用于前端预填充
    """
    try:
        from codes4.machine_design_guide import Modern_Machine_Designer
        import inspect
        import traceback
        
        # 创建临时实例以获取参数信息
        try:
            designer = Modern_Machine_Designer()
        except Exception as e:
            raise HTTPException(
                status_code=500, 
                detail=f"创建 Modern_Machine_Designer 实例失败: {str(e)}\n{traceback.format_exc()}"
            )
        
        # 获取所有参数
        try:
            all_params = designer.get_parameter_fields()
        except Exception as e:
            raise HTTPException(
                status_code=500,
                detail=f"获取参数列表失败: {str(e)}\n{traceback.format_exc()}"
            )
        
        # 按类型分组
        try:
            fixed_params = designer.get_parameters_by_type('fixed')
            free_params = designer.get_parameters_by_type('free')
            derived_params = designer.get_parameters_by_type('derived')
        except Exception as e:
            raise HTTPException(
                status_code=500,
                detail=f"按类型分组参数失败: {str(e)}\n{traceback.format_exc()}"
            )
        
        def format_calc(calc_func, args_list, param_name=""):
            """格式化calc函数为字符串，提取完整的源代码"""
            if calc_func is None:
                return None
            try:
                # 尝试获取lambda函数的源代码
                if hasattr(calc_func, '__code__'):
                    try:
                        # 获取完整的源代码
                        source_lines = inspect.getsourcelines(calc_func)
                        if source_lines and len(source_lines) > 0:
                            source = ''.join(source_lines[0])
                            # 提取lambda表达式部分
                            if 'lambda' in source:
                                # 找到lambda开始的位置
                                lambda_start = source.find('lambda')
                                if lambda_start != -1:
                                    # 提取从lambda开始到行尾的内容
                                    lambda_expr = source[lambda_start:].strip()
                                    # 移除可能的换行，但保留代码结构
                                    lines = lambda_expr.split('\n')
                                    if len(lines) > 0:
                                        # 取第一行，这通常包含完整的lambda
                                        first_line = lines[0].strip()
                                        # 如果第一行以逗号或右括号结束，去掉
                                        if first_line.endswith(','):
                                            first_line = first_line[:-1]
                                        # 尝试找到完整的lambda表达式
                                        # 查找匹配的括号
                                        paren_count = 0
                                        result = []
                                        for char in first_line:
                                            result.append(char)
                                            if char == '(':
                                                paren_count += 1
                                            elif char == ')':
                                                paren_count -= 1
                                                if paren_count == 0 and 'lambda' in ''.join(result):
                                                    break
                                        
                                        lambda_expr = ''.join(result).strip()
                                        # 移除末尾的逗号
                                        if lambda_expr.endswith(','):
                                            lambda_expr = lambda_expr[:-1]
                                        
                                        return lambda_expr
                                    return first_line.strip()
                    except (OSError, TypeError, ValueError, AttributeError) as e:
                        # 如果无法获取源代码，尝试其他方法
                        pass
                
                # 如果无法获取源代码，尝试从args构建lambda表达式
                if args_list and len(args_list) > 0:
                    arg_names = []
                    for arg in args_list:
                        try:
                            if hasattr(arg, 'name'):
                                arg_names.append(arg.name)
                            elif isinstance(arg, (int, float, str)):
                                arg_names.append(str(arg))
                            else:
                                arg_names.append(str(arg))
                        except:
                            arg_names.append("arg")
                    return f"lambda {', '.join(arg_names)}: ..."
                
                return str(calc_func)
            except Exception as e:
                return f"lambda: ...  # Error: {str(e)}"
        
        def format_args(args_list):
            """格式化args列表为参数名称列表"""
            if args_list is None:
                return []
            result = []
            for arg in args_list:
                try:
                    if hasattr(arg, 'name'):
                        result.append(arg.name)
                    elif isinstance(arg, (int, float, str, list)):
                        result.append(str(arg))
                    else:
                        result.append(str(arg))
                except:
                    result.append("arg")
            return result
        
        parameters_info = {
            "fixed": [],
            "free": [],
            "derived": []
        }
        
        # Fixed参数
        for param_name, param_obj in fixed_params.items():
            try:
                param_info = {
                    "name": param_name,
                    "displayName": getattr(param_obj, 'name', param_name),
                    "type": "fixed",
                    "unit": param_obj.unit or "mm",
                    "value": param_obj.value,
                    "comment": getattr(param_obj, 'comment', None) or getattr(param_obj, 'comnment', None) or ""
                }
                parameters_info["fixed"].append(param_info)
            except Exception as e:
                # 如果某个参数处理失败，记录错误但继续处理其他参数
                print(f"Warning: Failed to process fixed parameter '{param_name}': {e}")
                continue
        
        # Free参数
        for param_name, param_obj in free_params.items():
            try:
                param_info = {
                    "name": param_name,
                    "displayName": getattr(param_obj, 'name', param_name),
                    "type": "free",
                    "unit": param_obj.unit or "mm",
                    "value": param_obj.value,
                    "bounds": param_obj.bounds,
                    "calc_bounds": format_calc(getattr(param_obj, 'calc_bounds', None), getattr(param_obj, 'args', None), param_name),
                    "args": format_args(getattr(param_obj, 'args', None)),
                    "comment": getattr(param_obj, 'comment', None) or getattr(param_obj, 'comnment', None) or ""
                }
                parameters_info["free"].append(param_info)
            except Exception as e:
                # 如果某个参数处理失败，记录错误但继续处理其他参数
                print(f"Warning: Failed to process free parameter '{param_name}': {e}")
                continue
        
        # Derived参数
        for param_name, param_obj in derived_params.items():
            try:
                param_info = {
                    "name": param_name,
                    "displayName": getattr(param_obj, 'name', param_name),
                    "type": "derived",
                    "unit": param_obj.unit or "mm",
                    "value": param_obj.value,
                    "calc": format_calc(getattr(param_obj, 'calc', None), getattr(param_obj, 'args', None), param_name),
                    "args": format_args(getattr(param_obj, 'args', None)),
                    "comment": getattr(param_obj, 'comment', None) or getattr(param_obj, 'comnment', None) or ""
                }
                parameters_info["derived"].append(param_info)
            except Exception as e:
                # 如果某个参数处理失败，记录错误但继续处理其他参数
                print(f"Warning: Failed to process derived parameter '{param_name}': {e}")
                continue
        
        return {
            "parameters": parameters_info,
            "machine_class": designer.machine_class,
            "summary": {
                "fixed_count": len(parameters_info["fixed"]),
                "free_count": len(parameters_info["free"]),
                "derived_count": len(parameters_info["derived"]),
                "total_count": len(all_params)
            }
        }
        
    except HTTPException:
        raise
    except Exception as e:
        import traceback
        error_detail = f"获取参数列表时出错: {str(e)}\n{traceback.format_exc()}"
        print(error_detail)  # 在服务器日志中打印详细错误
        raise HTTPException(status_code=500, detail=f"获取参数列表时出错: {str(e)}")


@router.get("/fixed-parameters")
async def get_fixed_parameters() -> Dict[str, Any]:
    """
    获取所有fixed参数的元数据（向后兼容）
    """
    try:
        all_params = await get_all_parameters()
        return {
            "parameters": all_params["parameters"]["fixed"],
            "total_count": len(all_params["parameters"]["fixed"])
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"获取fixed参数列表时出错: {str(e)}")


@router.post("/generate")
async def generate_config(request: ConfigRequest) -> Dict[str, Any]:
    """
    根据输入参数生成JSON配置文件
    
    这个函数将调用codes4中的相关函数来生成完整的配置JSON
    目前返回一个示例结构，您可以根据实际需求修改
    """
    try:
        param_configs = request.parameters
        
        # 导入必要的模块
        from codes4.machine_design_guide import Modern_Machine_Designer, Parameter
        import ast
        
        # 创建Modern_Machine_Designer实例
        designer = Modern_Machine_Designer()
        
        # 存储用户配置的参数
        user_params = {}
        sensitivity_params = []
        
        # 处理每个参数配置
        for param_config in param_configs:
            param_name = param_config.name
            
            # 创建新的Parameter对象
            if param_config.type == "fixed":
                param = Parameter(
                    name=param_config.name,
                    type="fixed",
                    value=param_config.value,
                    unit=param_config.unit or "mm",
                    comnment=param_config.comment
                )
                # 如果参数已存在，更新它；否则需要动态添加（这里简化处理）
                existing_param = designer.get_parameter(param_name)
                if existing_param:
                    designer.set_parameter_value(param_name, param_config.value)
                user_params[param_name] = {
                    "type": "fixed",
                    "value": param_config.value,
                    "unit": param_config.unit or "mm"
                }
                
            elif param_config.type == "free":
                bounds = param_config.bounds if param_config.bounds else [0, 100]
                param = Parameter(
                    name=param_config.name,
                    type="free",
                    value=param_config.value,
                    bounds=bounds,
                    unit=param_config.unit or "mm",
                    comnment=param_config.comment
                )
                user_params[param_name] = {
                    "type": "free",
                    "value": param_config.value,
                    "bounds": bounds,
                    "unit": param_config.unit or "mm"
                }
                if param_config.sensitivityAnalysis:
                    sensitivity_params.append(param_name)
                    
            elif param_config.type == "derived":
                # 尝试解析calc函数
                calc_func = None
                if param_config.calc:
                    try:
                        # 尝试将字符串转换为lambda函数
                        # 注意：这里需要安全地执行用户代码
                        # 实际应用中应该使用更安全的方法
                        calc_func = eval(param_config.calc)  # 仅用于演示，生产环境需要更安全的实现
                    except Exception as e:
                        raise HTTPException(
                            status_code=400,
                            detail=f"参数 {param_name} 的计算函数解析失败: {str(e)}"
                        )
                
                param = Parameter(
                    name=param_config.name,
                    type="derived",
                    calc=calc_func,
                    unit=param_config.unit or "mm",
                    comnment=param_config.comment
                )
                user_params[param_name] = {
                    "type": "derived",
                    "calc": param_config.calc,
                    "unit": param_config.unit or "mm"
                }
        
        # 构建配置JSON
        config = {
            "metadata": {
                "version": "2.0",
                "created_at": datetime.now().isoformat(),
                "description": "ACMOP设计配置文件",
                "machine_class": designer.machine_class
            },
            "parameters": user_params,
            "fixed_parameters": {
                name: info["value"]
                for name, info in user_params.items()
                if info["type"] == "fixed"
            },
            "free_parameters": {
                name: {
                    "value": info["value"],
                    "bounds": info["bounds"]
                }
                for name, info in user_params.items()
                if info["type"] == "free"
            },
            "derived_parameters": {
                name: {
                    "calc": info["calc"]
                }
                for name, info in user_params.items()
                if info["type"] == "derived"
            },
            "winding_parameters": {
                "phase_number_m": designer.wily.phase_number_m,
                "stator_slot_number_Qs": designer.wily.stator_slot_number_Qs,
                "pole_pair_number_p": designer.wily.pole_pair_number_p,
                "suspension_pole_pair_number_ps": designer.wily.suspension_pole_pair_number_ps
            },
            "analysis_config": {
                "sensitivity_analysis": {
                    "enabled": len(sensitivity_params) > 0,
                    "parameters": sensitivity_params,
                    "ranges": {}
                },
                "optimization": {
                    "enabled": len([p for p in user_params.values() if p["type"] == "free"]) > 0,
                    "objectives": [],
                    "constraints": []
                }
            }
        }
        
        return config
        
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"生成配置时出错: {str(e)}")


# ==================== 可视化模式 API ====================

@router.post("/sensitivity-analysis")
async def get_sensitivity_analysis(request: AnalysisRequest) -> Dict[str, Any]:
    """
    获取敏感性分析结果
    
    根据提供的配置，返回敏感性分析数据用于可视化
    """
    try:
        config = request.config
        
        # TODO: 在这里调用实际的敏感性分析函数
        # 例如：
        # from codes4.acmop import run_sensitivity_analysis
        # results = run_sensitivity_analysis(config)
        
        # 目前返回示例数据
        parameters = config.get("design_parameters", {})
        param_names = list(parameters.keys())[:5]  # 取前5个参数作为示例
        
        sensitivity_data = []
        for param_name in param_names:
            # 生成示例敏感性数据
            base_value = parameters.get(param_name, 0)
            if isinstance(base_value, (int, float)) and base_value > 0:
                for i in range(-5, 6):
                    value = base_value * (1 + i * 0.1)
                    sensitivity_data.append({
                        "parameter": param_name,
                        "value": value,
                        "efficiency": 85 + (i * 0.5) + (hash(param_name) % 10) / 10,
                        "torque": 100 + (i * 2) + (hash(param_name) % 20) / 10,
                        "power_loss": 500 - (i * 10) + (hash(param_name) % 50) / 10
                    })
        
        return {
            "parameters": param_names,
            "sensitivityData": sensitivity_data,
            "summary": {
                "most_sensitive": param_names[0] if param_names else None,
                "analysis_date": datetime.now().isoformat()
            }
        }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"获取敏感性分析数据时出错: {str(e)}")


@router.get("/default-machine-designer")
async def get_default_machine_designer() -> Dict[str, Any]:
    """
    获取默认的 machine_designer.json 配置
    
    返回完整的默认配置，用于前端初始化
    """
    try:
        # 获取codes4目录路径
        codes4_dir = get_codes4_path()
        json_path = os.path.join(codes4_dir, "machine_designer.json")
        
        if not os.path.exists(json_path):
            raise HTTPException(
                status_code=404, 
                detail=f"默认配置文件不存在: {json_path}"
            )
        
        with open(json_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        return data
    except json.JSONDecodeError as e:
        raise HTTPException(
            status_code=500, 
            detail=f"JSON解析失败: {str(e)}"
        )
    except Exception as e:
        raise HTTPException(
            status_code=500, 
            detail=f"读取默认配置失败: {str(e)}"
        )


@router.get("/machine-designer-full")
async def get_machine_designer_full() -> Dict[str, Any]:
    """
    获取 machine_designer_full.json 配置
    
    返回完整的配置，包括所有元数据和几何信息
    """
    try:
        # 获取codes4目录路径
        codes4_dir = get_codes4_path()
        json_path = os.path.join(codes4_dir, "machine_designer_full.json")
        
        if not os.path.exists(json_path):
            raise HTTPException(
                status_code=404, 
                detail=f"完整配置文件不存在: {json_path}"
            )
        
        with open(json_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        return data
    except json.JSONDecodeError as e:
        raise HTTPException(
            status_code=500, 
            detail=f"JSON解析失败: {str(e)}"
        )
    except Exception as e:
        raise HTTPException(
            status_code=500, 
            detail=f"读取完整配置失败: {str(e)}"
        )


@router.post("/optimization-results")
async def get_optimization_results(request: AnalysisRequest) -> Dict[str, Any]:
    """
    获取优化结果
    
    根据提供的配置，返回优化结果数据用于可视化
    """
    try:
        config = request.config
        
        # TODO: 在这里调用实际的优化函数
        # 例如：
        # from codes4.acmop import run_optimization
        # results = run_optimization(config)
        
        # 目前返回示例数据
        # 生成示例Pareto前沿数据
        pareto_front = []
        for i in range(20):
            pareto_front.append({
                "objective1": 80 + i * 0.5 + (i % 3) * 0.2,  # 效率
                "objective2": 100 - i * 1.5 + (i % 2) * 0.5,  # 成本
                "objective3": 50 + i * 0.3 + (i % 4) * 0.1,  # 体积
                "solution_id": f"sol_{i}"
            })
        
        # 生成示例优化历史
        optimization_history = []
        for i in range(50):
            optimization_history.append({
                "iteration": i,
                "bestValue": 85 - i * 0.1 + (i % 5) * 0.05,
                "averageValue": 80 - i * 0.08 + (i % 5) * 0.03
            })
        
        # 生成示例最佳解
        best_solutions = [
            {
                "parameters": {
                    "power": config.get("design_parameters", {}).get("power", 50),
                    "slots": config.get("design_parameters", {}).get("slots", 12),
                    "poles": config.get("design_parameters", {}).get("poles", 4)
                },
                "performance": {
                    "efficiency": 92.5,
                    "torque": 150.3,
                    "power_loss": 350.2
                },
                "rank": 1
            },
            {
                "parameters": {
                    "power": config.get("design_parameters", {}).get("power", 50) * 0.95,
                    "slots": config.get("design_parameters", {}).get("slots", 12),
                    "poles": config.get("design_parameters", {}).get("poles", 4)
                },
                "performance": {
                    "efficiency": 91.8,
                    "torque": 148.7,
                    "power_loss": 380.5
                },
                "rank": 2
            }
        ]
        
        return {
            "paretoFront": pareto_front,
            "history": optimization_history,
            "bestSolutions": best_solutions,
            "objectives": ["效率 (%)", "成本", "体积"],
            "optimization_date": datetime.now().isoformat()
        }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"获取优化结果时出错: {str(e)}")


@router.get("/list-optimization-folders")
async def list_optimization_folders() -> Dict[str, Any]:
    """
    列出 _default 目录下的所有优化文件夹
    
    Returns:
        文件夹列表
    """
    try:
        default_dir = get_default_dir()
        
        if not os.path.exists(default_dir):
            return {"folders": []}
        
        # 列出所有子文件夹
        folders = []
        for item in os.listdir(default_dir):
            item_path = os.path.join(default_dir, item)
            if os.path.isdir(item_path):
                # 检查是否有 SwarmData.json
                swarm_data_path = os.path.join(item_path, "SwarmData.json")
                if os.path.exists(swarm_data_path):
                    folders.append(item)
        
        return {"folders": sorted(folders)}
        
    except Exception as e:
        import traceback
        raise HTTPException(
            status_code=500,
            detail=f"列出优化文件夹时出错: {str(e)}\n{traceback.format_exc()}"
        )


@router.get("/pareto-front")
async def get_pareto_front(
    folderName: str = None,
    path2SwarmData: str = None,
    path2MachineDesignerFull: str = None
) -> Dict[str, Any]:
    """
    获取Pareto前沿数据和优化配置
    
    Args:
        folderName: 优化文件夹名称（在_default目录下）
        path2SwarmData: SwarmData.json文件路径（相对于codes4目录）
        path2MachineDesignerFull: machine_designer_full.json文件路径（相对于codes4目录）
    
    Returns:
        Pareto前沿个体列表、优化配置等信息
    """
    try:
        codes4_dir = get_codes4_path()
        default_dir = get_default_dir()
        
        # 初始化变量
        machine_designer_data = {}
        select_fea_config_dict = None
        path2FEACsv = ""
        swarm_data_path = ""
        
        # 如果提供了文件夹名称，优先从该文件夹读取
        if folderName:
            folder_path = os.path.join(default_dir, folderName)
            if not os.path.exists(folder_path):
                raise HTTPException(
                    status_code=404,
                    detail=f"文件夹不存在: {folderName}"
                )
            
            # 读取该文件夹中的 machine_designer_full.json（如果存在）
            machine_designer_path_in_folder = os.path.join(folder_path, "machine_designer_full.json")
            if os.path.exists(machine_designer_path_in_folder):
                with open(machine_designer_path_in_folder, 'r', encoding='utf-8') as f:
                    machine_designer_data = json.load(f)
                select_fea_config_dict = machine_designer_data.get("select_fea_config_dict")
                # 获取 path2FEACsv，如果不存在则使用默认路径
                path2FEACsv = machine_designer_data.get("path2FEACsv")
                if not path2FEACsv:
                    # 检查 csv 目录是否存在
                    csv_dir = os.path.join(folder_path, "csv")
                    if os.path.exists(csv_dir):
                        path2FEACsv = csv_dir
            else:
                # 如果没有 machine_designer_full.json，使用默认路径
                csv_dir = os.path.join(folder_path, "csv")
                if os.path.exists(csv_dir):
                    path2FEACsv = csv_dir
            
            # 读取 SwarmData.json
            swarm_data_path = os.path.join(folder_path, "SwarmData.json")
        else:
            # 读取 machine_designer_full.json
            if path2MachineDesignerFull:
                machine_designer_path = os.path.join(codes4_dir, path2MachineDesignerFull)
            else:
                machine_designer_path = os.path.join(codes4_dir, "machine_designer_full.json")
            
            if not os.path.exists(machine_designer_path):
                raise HTTPException(
                    status_code=404,
                    detail=f"machine_designer_full.json 文件不存在: {machine_designer_path}"
                )
            
            with open(machine_designer_path, 'r', encoding='utf-8') as f:
                machine_designer_data = json.load(f)
            
            # 获取 select_fea_config_dict
            select_fea_config_dict = machine_designer_data.get("select_fea_config_dict", None)
            
            # 读取 SwarmData.json
            if path2SwarmData:
                swarm_data_path = os.path.join(codes4_dir, path2SwarmData, "SwarmData.json")
            else:
                path2SwarmData_rel = machine_designer_data.get("path2SwarmData", "../_default/SPMSM/")
                swarm_data_path = os.path.join(codes4_dir, path2SwarmData_rel, "SwarmData.json")
                path2FEACsv = machine_designer_data.get("path2FEACsv", "")
        
        # 获取 select_fea_config_dict（如果还没有）
        if not select_fea_config_dict:
            select_fea_config_dict = machine_designer_data.get("select_fea_config_dict", None)
            if not select_fea_config_dict:
                # 尝试使用默认配置
                select_fea_config_dict = "#0213 JMAG Bearingless Sub-hamonics"
        
        # 读取 machine_simulation.json 获取优化配置
        simulation_config_path = os.path.join(codes4_dir, "machine_simulation.json")
        if not os.path.exists(simulation_config_path):
            raise HTTPException(
                status_code=404,
                detail=f"machine_simulation.json 文件不存在: {simulation_config_path}"
            )
        
        with open(simulation_config_path, 'r', encoding='utf-8') as f:
            simulation_configs = json.load(f)
        
        # 获取对应的优化配置
        fea_config_dict = simulation_configs.get(select_fea_config_dict, {})
        if not fea_config_dict:
            raise HTTPException(
                status_code=404,
                detail=f"在 machine_simulation.json 中找不到配置: {select_fea_config_dict}"
            )
        
        # 提取 moo.* 配置
        moo_config = {
            k: v for k, v in fea_config_dict.items() 
            if k.startswith("moo.")
        }
        
        if not os.path.exists(swarm_data_path):
            # 如果文件不存在，返回空数据
            result = {
                "paretoFront": [],
                "allIndividuals": [],
                "mooConfig": moo_config,
                "objectives": [moo_config.get("moo.fitness_OA"), moo_config.get("moo.fitness_OB"), moo_config.get("moo.fitness_OC")],
                "select_fea_config_dict": select_fea_config_dict,
                "message": f"SwarmData.json 文件不存在: {swarm_data_path}"
            }
            return convert_numpy_types(result)
        
        with open(swarm_data_path, 'r', encoding='utf-8') as f:
            swarm_data = json.load(f)
        
        # 计算 Pareto 前沿 - 直接实现，避免导入 machine_design_guide
        pg = None
        try:
            import pygmo as pg
        except ImportError:
            # 如果没有 pygmo，使用简单的非支配排序实现
            pass
        
        # 提取目标函数值 f1, f2, f3
        fits = []
        swarm_data_items = list(swarm_data.items())
        
        for key, individual_data in swarm_data_items:
            # 提取 f1, f2, f3
            f1 = individual_data.get('f1', 0.0)
            f2 = individual_data.get('f2', 0.0)
            f3 = individual_data.get('f3', 0.0)
            
            # 处理空字典或 None 值
            if isinstance(f1, dict) and len(f1) == 0:
                f1 = 0.0
            if isinstance(f2, dict) and len(f2) == 0:
                f2 = 0.0
            if isinstance(f3, dict) and len(f3) == 0:
                f3 = 0.0
            
            f1 = float(f1) if f1 is not None else 0.0
            f2 = float(f2) if f2 is not None else 0.0
            f3 = float(f3) if f3 is not None else 0.0
            
            fits.append([f1, f2, f3])
        
        if len(fits) == 0:
            result = {
                "paretoFront": [],
                "allIndividuals": [],
                "mooConfig": moo_config,
                "objectives": [moo_config.get("moo.fitness_OA"), moo_config.get("moo.fitness_OB"), moo_config.get("moo.fitness_OC")],
                "select_fea_config_dict": select_fea_config_dict,
                "message": "SwarmData.json 中没有数据"
            }
            return convert_numpy_types(result)
        
        # 计算非支配排序
        if pg is not None:
            try:
                fronts, _, _, _ = pg.fast_non_dominated_sorting(fits)
                rank1_front = fronts[0] if len(fronts) > 0 else []
            except Exception as e:
                # 如果 pygmo 计算失败，使用简单实现
                rank1_front = simple_non_dominated_sorting(fits)
        else:
            # 使用简单的非支配排序实现
            rank1_front = simple_non_dominated_sorting(fits)
        
        # 构建 Pareto 前沿个体列表
        pareto_individuals = []
        for idx in rank1_front:
            if idx < len(swarm_data_items):
                individual_key, individual_data = swarm_data_items[idx]
            else:
                individual_key = f"ind_{idx}"
                individual_data = {}
            
            pareto_individuals.append(convert_numpy_types({
                "index": idx,
                "key": individual_key,
                "project_name": individual_data.get("project_name", ""),
                "individual_index": individual_data.get("individual_index", idx),
                "generation": individual_data.get("number_current_generation", 0),
                "objectives": {
                    "f1": individual_data.get("f1", 0),
                    "f2": individual_data.get("f2", 0),
                    "f3": individual_data.get("f3", 0)
                },
                "performance": {
                    k: v for k, v in individual_data.items() 
                    if k not in ["x_denorm_dict", "project_name", "individual_index", "number_current_generation", "f1", "f2", "f3"]
                },
                "parameters": individual_data.get("x_denorm_dict", {})
            }))
        
        # 构建所有个体列表（用于选择）
        all_individuals = []
        for idx, (key, individual_data) in enumerate(swarm_data_items):
            all_individuals.append(convert_numpy_types({
                "index": idx,
                "key": key,
                "project_name": individual_data.get("project_name", ""),
                "individual_index": individual_data.get("individual_index", idx),
                "generation": individual_data.get("number_current_generation", 0),
                "is_pareto": idx in rank1_front,
                "objectives": {
                    "f1": individual_data.get("f1", 0),
                    "f2": individual_data.get("f2", 0),
                    "f3": individual_data.get("f3", 0)
                },
                "parameters": individual_data.get("x_denorm_dict", {})  # 添加参数数据
            }))
        
        # 按索引排序
        all_individuals.sort(key=lambda x: x["index"])
        
        result = {
            "paretoFront": pareto_individuals,
            "allIndividuals": all_individuals,
            "mooConfig": moo_config,
            "objectives": [
                moo_config.get("moo.fitness_OA"),
                moo_config.get("moo.fitness_OB"),
                moo_config.get("moo.fitness_OC")
            ],
            "select_fea_config_dict": select_fea_config_dict,
            "path2FEACsv": path2FEACsv or machine_designer_data.get("path2FEACsv", ""),
            "path2SwarmData": os.path.dirname(swarm_data_path),
            "folderName": folderName if folderName else os.path.basename(os.path.dirname(swarm_data_path))
        }
        
        # 转换所有 numpy 类型为 Python 原生类型
        return convert_numpy_types(result)
        
    except Exception as e:
        import traceback
        raise HTTPException(
            status_code=500,
            detail=f"获取Pareto前沿数据时出错: {str(e)}\n{traceback.format_exc()}"
        )


@router.get("/generate-geometry-pdf")
async def generate_geometry_pdf(
    folderName: str,
    individualIndex: int
) -> Dict[str, Any]:
    """
    生成指定个体的几何 PDF 文件
    
    Args:
        folderName: 优化文件夹名称（在_default目录下）
        individualIndex: 个体索引（如 3383 对应 'ind3383'）
    
    Returns:
        包含 PDF 文件路径的字典
    """
    try:
        import sys
        codes4_dir = get_codes4_path()
        default_dir = get_default_dir()
        
        # 添加 codes4 目录到 Python 路径
        if codes4_dir not in sys.path:
            sys.path.insert(0, codes4_dir)
        
        # 导入必要的模块
        from machine_design_guide import Modern_Machine_Designer
        
        # 构建文件夹路径
        folder_path = os.path.join(default_dir, folderName)
        if not os.path.exists(folder_path):
            raise HTTPException(
                status_code=404,
                detail=f"文件夹不存在: {folderName}"
            )
        
        # 读取 machine_designer_full.json
        machine_designer_path = os.path.join(folder_path, "machine_designer_full.json")
        if not os.path.exists(machine_designer_path):
            raise HTTPException(
                status_code=404,
                detail=f"machine_designer_full.json 文件不存在: {machine_designer_path}"
            )
        
        # 加载 machine_designer 对象
        mmd = Modern_Machine_Designer.load_from_file_full(machine_designer_path)
        
        # 设置工作目录为文件夹路径（PDF 将生成在这里）
        original_cwd = os.getcwd()
        try:
            os.chdir(folder_path)
            
            # 生成带个体编号的 PDF 文件名
            pdf_filename = f"machine_geometry_ind{individualIndex}.pdf"
            svg_filename = f"machine_geometry_ind{individualIndex}.svg"
            
            # 读取 SwarmData.json 获取个体数据
            swarm_data_path = os.path.join(folder_path, "SwarmData.json")
            if not os.path.exists(swarm_data_path):
                raise HTTPException(
                    status_code=404,
                    detail=f"SwarmData.json 文件不存在: {swarm_data_path}"
                )
            
            import json
            with open(swarm_data_path, "r", encoding="utf-8") as f:
                swarm_data = json.load(f)
            
            # 查找个体键
            suffix = f"ind{individualIndex}"
            target_key = None
            for key in swarm_data.keys():
                if key.endswith(suffix):
                    target_key = key
                    break
            
            if target_key is None:
                raise HTTPException(
                    status_code=404,
                    detail=f"个体 {suffix} 在 SwarmData.json 中未找到"
                )
            
            # 提取并解码 x_denorm_dict
            from modern_machine_designer_utility import Swarm_Data_Analyzer
            x_denorm_dict_raw = swarm_data[target_key]["x_denorm_dict"]
            x_denorm_dict = Swarm_Data_Analyzer.decode_py_reduce_ordered_dict(x_denorm_dict_raw)
            
            # 更新几何参数
            mmd.update_geometric_parameters(x_denorm_dict=x_denorm_dict)
            
            # 手动调用绘图函数，使用自定义文件名
            lw = 0.1 if mmd.mm_r_ro.value < 15 else 0.5
            width_in_points = mmd.mm_r_so.value * 2.1
            height_in_points = mmd.mm_r_so.value * 2.1
            
            # 创建 CairoDrawer 并绘制
            from modern_machine_designer_utility import CairoDrawer
            drawer = CairoDrawer(width_in_points, height_in_points, filename=svg_filename)
            
            # 绘制各个组件
            mmd.machineGeometry['rotorCore'].draw(drawer, bool_draw_whole_model=True)
            mmd.machineGeometry['shaft'].draw(drawer)
            mmd.machineGeometry['rotorMagnet'].draw(drawer, bool_draw_whole_model=True)
            mmd.machineGeometry['statorCore'].draw(drawer, bool_draw_whole_model=True)
            mmd.machineGeometry['coils'].draw(drawer, bool_draw_whole_model=True)
            
            drawer.apply_stroke(lw=lw)
            drawer.surface.finish()
            
            # 转换 SVG 到 PDF（使用自定义文件名）
            import cairosvg
            svg_path = os.path.join(folder_path, svg_filename)
            pdf_path = os.path.join(folder_path, pdf_filename)
            
            if not os.path.exists(svg_path):
                raise HTTPException(
                    status_code=500,
                    detail=f"SVG 文件生成失败: {svg_path}"
                )
            
            # 转换 SVG 到 PDF
            cairosvg.svg2pdf(url=svg_path, write_to=pdf_path)
            
            # 清理 SVG 文件（可选）
            try:
                os.remove(svg_path)
            except:
                pass
            
            if not os.path.exists(pdf_path):
                raise HTTPException(
                    status_code=500,
                    detail=f"PDF 文件生成失败: {pdf_path}"
                )
            
            # 返回相对路径（相对于 backend 目录）
            backend_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
            relative_path = os.path.relpath(pdf_path, backend_dir)
            
            return {
                "pdfPath": relative_path.replace("\\", "/"),  # 统一使用正斜杠
                "absolutePath": pdf_path.replace("\\", "/"),
                "individualIndex": individualIndex,
                "filename": pdf_filename
            }
        finally:
            os.chdir(original_cwd)
            
    except Exception as e:
        import traceback
        raise HTTPException(
            status_code=500,
            detail=f"生成几何 PDF 时出错: {str(e)}\n{traceback.format_exc()}"
        )


class GeometryPdfRequest(BaseModel):
    folderName: str
    individualIndex: int
    parameters: Dict[str, float]

@router.post("/generate-geometry-pdf-with-params")
async def generate_geometry_pdf_with_params(
    request: GeometryPdfRequest
) -> Dict[str, Any]:
    """
    使用自定义参数生成指定个体的几何 PDF 文件
    
    Args:
        request: 包含 folderName, individualIndex 和 parameters 的请求对象
    
    Returns:
        包含 PDF 文件路径的字典
    """
    try:
        import sys
        codes4_dir = get_codes4_path()
        default_dir = get_default_dir()
        
        # 添加 codes4 目录到 Python 路径
        if codes4_dir not in sys.path:
            sys.path.insert(0, codes4_dir)
        
        # 导入必要的模块
        from machine_design_guide import Modern_Machine_Designer
        
        # 从请求中提取参数
        folderName = request.folderName
        individualIndex = request.individualIndex
        parameters = request.parameters
        
        # 构建文件夹路径
        folder_path = os.path.join(default_dir, folderName)
        if not os.path.exists(folder_path):
            raise HTTPException(
                status_code=404,
                detail=f"文件夹不存在: {folderName}"
            )
        
        # 读取 machine_designer_full.json
        machine_designer_path = os.path.join(folder_path, "machine_designer_full.json")
        if not os.path.exists(machine_designer_path):
            raise HTTPException(
                status_code=404,
                detail=f"machine_designer_full.json 文件不存在: {machine_designer_path}"
            )
        
        # 加载 machine_designer 对象
        mmd = Modern_Machine_Designer.load_from_file_full(machine_designer_path)
        
        original_cwd = os.getcwd()
        try:
            os.chdir(folder_path)
            
            # 生成带个体编号和时间戳的 PDF 文件名（避免覆盖）
            import time
            timestamp = int(time.time())
            pdf_filename = f"machine_geometry_ind{individualIndex}_fine_tune_{timestamp}.pdf"
            svg_filename = f"machine_geometry_ind{individualIndex}_fine_tune_{timestamp}.svg"
            
            # 使用自定义参数更新几何参数
            mmd.update_geometric_parameters(x_denorm_dict=parameters)
            
            # 手动调用绘图函数，使用自定义文件名
            lw = 0.1 if mmd.mm_r_ro.value < 15 else 0.5
            width_in_points = mmd.mm_r_so.value * 2.1
            height_in_points = mmd.mm_r_so.value * 2.1
            
            # 创建 CairoDrawer 并绘制
            from modern_machine_designer_utility import CairoDrawer
            drawer = CairoDrawer(width_in_points, height_in_points, filename=svg_filename)
            
            # 绘制各个组件
            mmd.machineGeometry['rotorCore'].draw(drawer, bool_draw_whole_model=True)
            mmd.machineGeometry['shaft'].draw(drawer)
            mmd.machineGeometry['rotorMagnet'].draw(drawer, bool_draw_whole_model=True)
            mmd.machineGeometry['statorCore'].draw(drawer, bool_draw_whole_model=True)
            mmd.machineGeometry['coils'].draw(drawer, bool_draw_whole_model=True)
            
            drawer.apply_stroke(lw=lw)
            drawer.surface.finish()
            
            # 转换 SVG 到 PDF（使用自定义文件名）
            import cairosvg
            svg_path = os.path.join(folder_path, svg_filename)
            pdf_path = os.path.join(folder_path, pdf_filename)
            
            if not os.path.exists(svg_path):
                raise HTTPException(
                    status_code=500,
                    detail=f"SVG 文件生成失败: {svg_path}"
                )
            
            # 转换 SVG 到 PDF
            cairosvg.svg2pdf(url=svg_path, write_to=pdf_path)
            
            # 清理 SVG 文件（可选）
            try:
                os.remove(svg_path)
            except:
                pass
            
            if not os.path.exists(pdf_path):
                raise HTTPException(
                    status_code=500,
                    detail=f"PDF 文件生成失败: {pdf_path}"
                )
            
            # 返回相对路径（相对于 backend 目录）
            backend_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
            relative_path = os.path.relpath(pdf_path, backend_dir)
            
            return {
                "pdfPath": relative_path.replace("\\", "/"),  # 统一使用正斜杠
                "absolutePath": pdf_path.replace("\\", "/"),
                "individualIndex": individualIndex,
                "filename": pdf_filename
            }
        finally:
            os.chdir(original_cwd)
            
    except Exception as e:
        import traceback
        raise HTTPException(
            status_code=500,
            detail=f"生成几何 PDF 时出错: {str(e)}\n{traceback.format_exc()}"
        )


# ==================== 辅助函数 ====================

def convert_numpy_types(obj):
    """
    递归地将 numpy 类型转换为 Python 原生类型
    用于确保 Pydantic 可以正确序列化数据
    """
    try:
        import numpy as np
        
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
    except ImportError:
        # 如果 numpy 未安装，跳过 numpy 类型检查
        pass
    
    if isinstance(obj, dict):
        return {key: convert_numpy_types(value) for key, value in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return [convert_numpy_types(item) for item in obj]
    else:
        return obj

def simple_non_dominated_sorting(fits):
    """
    简单的非支配排序实现（当 pygmo 不可用时使用）
    返回 Rank 1 (Pareto 前沿) 的个体索引列表
    """
    n = len(fits)
    if n == 0:
        return []
    
    # 判断个体 i 是否支配个体 j
    def dominates(i, j):
        # 对于最小化问题，如果 i 的所有目标值都 <= j，且至少有一个 < j，则 i 支配 j
        fit_i = fits[i]
        fit_j = fits[j]
        all_less_equal = all(fit_i[k] <= fit_j[k] for k in range(len(fit_i)))
        at_least_one_less = any(fit_i[k] < fit_j[k] for k in range(len(fit_i)))
        return all_less_equal and at_least_one_less
    
    # 找到所有非支配个体（Rank 1）
    rank1 = []
    for i in range(n):
        is_dominated = False
        for j in range(n):
            if i != j and dominates(j, i):
                is_dominated = True
                break
        if not is_dominated:
            rank1.append(i)
    
    return rank1








