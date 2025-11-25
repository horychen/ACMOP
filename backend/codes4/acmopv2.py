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

router = APIRouter(prefix="/api/acmopv2", tags=["ACMOP v2"])


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
        # 获取当前文件所在目录
        current_dir = os.path.dirname(os.path.abspath(__file__))
        json_path = os.path.join(current_dir, "machine_designer.json")
        
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


# ==================== 辅助函数 ====================

def get_codes4_path():
    """获取codes4目录路径"""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    return current_dir


# ==================== 健康检查 ===================

