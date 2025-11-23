from fastapi import APIRouter, HTTPException
from typing import Dict, List
import json
import os

router = APIRouter()

def get_codes4_path():
    """Get the path to codes4 directory"""
    # Get the backend directory (parent of app)
    current_dir = os.path.dirname(os.path.abspath(__file__))  # backend/app/routers
    backend_dir = os.path.dirname(os.path.dirname(current_dir))  # backend
    codes4_path = os.path.join(backend_dir, 'codes4')
    
    # Verify the path exists
    if not os.path.exists(codes4_path):
        raise FileNotFoundError(f"codes4 directory not found at: {codes4_path}")
    
    return codes4_path

@router.get("/specs")
async def get_available_specs() -> Dict[str, List[str]]:
    """Get list of available machine specifications"""
    try:
        codes4_path = get_codes4_path()
        specs_file = os.path.join(codes4_path, 'machine_specifications.json')
        
        if not os.path.exists(specs_file):
            raise FileNotFoundError(f"machine_specifications.json not found at: {specs_file}")
        
        with open(specs_file, 'r', encoding='utf-8') as f:
            specs = json.load(f)
        
        return {
            "specs": list(specs.keys())
        }
    except FileNotFoundError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except json.JSONDecodeError as e:
        raise HTTPException(status_code=500, detail=f"Error parsing JSON: {str(e)}")
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error loading specifications: {str(e)}")

@router.get("/fea-configs")
async def get_available_fea_configs() -> Dict[str, List[str]]:
    """Get list of available FEA configuration dictionaries"""
    try:
        codes4_path = get_codes4_path()
        fea_file = os.path.join(codes4_path, 'machine_simulation.json')
        
        if not os.path.exists(fea_file):
            raise FileNotFoundError(f"machine_simulation.json not found at: {fea_file}")
        
        with open(fea_file, 'r', encoding='utf-8') as f:
            fea_configs = json.load(f)
        
        return {
            "fea_configs": list(fea_configs.keys())
        }
    except FileNotFoundError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except json.JSONDecodeError as e:
        raise HTTPException(status_code=500, detail=f"Error parsing JSON: {str(e)}")
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error loading FEA configurations: {str(e)}")

@router.get("/specs/{spec_name}")
async def get_spec_details(spec_name: str) -> Dict:
    """Get detailed information for a specific machine specification"""
    try:
        codes4_path = get_codes4_path()
        specs_file = os.path.join(codes4_path, 'machine_specifications.json')
        
        if not os.path.exists(specs_file):
            raise FileNotFoundError(f"machine_specifications.json not found at: {specs_file}")
        
        with open(specs_file, 'r', encoding='utf-8') as f:
            specs = json.load(f)
        
        if spec_name not in specs:
            raise HTTPException(status_code=404, detail=f"Specification '{spec_name}' not found")
        
        return specs[spec_name]
    except HTTPException:
        raise
    except FileNotFoundError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except json.JSONDecodeError as e:
        raise HTTPException(status_code=500, detail=f"Error parsing JSON: {str(e)}")
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error loading specification: {str(e)}")

@router.get("/fea-configs/{config_name}")
async def get_fea_config_details(config_name: str) -> Dict:
    """Get detailed information for a specific FEA configuration"""
    try:
        codes4_path = get_codes4_path()
        fea_file = os.path.join(codes4_path, 'machine_simulation.json')
        
        if not os.path.exists(fea_file):
            raise FileNotFoundError(f"machine_simulation.json not found at: {fea_file}")
        
        with open(fea_file, 'r', encoding='utf-8') as f:
            fea_configs = json.load(f)
        
        if config_name not in fea_configs:
            raise HTTPException(status_code=404, detail=f"FEA configuration '{config_name}' not found")
        
        return fea_configs[config_name]
    except HTTPException:
        raise
    except FileNotFoundError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except json.JSONDecodeError as e:
        raise HTTPException(status_code=500, detail=f"Error parsing JSON: {str(e)}")
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error loading FEA configuration: {str(e)}")

@router.post("/initialize")
async def initialize_design(config: Dict):
    """Initialize a design process with given configuration"""
    try:
        # This will be implemented to call the backend design process
        # For now, just return success
        return {
            "status": "success",
            "message": "Design initialization received",
            "config": config
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error initializing design: {str(e)}")

@router.get("/visualization/geometric-params")
async def get_geometric_params_for_visualization(spec_name: str) -> Dict:
    """Get geometric parameters from a machine specification for visualization"""
    try:
        codes4_path = get_codes4_path()
        specs_file = os.path.join(codes4_path, 'machine_specifications.json')
        
        if not os.path.exists(specs_file):
            raise FileNotFoundError(f"machine_specifications.json not found at: {specs_file}")
        
        with open(specs_file, 'r', encoding='utf-8') as f:
            specs = json.load(f)
        
        if spec_name not in specs:
            raise HTTPException(status_code=404, detail=f"Specification '{spec_name}' not found")
        
        spec = specs[spec_name]
        
        # Extract geometric parameters
        gp_user = spec.get("GP-user", {})
        
        # Convert to a format suitable for visualization
        geometric_params = {}
        for key, param in gp_user.items():
            geometric_params[key] = {
                "type": param.get("type", "unknown"),
                "value": param.get("value"),
                "bounds": param.get("bounds")
            }
        
        # Also include basic machine topology info
        return {
            "spec_name": spec_name,
            "machine_type": spec.get("machine_type", "Unknown"),
            "slots": spec.get("Qs", 0),
            "poles": spec.get("p", 0),
            "pole_pairs": spec.get("ps", 0),
            "geometric_parameters": geometric_params
        }
    except HTTPException:
        raise
    except FileNotFoundError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except json.JSONDecodeError as e:
        raise HTTPException(status_code=500, detail=f"Error parsing JSON: {str(e)}")
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error loading geometric parameters: {str(e)}")

