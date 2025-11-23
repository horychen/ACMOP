from fastapi import APIRouter, HTTPException
from typing import List
import os
import json
from app.schemas.project import ProjectSpec

router = APIRouter()

# Helper to load specifications (similar to acmop.py loading)
def load_specifications():
    # Assuming machine_specifications.json is in backend/codes4/
    # Adjust path as necessary based on actual location
    base_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', 'backend', 'codes4'))
    spec_path = os.path.join(base_path, 'machine_specifications.json')
    
    if not os.path.exists(spec_path):
        # Fallback or error
        return {}
        
    with open(spec_path, 'r') as f:
        return json.load(f)

@router.get("/", response_model=List[str])
async def list_projects():
    specs = load_specifications()
    return list(specs.keys())

@router.get("/{project_name}")
async def get_project_details(project_name: str):
    specs = load_specifications()
    if project_name not in specs:
        raise HTTPException(status_code=404, detail="Project not found")
    return specs[project_name]
