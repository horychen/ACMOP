from fastapi import APIRouter, HTTPException
from typing import List
import os
import json
from app.schemas.project import ProjectSpec

router = APIRouter()

# Helper to load specifications (similar to acmop.py loading)
def load_specifications():
    """Load machine specifications from codes4 directory"""
    # Get the backend directory (parent of app)
    current_dir = os.path.dirname(os.path.abspath(__file__))  # backend/app/routers
    backend_dir = os.path.dirname(os.path.dirname(current_dir))  # backend
    codes4_path = os.path.join(backend_dir, 'codes4')
    spec_path = os.path.join(codes4_path, 'machine_specifications.json')
    
    if not os.path.exists(spec_path):
        # Return empty dict if file doesn't exist
        return {}
        
    with open(spec_path, 'r', encoding='utf-8') as f:
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
