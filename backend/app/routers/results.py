from fastapi import APIRouter, HTTPException
import os
from fastapi import APIRouter, HTTPException
import os
import json
from typing import List, Dict, Any

router = APIRouter()

# Helper to find swarm data
@router.get("/swarm/{project_name}")
async def get_swarm_data(project_name: str, project_loc: str = "../_default/"):
    # Construct path: backend/_default/{project_name}/{project_name}.json
    # We need to be careful about relative paths. 
    # Assuming project_loc is relative to backend/codes4/ or backend/app/
    # Let's use absolute path based on where we know the file is.
    
    # Base path for _default directory (assuming it's at backend/_default)
    base_default_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '_default'))
    
    # Try finding the file
    # Pattern 1: {project_name}/{project_name}.json
    path1 = os.path.join(base_default_dir, project_name.replace(' ', '_'), f"{project_name}.json")
    # Pattern 2: {project_name}/swarm_data.json
    path2 = os.path.join(base_default_dir, project_name.replace(' ', '_'), "swarm_data.json")
    
    target_path = None
    if os.path.exists(path1):
        target_path = path1
    elif os.path.exists(path2):
        target_path = path2
    else:
        # Fallback to checking if project_loc was provided differently (not implemented yet)
        pass

    if not target_path or not os.path.exists(target_path):
         raise HTTPException(status_code=404, detail=f"Swarm data not found for {project_name}")

    with open(target_path, 'r') as f:
        content = f.read()
        
    # Handle the weird format starting with comma
    if content.strip().startswith(','):
        # The file might start with whitespace then a comma, or just a comma.
        # User's script uses buf[1:], assuming the comma is at index 0.
        # To be robust, let's find the first comma.
        first_comma_index = content.find(',')
        if first_comma_index != -1:
            try:
                # content[first_comma_index+1:] skips the comma
                data = json.loads('{' + content[first_comma_index+1:] + '}')
            except json.JSONDecodeError:
                 # Try normal load if that fails
                try:
                    data = json.loads(content)
                except:
                    raise HTTPException(status_code=500, detail="Failed to parse swarm data JSON")
        else:
             # Should not happen if strip().startswith(',') is true
             raise HTTPException(status_code=500, detail="Invalid JSON format detected")
    else:
        data = json.loads(content)
        
    return data
