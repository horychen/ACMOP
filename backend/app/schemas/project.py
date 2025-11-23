from pydantic import BaseModel
from typing import Dict, Any, Optional

class ProjectSpec(BaseModel):
    # Define fields based on machine_specifications.json structure
    # For now, allowing arbitrary dict as the structure might vary
    pass
