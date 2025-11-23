from fastapi import APIRouter, BackgroundTasks, HTTPException
from pydantic import BaseModel
from app.services.optimization_service import optimization_service

router = APIRouter()

class OptimizationRequest(BaseModel):
    spec_name: str
    fea_config: str
    project_loc: str

@router.post("/start")
async def start_optimization(request: OptimizationRequest, background_tasks: BackgroundTasks):
    try:
        # Run in background as it's a long running process
        background_tasks.add_task(
            optimization_service.run_optimization, 
            request.spec_name, 
            request.fea_config, 
            request.project_loc
        )
        return {"message": "Optimization started in background"}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router.post("/winding")
async def run_winding(request: OptimizationRequest):
    try:
        result = optimization_service.run_winding_part(
            request.spec_name, 
            request.fea_config, 
            request.project_loc
        )
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
