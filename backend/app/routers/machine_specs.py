from fastapi import APIRouter, HTTPException
from dataclasses import asdict
from backend.codes4.user_minitureMachine import MotorSpecs

router = APIRouter()

@router.get("/machine-specs")
async def get_machine_specs():
    try:
        specs = MotorSpecs()
        return asdict(specs)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
