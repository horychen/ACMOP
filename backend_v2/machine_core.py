from pydantic import BaseModel
from typing import Optional, List

class StatorParameters(BaseModel):
    OD: float = 13.0
    ID: float = 8.3
    toothDepth: float = 2.0
    toothWidth: float = 1.2
    tooth_shape: str = "closed"
    liner: float = 0.1

class RotorParameters(BaseModel):
    OD: float = 8.0
    ID: float = 3.0
    airGap: float = 0.15
    magnetDepth: float = 1.0

class WindingParameters(BaseModel):
    slot_count: int = 12
    pole_count: int = 10
    awg: int = 31
    l_stack: float = 50.0
    rated_speed: float = 3000.0
    rated_current_density: float = 5.0

class MotorParameters(BaseModel):
    stator: StatorParameters = StatorParameters()
    rotor: RotorParameters = RotorParameters()
    winding: WindingParameters = WindingParameters()
