from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

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

app = FastAPI(title="ACMOP Backend V2")

# Allow CORS for Next.js frontend
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000", "http://localhost:3001"], # Add Next.js default ports
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Global in-memory state for demonstration, later this should be tied to a project/session
motor_params = MotorParameters()
# Initialize with some default test values mimicking user_input
motor_params.stator.OD = 13.0
motor_params.stator.ID = 8.3
motor_params.stator.toothDepth = 2.0
motor_params.stator.toothWidth = 1.2
motor_params.stator.tooth_shape = "closed"

from machine import user_input

@app.get("/api/motor-parameters")
def get_motor_parameters():
    geo = user_input.get('geometry', {})
    win = user_input.get('winding', {})
    
    r_so = geo.get('r_stator_outer', 6.5)
    r_ro = geo.get('r_rotor_outer', 4.0)
    d_gap = geo.get('d_air_gap', 0.15)
    d_tooth = geo.get('d_tooth', 2.0)
    w_tooth = geo.get('w_tooth', 1.2)
    d_mag = geo.get('d_magnet', 3.0)
    r_shaft = geo.get('r_shaft', 0.0)
    
    r_si = r_ro + d_gap
    yoke = r_so - r_si - d_tooth
    
    params = {
        "stator": {
            "OD": r_so * 2,
            "ID": r_si * 2,
            "toothDepth": d_tooth,
            "toothWidth": w_tooth,
            "yoke": yoke,
            "liner": 0.1
        },
        "rotor": {
            "OD": r_ro * 2,
            "ID": r_shaft * 2,
            "magnetDepth": d_mag,
            "airGap": d_gap
        },
        "winding": {
            "awg": 31,
            "J": win.get('rated_current_density', 14)
        }
    }
    return params

@app.post("/api/motor-parameters")
def update_motor_parameters(params: MotorParameters):
    # This was originally updating a local dict, for now let's just mock it
    # since actual updating requires modifying user_input
    return params

from fastapi.responses import FileResponse
import os

@app.get("/api/stator-svg")
def get_stator_svg():
    file_path = os.path.join(os.path.dirname(__file__), "stator_v2.svg")
    if os.path.exists(file_path):
        return FileResponse(file_path, media_type="image/svg+xml")
    return {"error": "SVG file not found"}

if __name__ == "__main__":
    import uvicorn
    # Make sure this runs on a different port than the frontend or V1 backend
    uvicorn.run(app, host="0.0.0.0", port=8001)
