from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from machine_core import MotorParameters

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

@app.get("/api/motor-parameters")
def get_motor_parameters():
    return motor_params

@app.post("/api/motor-parameters")
def update_motor_parameters(params: MotorParameters):
    global motor_params
    motor_params = params
    return motor_params

if __name__ == "__main__":
    import uvicorn
    # Make sure this runs on a different port than the frontend or V1 backend
    uvicorn.run(app, host="0.0.0.0", port=8000)
