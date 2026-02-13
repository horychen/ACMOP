import sys
import os
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

# Add project root and codes4 to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'backend', 'codes4')))

from app.routers import project, optimization, results, design, acmop, machine_specs, debug

app = FastAPI(title="ACMOP Backend", version="1.0.0")

# CORS configuration
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000"],  # Next.js default port
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Include routers
app.include_router(project.router, prefix="/api/projects", tags=["Projects"])
app.include_router(optimization.router, prefix="/api/optimization", tags=["Optimization"])
app.include_router(results.router, prefix="/api/results", tags=["Results"])
app.include_router(design.router, prefix="/api/design", tags=["Design"])
app.include_router(acmop.router)  # ACMOP v2 API
app.include_router(machine_specs.router, prefix="/api", tags=["Machine Specs"])
app.include_router(debug.router, prefix="/api/debug", tags=["Debug"])

@app.get("/api/health")
async def health():
    return {"status": "ok", "time": os.path.getmtime(__file__)}

@app.get("/")
async def root():
    return {"message": "Welcome to ACMOP Backend API"}
