import sys, os
import asyncio
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'codes4')))
from app.routers.debug import get_debug_points, get_debug_geometry, get_inspection_data

async def main():
    try:
        print("Testing points...")
        res1 = await get_debug_points()
        print("Points OK. Testing geometry...")
        res2 = await get_debug_geometry()
        print("Geometry OK. Testing inspection...")
        res3 = await get_inspection_data()
        print("Inspection OK.")
    except Exception as e:
        import traceback
        traceback.print_exc()

asyncio.run(main())
