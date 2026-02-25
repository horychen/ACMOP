import sys, os
import asyncio
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'codes4')))

from app.routers.acmop import get_all_parameters
from app.routers.machine_specs import get_machine_specs

async def main():
    print("Testing all-parameters...")
    try:
        res1 = await get_all_parameters()
        print("all-parameters OK.")
    except Exception as e:
        import traceback
        traceback.print_exc()

    print("\nTesting machine-specs...")
    try:
        res2 = await get_machine_specs()
        print("machine-specs OK.")
    except Exception as e:
        import traceback
        traceback.print_exc()

asyncio.run(main())
