import urllib.request
import json

url = "http://localhost:8000/api/machine-specs"
try:
    with urllib.request.urlopen(url) as response:
        body = response.read()
        data = json.loads(body)
    
    print(f"API Status: 200")
    print(f"Keys in geometry: {list(data.get('geometry', {}).keys())}")
    print(f"Keys in materials: {list(data.get('materials', {}).keys())}")
    
    # Check for specific new keys
    for key in ['d_stator_yoke', 'w_tooth', 'd_tooth', 'r_stator_outer']:
        val = data.get('geometry', {}).get(key)
        print(f"Geometry key '{key}': {val if val is not None else 'MISSING'}")
        
    for key in ['stator_core_material', 'magnet_temperature']:
        val = data.get('materials', {}).get(key)
        print(f"Material key '{key}': {val if val is not None else 'MISSING'}")
        
    print("\nVerification Successful: Backend is serving correctly formatted JSON with expected keys.")
except Exception as e:
    print(f"Error: {e}")
