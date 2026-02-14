import requests

try:
    r = requests.get('http://127.0.0.1:8000/api/machine-specs')
    r.raise_for_status()
    data = r.json()
    geom = data.get('geometry', {})
    winding = data.get('winding', {})
    
    print(f"r_shaft in geom keys: {'r_shaft' in geom}")
    if 'r_shaft' in geom:
        print(f"r_shaft: {geom['r_shaft']}")
        
    w_keys = ['slot_count', 'pole_count', 'wire_diameter_with_insulation']
    for k in w_keys:
        print(f"winding.{k}: {winding.get(k)}")
            
except Exception as e:
    print(f"Error: {e}")
