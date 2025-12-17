from fastapi import APIRouter, HTTPException
import os
from fastapi import APIRouter, HTTPException
import os
import json
from typing import List, Dict, Any

router = APIRouter()

# Helper to find swarm data
@router.get("/swarm/{project_name}")
async def get_swarm_data(project_name: str, project_loc: str = "../_default/"):
    # Construct path: backend/_default/{project_name}/{project_name}.json
    # We need to be careful about relative paths. 
    # Assuming project_loc is relative to backend/codes4/ or backend/app/
    # Let's use absolute path based on where we know the file is.
    
    # Base path for _default directory (assuming it's at backend/_default)
    base_default_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '_default'))
    
    # Try finding the file
    # Pattern 1: {project_name}/{project_name}.json
    path1 = os.path.join(base_default_dir, project_name.replace(' ', '_'), f"{project_name}.json")
    # Pattern 2: {project_name}/swarm_data.json
    path2 = os.path.join(base_default_dir, project_name.replace(' ', '_'), "swarm_data.json")
    
    target_path = None
    if os.path.exists(path1):
        target_path = path1
    elif os.path.exists(path2):
        target_path = path2
    else:
        # Fallback to checking if project_loc was provided differently (not implemented yet)
        pass

    if not target_path or not os.path.exists(target_path):
         raise HTTPException(status_code=404, detail=f"Swarm data not found for {project_name}")

    with open(target_path, 'r') as f:
        content = f.read()
        
    # Handle the weird format starting with comma
    if content.strip().startswith(','):
        # The file might start with whitespace then a comma, or just a comma.
        # User's script uses buf[1:], assuming the comma is at index 0.
        # To be robust, let's find the first comma.
        first_comma_index = content.find(',')
        if first_comma_index != -1:
            try:
                # content[first_comma_index+1:] skips the comma
                data = json.loads('{' + content[first_comma_index+1:] + '}')
            except json.JSONDecodeError:
                 # Try normal load if that fails
                try:
                    data = json.loads(content)
                except:
                    raise HTTPException(status_code=500, detail="Failed to parse swarm data JSON")
        else:
             # Should not happen if strip().startswith(',') is true
             raise HTTPException(status_code=500, detail="Invalid JSON format detected")
    else:
        data = json.loads(content)
        
    return data

@router.get("/csv/list/{project_name}")
async def list_csv_files(project_name: str):
    """List all CSV files in the project's csv directory."""
    base_default_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '_default'))
    # Assuming the structure is backend/_default/{project_name}/csv
    # We need to handle the project name carefully. The folder might be named with underscores instead of spaces.
    project_dir_name = project_name.replace(' ', '_')
    csv_dir = os.path.join(base_default_dir, project_dir_name, "csv")
    
    if not os.path.exists(csv_dir):
        # Try finding the project directory first if the exact name match fails
        if not os.path.exists(os.path.join(base_default_dir, project_dir_name)):
             raise HTTPException(status_code=404, detail=f"Project directory not found for {project_name}")
        return [] # Return empty list if csv dir doesn't exist but project does

    csv_files = [f for f in os.listdir(csv_dir) if f.endswith('.csv')]
    return csv_files

@router.get("/csv/content/{project_name}/{filename}")
async def get_csv_content(project_name: str, filename: str):
    """Get the content of a specific CSV file."""
    base_default_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '_default'))
    project_dir_name = project_name.replace(' ', '_')
    csv_path = os.path.join(base_default_dir, project_dir_name, "csv", filename)
    
    if not os.path.exists(csv_path):
        raise HTTPException(status_code=404, detail=f"CSV file not found: {filename}")
        
    with open(csv_path, 'r') as f:
        content = f.read()
        
    return {"content": content}

@router.get("/csv/list-from-path")
async def list_csv_files_from_path(path: str):
    """List all CSV files from a given path (path2FEACsv)."""
    if not path:
        raise HTTPException(status_code=400, detail="Path parameter is required")
    
    # Normalize the path
    csv_dir = os.path.normpath(path)
    
    if not os.path.exists(csv_dir):
        raise HTTPException(status_code=404, detail=f"CSV directory not found: {csv_dir}")
    
    if not os.path.isdir(csv_dir):
        raise HTTPException(status_code=400, detail=f"Path is not a directory: {csv_dir}")
    
    try:
        csv_files = [f for f in os.listdir(csv_dir) if f.endswith('.csv')]
        return csv_files
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error reading directory: {str(e)}")

@router.get("/csv/content-from-path")
async def get_csv_content_from_path(path: str, filename: str):
    """Get the content of a specific CSV file from a given path."""
    if not path or not filename:
        raise HTTPException(status_code=400, detail="Path and filename parameters are required")
    
    # Normalize the path and join with filename
    csv_dir = os.path.normpath(path)
    csv_path = os.path.join(csv_dir, filename)
    
    # Security check: ensure the file is within the csv_dir
    csv_path = os.path.normpath(csv_path)
    if not csv_path.startswith(os.path.normpath(csv_dir)):
        raise HTTPException(status_code=403, detail="Invalid file path")
    
    if not os.path.exists(csv_path):
        raise HTTPException(status_code=404, detail=f"CSV file not found: {filename}")
    
    try:
        with open(csv_path, 'r', encoding='utf-8') as f:
            content = f.read()
        return {"content": content}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error reading file: {str(e)}")

@router.get("/pdf")
async def get_pdf(path: str):
    """Get a PDF file by path (relative to backend directory)."""
    from fastapi.responses import Response
    
    import os
    current_dir = os.path.dirname(os.path.abspath(__file__))
    backend_dir = os.path.abspath(os.path.join(current_dir, '..', '..'))
    
    # 构建完整路径（相对于 backend 目录）
    pdf_path = os.path.join(backend_dir, path)
    pdf_path = os.path.normpath(os.path.abspath(pdf_path))
    
    # 安全检查：确保路径在 backend 目录下
    backend_dir_normalized = os.path.normpath(os.path.abspath(backend_dir))
    if not pdf_path.startswith(backend_dir_normalized):
        raise HTTPException(
            status_code=403,
            detail="Access denied: Path outside backend directory"
        )
    
    if not os.path.exists(pdf_path):
        raise HTTPException(
            status_code=404,
            detail=f"PDF file not found: {pdf_path}"
        )
    
    if not pdf_path.lower().endswith('.pdf'):
        raise HTTPException(
            status_code=400,
            detail="File is not a PDF"
        )
    
    # 读取 PDF 文件内容
    try:
        with open(pdf_path, 'rb') as f:
            content = f.read()
        
        if not content:
            raise HTTPException(status_code=500, detail="PDF file is empty")
        
        # 使用 Response 而不是 FileResponse，并设置 Content-Disposition 为 inline
        # 这样浏览器会内联显示 PDF 而不是下载
        return Response(
            content=content,
            media_type='application/pdf',
            headers={
                'Content-Disposition': f'inline; filename="{os.path.basename(pdf_path)}"',
                'Content-Type': 'application/pdf',
                'Cache-Control': 'no-cache'
            }
        )
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error reading PDF file: {str(e)}")


@router.get("/pdf/machine-geometry")
async def get_machine_geometry_pdf():
    """Get the machine_geometry.pdf file."""
    from fastapi.responses import Response
    
    # Try to find the PDF in codes4 directory
    import os
    current_dir = os.path.dirname(os.path.abspath(__file__))
    # results.py is in backend/app/routers/, so go up to backend/, then to codes4
    backend_dir = os.path.abspath(os.path.join(current_dir, '..', '..'))
    codes4_dir = os.path.join(backend_dir, 'codes4')
    pdf_path = os.path.join(codes4_dir, 'machine_geometry.pdf')
    
    # Normalize the path
    pdf_path = os.path.normpath(os.path.abspath(pdf_path))
    
    # Debug: log the path being checked
    if not os.path.exists(pdf_path):
        # Try alternative path: maybe codes4 is at the same level as app
        alt_codes4_dir = os.path.join(os.path.dirname(backend_dir), 'codes4')
        alt_pdf_path = os.path.join(alt_codes4_dir, 'machine_geometry.pdf')
        alt_pdf_path = os.path.normpath(os.path.abspath(alt_pdf_path))
        
        if os.path.exists(alt_pdf_path):
            pdf_path = alt_pdf_path
        else:
            raise HTTPException(
                status_code=404, 
                detail=f"PDF file not found. Checked: {pdf_path}, {alt_pdf_path}"
            )
    
    try:
        # Read the PDF file synchronously
        with open(pdf_path, 'rb') as f:
            content = f.read()
        
        if not content:
            raise HTTPException(status_code=500, detail="PDF file is empty")
        
        return Response(
            content=content,
            media_type='application/pdf',
            headers={
                'Content-Disposition': 'inline; filename="machine_geometry.pdf"',
                'Content-Type': 'application/pdf',
                'Cache-Control': 'no-cache'
            }
        )
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error reading PDF file: {str(e)}")