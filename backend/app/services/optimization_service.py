import sys
import os
import logging

# Ensure codes4 is in path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', 'backend', 'codes4')))

try:
    import acmop
except ImportError:
    # Handle case where imports fail (e.g. during testing without full environment)
    acmop = None

class OptimizationService:
    def __init__(self):
        self.logger = logging.getLogger(__name__)

    def run_optimization(self, spec_name: str, fea_config: str, project_loc: str):
        if not acmop:
            raise ImportError("acmop module not found")
        
        # Initialize wrapper
        mop = acmop.AC_Machine_Optiomization_Wrapper(
            select_spec=spec_name,
            select_fea_config_dict=fea_config,
            project_loc=project_loc,
            bool_show_GUI=False 
        )
        
        # Run optimization (Part 4)
        mop.part_optimization()
        return {"status": "Optimization started"}

    def run_winding_part(self, spec_name: str, fea_config: str, project_loc: str):
        if not acmop:
            raise ImportError("acmop module not found")
            
        mop = acmop.AC_Machine_Optiomization_Wrapper(
            select_spec=spec_name,
            select_fea_config_dict=fea_config,
            project_loc=project_loc,
            bool_show_GUI=False
        )
        
        mop.part_winding()
        return {"status": "Winding part executed"}

optimization_service = OptimizationService()
