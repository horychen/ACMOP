#%% 基本性能计算：电磁负荷、转子外径
from pylab import np, plt
from dataclasses import dataclass

@dataclass
class UserInput_FSPM():
    mec_power: float
    speed_rpm: float
    guess_air_gap_flux_density_B: float = 0.8
    guess_stator_yoke_flux_density_Bsy: float = 0.3
    mm_stack_length: float = 10
    Pa_TangentialStress: float = 6000 # = 0.5*A*B
    m: int = 3
    Qs: int = 24 # 12
    pm: int = 22 # 10

    def __post_init__(self):
        self.Zs = self.Qs * 2
        #     self.required_torque = self.mec_power/(2*np.pi*self.speed_rpm/60)
        #     self.guess_air_gap_flux_density_B = 1.2
        #     guess_linear_current_density_A = self.Pa_TangentialStress*2/self.guess_air_gap_flux_density_B
        #     # (6.2)
        #     self.RotorActiveVolumn_Vr = self.required_torque / (2*self.Pa_TangentialStress)
        #     self.RotorActiveCrossSectArea_Sr = self.RotorActiveVolumn_Vr / (self.mm_stack_length*1e-3)
        #     self.RotorOuterRadius_r_or = np.sqrt(self.RotorActiveCrossSectArea_Sr/np.pi)
        #     self.mm_r_ro = self.RotorOuterRadius_r_or*1e3

        number_of_cells_per_phase_Nc = self.Qs/self.m
        stator_angular_slot_span_theta_s = 2*np.pi / self.Qs
        'The angular difference between stator and rotor teeth should be multiple of 2*pi_elec / (2*self.m)'
        for i in range(1, 2*self.m):
            rotor_angular_slot_span_theta_r_verPositive = stator_angular_slot_span_theta_s/(1+i/2/self.m)
            rotor_angular_slot_span_theta_r_verNegative = stator_angular_slot_span_theta_s/(1-i/2/self.m)
            print(  'pm should be:',
                    2*np.pi/rotor_angular_slot_span_theta_r_verPositive, \
                    2*np.pi/rotor_angular_slot_span_theta_r_verNegative, \
                    number_of_cells_per_phase_Nc*(self.m+0.5*i),\
                    number_of_cells_per_phase_Nc*(self.m-0.5*i) )

# slice = UserInput_FSPM(mec_power=100, speed_rpm=200, Pa_TangentialStress=12000)
# slice = UserInput_FSPM(mec_power=50, speed_rpm=400, Pa_TangentialStress=3000)
slice = UserInput_FSPM(mec_power=30, speed_rpm=400, Pa_TangentialStress=3000)
required_torque = slice.mec_power/(2*np.pi*slice.speed_rpm/60)
guess_linear_current_density_A = slice.Pa_TangentialStress*2/slice.guess_air_gap_flux_density_B

# (6.2)
RotorActiveVolumn_Vr = required_torque / (2*slice.Pa_TangentialStress)
RotorActiveCrossSectArea_Sr = RotorActiveVolumn_Vr / (slice.mm_stack_length*1e-3)
RotorOuterRadius_r_or = np.sqrt(RotorActiveCrossSectArea_Sr/np.pi)
mm_r_ro = RotorOuterRadius_r_or*1e3

aspect_ratio__rotor_axial_to_diameter_ratio = 2*mm_r_ro/slice.mm_stack_length


#% 磁路法计算：永磁体大小？

mm_air_gap_length = 3
mm_r_si     = mm_r_ro + mm_air_gap_length
mm_r_airgap = mm_r_ro + mm_air_gap_length*0.5

rotor_to_stator_tooth_width_ratio = 1.4
fea_config_dict = { "WindingFill": 0.3,
                    "Js": 4e6, # A/m^2
                    }

mm_stator_tooth_width_w_UCoreWidth = 2*np.pi*mm_r_si / slice.Qs / 4
mm_rotor_tooth_width_w_rt = rotor_to_stator_tooth_width_ratio * mm_stator_tooth_width_w_UCoreWidth # Empirical: slightly wider (5.18) Habil-Gruber

mmf_per_phase = guess_linear_current_density_A * 2*np.pi*mm_r_airgap / slice.m
mmf_per_slot  = guess_linear_current_density_A * 2*np.pi*mm_r_airgap / slice.Qs
# stator_current_per_phase = mmf_per_phase /  2 * N * I

stator_slot_area = mmf_per_slot / fea_config_dict['Js'] / fea_config_dict['WindingFill'] 

# stator slot depth
stator_inner_radius_r_is_eff = mm_r_si
temp = (2*np.pi*stator_inner_radius_r_is_eff - slice.Qs*mm_stator_tooth_width_w_UCoreWidth*3) # mm_stator_tooth_width_w_UCoreWidth * 3 = stator non-slot part width
stator_tooth_depth_d_st = ( np.sqrt(temp**2 + 4*np.pi*stator_slot_area*slice.Qs) - temp ) / (2*np.pi)
print(stator_tooth_depth_d_st)

mm_d_sy = slice.guess_air_gap_flux_density_B * mm_stator_tooth_width_w_UCoreWidth / slice.guess_stator_yoke_flux_density_Bsy


# mm_r_so = mm_r_si + mm_d_st + mm_d_sy
# split_ratio = mm_r_ro / mm_r_so

# %%

# %%
