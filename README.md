# ACMOP

> Alternating Current Machine Optimization Project

## Requirements (Obsolete):

- JMAG Designer 17.1.05i
> Higher version might not work as the invoking method names might be changed.

- Anaconda3-2020.02-Windows-x86_64.exe with python 3.7.6
> If you use newer version of anaconda3 after 2020.02, pacakge pygmo might not work, as far as I know.

- pygmo:
    - conda config --add channels conda-forge
    - conda install pygmo

- pyx, pyfemm, streamlit, jsonpickle and others (if any) can be installed via pip

- texlive, or other latex compiler (if you want pdf report for motor design)

## Installation:

- 1. Get [Anaconda3-2021.05-Windows-x86_64.exe 477.2M 2021-05-13 22:08:48](https://repo.anaconda.com/archive/)

- 2. (OPTIONAL) Create a virtual environment by `conda create -n your-env-name python=3.8.8 pygmo`, and then activate your virtual env by `conda activate your-env-name`.

- 3. Install the rest dependencies from PyPI: `pip install pyx, pyfemm, jsonpickle, recordtype, pint`

## Features
- Restartable. Upon interrupts, this program is able to re-start from the file "swarm_data.txt" of the current run (or even from a different run).

- Obtian the machine_design_variant python object after the optimization.
```
# Recover a design from its jsonpickle-object file
# cd to where __p2ps1-Q12y3-0001.json is located.
import sys; sys.path.insert(0, r'D:\DrH\acmop\codes3')
import utility_json; variant = utility_json.from_json_recursively('p2ps1-Q12y3-0001')
variant.build_jmag_project(variant.project_meta_data)
```

## TODO (FEMM)
- There is no need to draw coils smaller than slots in FEMM and this has been causing too many meshes inside the slots, slowing down FEA in FEMM.
- The sliding band is drawing inside the region of mechanical air gap. Since the sleeve band is not even modeled in FEMM, maybe we should draw sliding band at the middle of the magnetic air gap. I am not sure if this will improve the meshing quality in the air gap.
- Magnet eddy current is not yet considered.
- Run multiple FEMM instances in parallel to sovle for transient FEA results.

## TODO (JMAG)

- Currently, DPNV winding is implemented as a dual three phase system excited by current source inverter. Therefore, the JMAG circuit's (current source) excitation will give wrong terminal voltage results and the no voltage property is not observed (but the electromagnetic performance is correct because the current is correct).
- Variable d_rp and d_rs can be larger than d_pm. This is not allowed for using sleeve to retain the magnet. Should add new variable like d_pm_to_iron instead.
- Upon restarting, the constraints are not applied to the archive. High ripple design or low FRW may take up space in Pareto front.
- The constraints are fake constraints in current optimization code implementation. For example, if the optimization is re-started from the archive, those designs with FRW<0.75 could dominate other designs even though the constraint should be FRW>0.75.
- The iron loss data is only valid from 50 Hz up to 600 Hz. See http://www.femm.info/wiki/SPMLoss (Search for 530).
- Make the 2.5D plot a function for group member to use (slack writing tools)
- Revise build_str_results function to include rotor weight after "Average Force Mag".
- Ashad: Following Pyrhonen's book to derive a BPMSM initial design (we use Bianchi 2006 for now--which is simplified) for your 160,000 rpm BPMSM.
- Take care of the dummy variables in the initial design script for IM and PMSM. (At least let the user know they don't need to care about those variables.)
- New re-starting feature: use the optimal design as the new template design, and use a new bound according to the new template design.
- Ashad: I think the correct way of calculating error angle and magnitude for Q12p4 is by taking 5 electrical periods as suspension winding has upto 1/5 subharmonics. p^aster=5? See harmonics calculation.pdf. Jiahao: I am adding this comment from Ashad to Readme.md file of bopt-python for future development reference---the simulation setting is dependent on winding layout.
- Flux weakening and strenthening have great impact on the suspension force performance/capability. This needs more research.
- There is still bug in displacement power factor calculation:
> Different Q6 p1 ps2 designs tested for transient decay phenomenon in 2nd TSS.
> A: 3 cycles + 0.5 cycle
> B: 3 cycles + 0.5 cycle + 10 cycles
>                                         [kNm/m3]  [1]    [%]     [%]     [deg]  [%]     [$]     [1]
>                                           TRV     FRW    Trip    Em      Ea     eta     Cost    Power factor
> proj3297-SPMSM_IDQ6p1s1 | PMSM Q06p1y2 A 29.736, 0.926, 20.542, 10.372, 5.466, 93.806, 199.509, 0.947
> proj3297-SPMSM_IDQ6p1s1 | PMSM Q06p1y2 B 29.790, 0.928, 20.549, 13.810, 7.654, 93.815, 199.374, 0.947
> 
> proj3303-SPMSM_IDQ6p1s1 | PMSM Q06p1y2 A 28.749, 0.891, 23.464, 9.416, 4.963, 93.653, 200.987, 0.946
> proj3303-SPMSM_IDQ6p1s1 | PMSM Q06p1y2 B 28.796, 0.892, 23.549, 12.750, 7.323, 93.662, 200.862, 0.947
> 
> proj3376-SPMSM_IDQ6p1s1 | PMSM Q06p1y2 A 29.588, 0.943, 23.687, 12.583, 6.161, 93.767, 199.832, 0.945
> proj3376-SPMSM_IDQ6p1s1 | PMSM Q06p1y2 B 29.628, 0.945, 23.724, 16.675, 9.147, 93.773, 199.733, 0.945
> 
> proj3527-SPMSM_IDQ6p1s1 | PMSM Q06p1y2 A 29.588, 0.943, 23.687, 12.583, 6.161, 93.767, 199.832, 0.945
> proj3527-SPMSM_IDQ6p1s1 | PMSM Q06p1y2 B 29.628, 0.945, 23.724, 16.675, 9.147, 93.773, 199.733, 0.945
> 
> proj3534-SPMSM_IDQ6p1s1 | PMSM Q06p1y2 A 29.787, 0.933, 20.818, 10.372, 5.524, 93.766, 198.617, 0.947
> proj3534-SPMSM_IDQ6p1s1 | PMSM Q06p1y2 B 29.842, 0.935, 20.839, 13.712, 7.580, 93.775, 198.482, 0.947
> 
> proj3548-SPMSM_IDQ6p1s1 | PMSM Q06p1y2 A 28.703, 0.915, 17.838, 13.290, 6.590, 93.914, 193.032, 0.951
> proj3548-SPMSM_IDQ6p1s1 | PMSM Q06p1y2 B 28.752, 0.917, 17.863, 17.776, 9.896, 93.922, 192.902, 0.951
> 
> proj3747-SPMSM_IDQ6p1s1 | PMSM Q06p1y2 A 29.196, 0.950, 21.184, 10.667, 6.124, 93.467, 202.410, 0.933
> proj3747-SPMSM_IDQ6p1s1 | PMSM Q06p1y2 B 29.244, 0.951, 21.236, 13.908, 6.860, 93.476, 202.290, 0.933
> 
> proj3748-SPMSM_IDQ6p1s1 | PMSM Q06p1y2 A 29.277, 0.917, 21.198, 10.592, 5.719, 93.760, 197.828, -0.657
> proj3748-SPMSM_IDQ6p1s1 | PMSM Q06p1y2 B 29.329, 0.919, 21.252, 13.779, 8.086, 93.769, 197.693, 0.949
> 
> proj3850-SPMSM_IDQ6p1s1 | PMSM Q06p1y2 A 30.820, 1.019, 26.083, 11.041, 6.136, 93.786, 192.525, 0.945
> proj3850-SPMSM_IDQ6p1s1 | PMSM Q06p1y2 B 30.867, 1.021, 26.070, 14.756, 8.258, 93.794, 192.413, 0.945
> 
> proj3895-SPMSM_IDQ6p1s1 | PMSM Q06p1y2 A 29.329, 0.953, 21.010, 10.600, 6.140, 93.487, 202.079, 0.932
> proj3895-SPMSM_IDQ6p1s1 | PMSM Q06p1y2 B 29.379, 0.955, 21.111, 13.794, 6.846, 93.496, 201.957, 0.932

- More raw data should be added to json file. For example, the cost of steel, cost of PM and cost of copper.
- There is a bug when building x_denorm of the initial design of induction motor, where the stator yoke depth is negative (the original stator outer radius is smaller than the computed one).

- [TODO] Iron loss must use half cycle FFT as one fourth cycle could be wrong.

- [TODO] Restarting optimziation with constaints applied.

- [TODO] Double check winding's initial excitation position.
