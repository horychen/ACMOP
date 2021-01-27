import os
# obsolete
def where_am_i(fea_config_dict):
    dir_interpreter = os.path.abspath('')
    print(dir_interpreter)
    if dir_interpreter[0] == 'D':
        print('you are on Legion Y730')
        dir_parent = 'D:/OneDrive - UW-Madison/c/'
        dir_codes = dir_parent + 'codes3/'
        dir_lib = dir_parent + 'codes3/'
        dir_femm_files = 'D:/femm42/' # .ans files are too large to store on OneDrive anymore
        dir_project_files = 'D:/JMAG_Files/'
        pc_name = 'Y730'
    elif dir_interpreter[0] == 'I':
        print('you are on Severson02')
        dir_parent = 'I:/jchen782/c/'
        dir_codes = dir_parent + 'codes3/'
        dir_lib = dir_parent + 'codes3/'
        dir_femm_files = 'I:/jchen782/FEMM/'
        dir_project_files = 'I:/jchen782/JMAG/'
        pc_name = 'Severson02'
    elif dir_interpreter[0] == 'K':
        print('you are on Severson01')
        dir_parent = 'K:/jchen782/c/'
        dir_codes = dir_parent + 'codes3/'
        dir_lib = dir_parent + 'codes3/'
        dir_femm_files = 'K:/jchen782/FEMM/'
        dir_project_files = 'K:/jchen782/JMAG/'
        pc_name = 'Severson01'
    else:
        # elif 'Chen' in 'dir_interpreter':
        print('you are on T440p')
        # dir_parent = 'C:/Users/Hory Chen/OneDrive - UW-Madison/c/'
        dir_parent = 'd:/c/'
        dir_codes = dir_parent + 'codes3/'
        dir_lib = dir_parent + 'codes3/'
        dir_femm_files = 'd:/femm42/'
        dir_project_files = 'd:/JMAG_Files/'
        pc_name = 'T440p'

    os.chdir(dir_codes)

    fea_config_dict['dir_parent']            = dir_parent
    fea_config_dict['dir_lib']               = dir_lib
    fea_config_dict['dir_codes']             = dir_codes
    fea_config_dict['dir_femm_files']        = dir_femm_files
    fea_config_dict['dir_project_files']     = dir_project_files
    # fea_config_dict['dir_initial_design']    = dir_initial_design
    # fea_config_dict['dir_csv_output_folder'] = dir_csv_output_folder
    fea_config_dict['pc_name']               = pc_name
    fea_config_dict['dir_interpreter']       = dir_interpreter

    if pc_name == 'Y730':
        fea_config_dict['delete_results_after_calculation'] = False # save disk space
        if fea_config_dict['Restart'] == False:
            fea_config_dict['OnlyTableResults'] = True  # save disk space for my PC
    # However, we need field data for iron loss calculation
    fea_config_dict['OnlyTableResults'] = False 

# self-introspective
def where_am_i_v2(fea_config_dict, bool_post_processing=False):
    def get_pc_name():
        import platform
        import socket
        n1 = platform.node()
        n2 = socket.gethostname()
        n3 = os.environ["COMPUTERNAME"]
        if n1 == n2 == n3:
            return n1
        elif n1 == n2:
            return n1
        elif n1 == n3:
            return n1
        elif n2 == n3:
            return n2
        else:
            raise Exception("Computer names are not equal to each other.")

    dir_interpreter = os.path.abspath('') + '/'
    if 'DrH' not in dir_interpreter:
        # print(dir_interpreter)
        raise Exception("Revert this back before you run on server.")
    # 先用这个救救急，到服务器的时候这里要报错
    # dir_interpreter = r'D:\DrH\bopt-python\codes3/'

    dir_parent = dir_interpreter #+ '../' # '../' <- relative path is not always a good option

    # if not bool_post_processing:
    #     if 'codes3' not in dir_interpreter:
    #         raise Exception('Do not run the script from your current directory: %s.\nPlease cd to ./codes3 before running python scripts.'%(dir_interpreter))
    dir_codes = dir_parent + 'codes3/'
    # dir_femm_files = dir_parent + 'femm_files/'
    pc_name = get_pc_name()
    os.chdir(dir_codes)
    # print(dir_parent, dir_codes, pc_name, sep='\n'); quit()
    fea_config_dict['dir.interpreter']       = dir_interpreter 
    fea_config_dict['dir.parent']            = dir_parent
    fea_config_dict['dir.codes']             = dir_codes
    fea_config_dict['dir.femm_files']        = dir_parent
    fea_config_dict['pc_name']               = pc_name

    fea_config_dict['delete_results_after_calculation'] = False # True for saving disk space (but you lose voltage profile and element data)
    if fea_config_dict['designer.Restart'] == False:
        fea_config_dict['designer.OnlyTableResults'] = True  # save disk space for my PC
    # However, we need field data for iron loss calculation
    fea_config_dict['designer.OnlyTableResults'] = False 


if __name__ == '__main__':
    # add path to fea_config_dict
    where_am_i_v2(fea_config_dict)
