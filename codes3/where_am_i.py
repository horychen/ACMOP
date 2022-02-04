import os

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

    # dir_interpreter = os.path.abspath('') + '/'
    # if 'DrH' not in dir_interpreter:
    #     # print(dir_interpreter)
    #     raise Exception("Revert this back before you run on server.")
    # 先用这个救救急，到服务器的时候这里要报错
    # dir_interpreter = r'D:\DrH\bopt-python\codes3/'
    # dir_parent = dir_interpreter #+ '../' # '../' <- relative path is not always a good option
    # if not bool_post_processing:
    #     if 'codes3' not in dir_interpreter:
    #         raise Exception('Do not run the script from your current directory: %s.\nPlease cd to ./codes3 before running python scripts.'%(dir_interpreter))

    # A = os.path.join(os.path.dirname(__file__), '..')
    # # A is the parent directory of the directory where program resides.
    # B = os.path.dirname(os.path.realpath(__file__))
    # # B is the canonicalised (?) directory where the program resides.
    # C = os.path.abspath(os.path.dirname(__file__))
    # # C is the absolute path of the directory where the program resides.
    # print(__file__)
    # print(os.path.join(os.path.dirname(__file__), '..'))
    # print(os.path.dirname(os.path.realpath(__file__)))
    # print(os.path.abspath(os.path.dirname(__file__)))

    dir_parent = os.path.join(os.path.dirname(__file__), '..') + '/'
    dir_codes  = os.path.abspath(os.path.dirname(__file__))
    # dir_femm_files = dir_parent + 'femm_files/'
    pc_name = get_pc_name()
    os.chdir(dir_codes)
    print('[where_am_i.py] CD to:', dir_codes)

    # print(dir_parent, dir_codes, pc_name, sep='\n'); quit()
    # fea_config_dict['dir.interpreter']       = dir_interpreter 
    fea_config_dict['dir.parent']            = dir_parent
    # fea_config_dict['dir.codes']             = dir_codes
    # fea_config_dict['dir.femm_files']        = dir_parent
    fea_config_dict['pc_name']               = pc_name

    fea_config_dict['delete_results_after_calculation'] = False # True for saving disk space (but you lose voltage profile and element data)
    if 'designer.OnlyTableResults' in fea_config_dict.keys():
        # if fea_config_dict['designer.Restart'] == False:
        #     fea_config_dict['designer.OnlyTableResults'] = True  # save disk space for my PC
        # However, we need field data for iron loss calculation
        fea_config_dict['designer.OnlyTableResults'] = False 


if __name__ == '__main__':
    # add path to fea_config_dict
    where_am_i_v2(fea_config_dict)
