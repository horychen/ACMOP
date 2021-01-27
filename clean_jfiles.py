
import os
import shutil
for el in os.listdir():
    # if '_pemd' in el:
    # if '_Q24p1y9_restart_from_optimal_and_reevaluate_wo_csv' in el:
    # if '_pemd2020_prolonged - Copy' in el:

    if os.path.isdir(el):
        for _el in os.listdir(el):
            if os.path.isdir(el+'/'+_el):
                # print( _el)

                for __el in os.listdir(el+'/'+_el):
                    if 'jfiles' in __el:
                        target = el+'/'+_el + '/' +__el
                        # if '_pemd2020/PMSM_Q24p1y9_A/proj99999.jfiles' in target:
                        # shutil.rmtree(target)
                        print('delete', target)

# use os.walk instead
def get_swarm_group(self, folder_of_collection):
    path = self.path2boptPython + folder_of_collection
    dict_path2SwarmDataOfTheSpecification = dict()
    dict_settingsOfTheSpecification = dict()
    list_specifications = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if 'swarm_data.txt' in file:
                # print(root)
                normpath = os.path.normpath(root)
                specification = normpath.split(os.sep)[-1].replace('_', ' ') # convert folder name to spec name
                list_specifications.append(specification)
                dict_path2SwarmDataOfTheSpecification[specification] = root + '/'
                # print(specification)
            if 'settings.txt' in file:
                with open(root+'/'+file, 'r') as f:
                    buf = f.read()
                    lst = buf.split('|')
                    specification = lst[0].strip()
                    select_fea_config_dict = lst[1].strip()
                    dict_settingsOfTheSpecification[specification] = select_fea_config_dict
    return list_specifications, dict_path2SwarmDataOfTheSpecification, dict_settingsOfTheSpecification
