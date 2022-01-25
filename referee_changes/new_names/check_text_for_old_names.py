def check_test(ids_to_check): 
    import pandas as pd
    import numpy as np

    df = pd.read_csv('./referee_changes/new_names/renamed.csv')
    old_names = np.array(df['ID'])
    new_names = np.array(df['preferred_name'])

    for test_id in ids_to_check: 

        for old_id, new_id in zip(old_names, new_names): 
            if test_id==old_id: 
                if test_id != new_id: 
                    print(test_id, '###', new_id)
                    break

        if test_id not in old_names: 
            print("ERROR: "+test_id)

def edit_whole_file(): 
    import pandas as pd
    import numpy as np

    df = pd.read_csv('./referee_changes/new_names/renamed.csv')
    old_names = np.array(df['ID'])
    new_names = np.array(df['preferred_name'])

    changes = 0

    with open('./referee_changes/check_changed/raw_tex.txt', 'r') as f:
        new_list_of_lines = []

        for line in f:
            for old_id, new_id in zip(old_names, new_names): 
                old_id_space = old_id.replace('_',' ')
                new_id_space = new_id.replace('_', ' ')
                if old_id in line: 
                    if old_id != new_id: 
                        print(old_id,'####', new_id, '@@@@', line)
                        #line = line.replace(old_id, new_id)
                        changes += 1
                if old_id_space in line: 
                    if old_id_space != new_id_space: 
                        #line = line.replace(old_id_space, new_id_space)
                        print(old_id_space,'####', new_id_space, '@@@@', line)
                        #changes += 1

            new_list_of_lines.append(line)
    '''
    with open('./referee_changes/check_changed/raw_tex.txt', 'w') as f:
        for new_line in new_list_of_lines: 
            f.write(new_line)
    '''
    print('number of changes: ', changes)

#edit_whole_file()

ids_to_check =['GDR1_2162872928138757248', '2MASS_J20523394+4429168', 
               'FHK_577', 'GDR1_2162251394832460672', 'GDR1_2163140315626343168',
               'FHK_348'] # issues noted


fig_14_ids = ['FHK_176', 'FHK_286'] # good

fig_8_ids = ['FHK_32', '2MASS_J20512059+4420322'] # issues noted

fig_7_ids = ['FHK_163', '2MASS_J20582381+4353114'] # issues noted

check_test(fig_7_ids)