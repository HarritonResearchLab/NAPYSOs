def new_key_for_lah(first_key,results_key,new_file):
    # Import(s)
    import pandas as pd

    # Action

    # That ~324 object result file
    results_df = pd.read_csv(results_key)
    results_ids = list(results_df['ID'])

    # That 696 key file
    old_df = pd.read_csv(first_key)
    old_ids = list(old_df['ID'])
    old_ras = list(old_df['RA'])
    old_decs = list(old_df['DEC'])

    # RAS and DECS for results objects
    new_ras = []
    new_decs = []
    for id in results_ids: 
        for id_index, old_id in enumerate(old_ids): 
            if id==old_id: 
                new_ras.append(float(old_ras[id_index]))
                new_decs.append(float(old_decs[id_index]))
                break

    results_df.insert(loc=1,column='RA',value=list(zip(new_ras)))
    results_df.insert(loc=2,column='DEC',value=list(zip(new_decs)))

    results_df.to_csv(new_file,index=False)

def remove_weird_stuff(file,n_file):
    # Import(s)

    # Action
    lines = []
    with open(file,'r') as f:
        for line in f: 
            line = line.replace('"(','')
            line = line.replace(')"','')
            line = line.replace(",,",',')
            lines.append(line.replace('\n',''))
    
    with open(n_file,'a') as f:
        for line in lines:
            f.write(line+'\n')

first_key = '/home/thaddaeus/FMU/HRL/LAH2.0/data/695/key.csv'
results_key = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/iota/q_and_lm/second/results.csv'
new_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/iota/q_and_lm/second/results_RA_AND_DEC.csv'
#new_key_for_lah(first_key=first_key,results_key=results_key,new_file=new_file)

remove_weird_stuff(file=new_file,n_file='/home/thaddaeus/FMU/HRL/LAH2.0/efforts/iota/q_and_lm/second/1results_RA_AND_DEC.csv')