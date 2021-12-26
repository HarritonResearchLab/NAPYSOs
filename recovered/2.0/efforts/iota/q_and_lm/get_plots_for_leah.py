def get_plots(directory,new_dir,outfile):
    # Import(s)
    from os import listdir
    from os.path import join
    import math
    import shutil
    
    # Action
    onlyfiles = sorted([f for f in listdir(directory)])
    last_index = int(math.floor(len(onlyfiles)/2)-1)
    
    with open(outfile,'a') as f:
        for fn_index, file_name in enumerate(onlyfiles):
            if fn_index <= last_index: 
                object_id = file_name.split(':')[2].replace('.png','')
                f.write(object_id+'\n')
                shutil.copyfile(join(directory,file_name),join(new_dir,file_name))


old_dir = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/iota/q_and_lm/second/plots/labeled_by_q'
new_dir = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/iota/q_and_lm/second/plots/q_plots_for_leah'
out_file = '/home/thaddaeus/FMU/HRL/LAH2.0/efforts/iota/q_and_lm/second/plots/leah_class_file.csv'
get_plots(directory=old_dir,new_dir=new_dir,outfile=out_file)