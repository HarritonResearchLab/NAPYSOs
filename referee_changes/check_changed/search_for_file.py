import os

path_list = [os.path.join(dirpath,filename) for dirpath, _, filenames in os.walk('.') for filename in filenames if filename.endswith('.py')]


# the file is recovered\2.0\efforts\august21\aug23-29\froebrich\froebrich_periods.py

search_string = '2MASS_J20582381+4353114'

for file_path in path_list: 
    with open(file_path, 'r') as f: 
        for line in f: 
            if search_string in line: 
                print(file_path)
                break 