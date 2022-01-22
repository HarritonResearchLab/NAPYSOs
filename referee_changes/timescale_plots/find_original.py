import os
rootdir = r'C:\Users\Research\Documents\GitHub\NAPYSOs\recovered\2.0'

for subdir, dirs, files in os.walk(rootdir):
    for file in files:
        if file.split('.')[-1] == 'py': 
            joined_path = os.path.join(subdir, file)
            with open(joined_path, 'r') as f: 
                for line in f: 
                    if 'All Timescales' in line: 
                        print(joined_path)
                        break