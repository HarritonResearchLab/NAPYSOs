def display_and_record(key):
    # Import(s)
    import matplotlib.pyplot as plt
    import pandas as pd
    import matplotlib.image as mpimg
    import time

    # Action
    #key_df = pd.read_csv(key)
    #ids = list(key_df['ID'])

    ids = ['FHK_445','FHK_449']
    classifications = []
    for id in ids: 
        img_path = '+++.png'.replace('+++',id)
        img = mpimg.imread(img_path)
        imgplot = plt.imshow(img)
        plt.show()
        plt.clf()
        img_class = input("Image class: ")
        classifications.append(img_class)
        
    for id, classif in zip(ids,classifications):
        print(id+': '+classif)

    for i in range(12):
        bap = input("ID ****** class: ")

display_and_record(key='')
