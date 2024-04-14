""" JADES_readTrilegal.py

    Reads in multiple trilegal files for specific filter

    3/10/24
"""
import pickle
import numpy as np

def read_Trilegal(filt,file_list, file_path_pic):
    #full list of several trilegal files data for specific filter
    full_list=[]
    #for every file
    for file_path in file_list:
        #open file
        file=open(file_path)
        print('Reading '+str(file_path))
        #list to append looped file data
        file_data=[]

        #while there are lines to read
        while (True):
            #read each line
            line=file.readline()
            
            #if first line in file
            if line.startswith('#Gc'):
                #split header
                header=line.split()
                
                #read next line
                line=file.readline()

            #break loop at end
            if line.startswith('#TRILEGAL'):
                print('closed')
                file.close()
                break
            
            #split each line to get values from inputted filter
            star=line.split()
            
            file_data.append(float(star[header.index(filt)]))

        #append list to full list
        full_list.append(file_data)
    
    trilegal_arr=np.array(full_list)
    
    with open(file_path_pic,'wb') as file:
        pickle.dump(trilegal_arr,file)

    with open(file_path_pic,'rb') as file:
        data=pickle.load(file)
    print(data)

file_path1=r'c:\Users\kasia\Downloads\smaller_trilegal_trial.dat'
file_path2=r'c:\Users\kasia\Downloads\smaller_trilegal_trial.dat'
file_path_pic='trilegal_test.pickle'
filt='F090W'
read_Trilegal(filt,[file_path1,file_path2],file_path_pic)