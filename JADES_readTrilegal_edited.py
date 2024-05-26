""" JADES_readTrilegal.py

    Reads in multiple trilegal files for specific filter

    3/10/24
"""
import pickle
import numpy as np
import sys

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


    max_length = max(len(data) for data in full_list)
    
    # Pad shorter lists with zeros
    padded_list = []
    for data in full_list:
        padding = max_length - len(data)
        padded_data = np.pad(data, (0, padding), mode='constant')
        padded_list.append(padded_data)

    trilegal_arr = np.array(padded_list)

    # Set printing options to display full array
    np.set_printoptions(threshold=sys.maxsize)   
    #trilegal_arr=np.array(full_list)
    
    with open(file_path_pic,'wb') as file:
        pickle.dump(trilegal_arr,file)

    with open(file_path_pic,'rb') as file:
        data=pickle.load(file)
    print(data)


file_path1='/Users/yuqifang/Desktop/JADES/trilegal/trilegal_output/02/0a2ae3a954d6dde3213d351f5897237f2d337077d37c2ef992f91e08bc0b46090.dat'
file_path2='/Users/yuqifang/Desktop/JADES/trilegal/trilegal_output/02/0a2ae3a954d6dde3213d351f5897237f2d337077d37c2ef992f91e08bc0b46091.dat'
file_path3='/Users/yuqifang/Desktop/JADES/trilegal/trilegal_output/02/0a2ae3a954d6dde3213d351f5897237f2d337077d37c2ef992f91e08bc0b46092.dat'
file_path4='/Users/yuqifang/Desktop/JADES/trilegal/trilegal_output/02/0a2ae3a954d6dde3213d351f5897237f2d337077d37c2ef992f91e08bc0b46093.dat'
file_path5='/Users/yuqifang/Desktop/JADES/trilegal/trilegal_output/02/0a2ae3a954d6dde3213d351f5897237f2d337077d37c2ef992f91e08bc0b46094.dat'
file_path6='/Users/yuqifang/Desktop/JADES/trilegal/trilegal_output/02/0a2ae3a954d6dde3213d351f5897237f2d337077d37c2ef992f91e08bc0b46095.dat'
file_path7='/Users/yuqifang/Desktop/JADES/trilegal/trilegal_output/02/0a2ae3a954d6dde3213d351f5897237f2d337077d37c2ef992f91e08bc0b46096.dat'
file_path8='/Users/yuqifang/Desktop/JADES/trilegal/trilegal_output/02/0a2ae3a954d6dde3213d351f5897237f2d337077d37c2ef992f91e08bc0b46097.dat'
file_path9='/Users/yuqifang/Desktop/JADES/trilegal/trilegal_output/02/0a2ae3a954d6dde3213d351f5897237f2d337077d37c2ef992f91e08bc0b46098.dat'
file_path10='/Users/yuqifang/Desktop/JADES/trilegal/trilegal_output/02/0a2ae3a954d6dde3213d351f5897237f2d337077d37c2ef992f91e08bc0b46099.dat'




file_path_pic='trilegal_test.pickle'
filt='F115W'
read_Trilegal(filt,[file_path1,file_path2,file_path3,file_path4,file_path5,file_path6,file_path7,file_path8,file_path9,file_path10],file_path_pic)












