""" JADES_readTrilegal.py

    reads in multiple trilegal files for filter list and dumps into pickle files for each filter

    6/8/24
"""
import pickle
import numpy as np
import os

def read_Trilegal(filters_list,directory_path,Ver):
    """read_Trilegal()
        #reads in multiple trilegal files for filter list and dumps into pickle files for each filter

    Args:
        filters_list ('list'): list of strings specifying filters to read
        directory_path ('string'): path to trilegal data
        Ver ('string'): data version to read
    """

    #create list of file paths
    dir_list = os.listdir(directory_path)
    trilegal_list = [f for f in dir_list if not f.startswith('.')]
    print("Files and directories in '", trilegal_list, "' :")
    # prints all files
    print(trilegal_list)

    #For all folders in trilegal.zip (all filter simulations)
    for i in range(len(trilegal_list)):
        full_list=[]

        #Find filter path of filter
        filter_path=str(directory_path)+'\\'+ str(trilegal_list[i])
        print('Reading ', trilegal_list[i])
        #Create list of all simulations
        sims_list=os.listdir(filter_path)
        
        #For every simulation
        for j in range(len(sims_list)):
            sim_path=str(filter_path)+'\\'+ str(sims_list[j])
            #print(sim_path)
            #open file
            file=open(sim_path)
            #print('Reading '+str(file_path))
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
                    file.close()
                    break
            
                #split each line to get values from inputted filter
                star=line.split()
            
                file_data.append(float(star[header.index(filters_list[i])]))
            #append list to full list
            full_list.append(file_data)

        max_length = max(len(data) for data in full_list)
    
        # Pad shorter lists with zeros
        padded_list = []
        for data in full_list:
            padding = max_length - len(data)
            padded_data = np.pad(data, (0, padding), mode='constant', constant_values=(np.nan))
            padded_list.append(padded_data)

        trilegal_arr=np.array(padded_list)
        #dump in pickle file
        with open(r'd:\JADES\ISL\Trilegal_pickles\trilegal_'+str(filters_list[i])+'_'+str(Ver)+'.pickle','wb') as file:
            pickle.dump(trilegal_arr,file)


#For v1-------------------
#ver='v1'
#dir_path= r'd:\JADES\ISL\trilegal_v1'
#filt_list_full=['F090W','F115W','F150W','F200W','F277W','F335M','F356W','F410M','F444W']
#---------------------------
#For v2------------------
#dir_path= r'd:\JADES\ISL\trilegal'
#filt_list_full=['F090W', 'F115W', 'F150W', 'F182M', 'F200W', 'F210M', 'F277W', 'F335M', 'F356W', 'F410M', 'F430M', 'F444W', 'F460M', 'F480M']
#-------------------------
#For v3------------------
ver='v3'
dir_path= r'd:\JADES\ISL\trilegal_v3'
filt_list_full=['F090W','F115W','F150W','F182M','F200W','F210M','F277W','F335M','F356W','F410M','F444W']
#-------------------------

read_Trilegal(filt_list_full,dir_path,ver)













