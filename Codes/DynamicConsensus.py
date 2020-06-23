# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 12:34:06 2020
@author: vxr131730 - Venkatraman Renganathan

This script simulates Consensus Algorithm
This script is tested in Python 3.7, Windows 10, 64-bit
(C) Venkatraman Renganathan, 2020.  Email: vrengana@utdallas.edu
    Karthik Ganapathy, 2020. Email: karthik.ganapathy@utdallas.edu

This program is a free software: you can redistribute it and/or modify it
under the terms of the GNU lesser General Public License, either version 
3.7, or any later version. This program is distributed in the hope that it 
will be useful, but WITHOUT ANY WARRANTY. 
"""

###############################################################################
###############################################################################

# Import all the required libraries
import numpy as np
import matplotlib.pyplot as plt

###############################################################################
###############################################################################

def SpoofResilientWMSR(N, T, F, x0):
    """
    Function spoofing_wmsr updates the information state of each
    vehicles after sorting & removing extreme values from its in-neighbors
    Input Parameters:
    N  : Number of agents
    T  : Time Span 
    F  : Number of malicious agents
    x0 : Initial value of agents
    """   
    # Define Data
    x = np.zeros((T+1,N))
    # Degree Matrix
    D = np.diag([2, 3, 4, 4, 4, 3, 2])
    # Adjacency Matrix
    A = np.array([[0, 1, 1, 0, 0, 0, 0],
                  [1, 0, 1, 1, 0, 0, 0],
                  [1, 1, 0, 1, 1, 0, 0],
                  [0, 1, 1, 0, 1, 1, 0],
                  [0, 0, 1, 1, 0, 1, 1],
                  [0, 0, 0, 1, 1, 0, 1],
                  [0, 0, 0, 0, 1, 1, 0]])
    # Laplacian matrix
    L = D - A
    
    # Set values of all vehicles at time = 0 to x_0
    x[0,:] = x0
    # Set malicious node value to random value
    x[:,3] = 50 + np.random.randn(T+1)
    
    for k in range(T):        
        for i in range(N):
            # Node #4 is malicious, so don't update that index
            if i == 3:
                continue
            
            # Extract the i th row
            L_i_row = L[:,i]            
            beforeSort = np.column_stack((x[k,:].T, L_i_row.T))
            
            # Extract only in-neighbors                            
            beforeSort = beforeSort[np.where(beforeSort[:,1]< 0),0][0,:]
                    
            # Removing larger values - sort descendingly
            ascendSort = np.sort(beforeSort)
            ascendSort = ascendSort[::-1]            
            indices    = np.nonzero(ascendSort > x[k,i])
            if not indices:
                if(len(indices) > F):
                    # if # of values larger than x(i) > F, delete F larger ones
                    for j in range(F):
                        np.delete(ascendSort, indices[j]) 
                        print("Hi mama")
                else:
                    # else delete all larger values
                    np.delete(ascendSort, indices)
                    print("Hi mama")
                    
            # Removing smaller values          
            ascendSort = np.sort(ascendSort)
            indices    = np.nonzero(ascendSort < x[k,i]) 
            if not indices:
                if(len(indices) > F):
                    for j in range(F):
                        np.delete(ascendSort, indices[j],axis=0)
                else:
                    # else delete all smaller values
                    np.delete(ascendSort, indices,axis=0)            
            
            weight = 1/(len(ascendSort)+1)
            
            # WMSR Update
            x[k+1,i] = np.sum(weight*ascendSort) + weight* x[k,i]         
    return x
                

###############################################################################
###############################################################################
###############################################################################

def main():    
    
    # Close any existing figure
    plt.close('all')   
    
    N  = 7  # Number of agents
    T  = 30 # Time Span 
    F  = 1  # Number of malicious agents    
    x0 = 50*np.random.rand(N)
    # Get the spoof resilient consensus updates for all time steps
    x  = SpoofResilientWMSR(N, T, F, x0)    
    # Plot the values
    timeVector = np.arange(T+1)
    plt.plot(timeVector, x)
    # Code to make the plot more readable
    plt.xlabel("Time Steps", fontsize=18)
    plt.ylabel("Agents States", fontsize=18)    
    plt.tick_params(axis='both', labelsize=18) 
    

    

###############################################################################

if __name__ == '__main__':
    main()
    
###############################################################################
###############################################################################
###################### END OF THE FILE ########################################
###############################################################################
###############################################################################