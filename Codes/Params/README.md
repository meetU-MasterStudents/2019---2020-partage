Parameter files needed for project
==================================

Statistical potential file `dope.par` 
-------------------------------------
Statistical potential values. Contains energy value for a pair of specific atoms at a given bin distance of 0.5 Angstroms width e.g : first value is for a pair of atoms at a distance between 0.25 and 0.75 A, second value is for a pair of atoms at a distance between 0.75 and 1.25 A and so on until reach 15.0 A.  
  
Format:  
 
     +----------------------- First  Residu pair   
     |  +-------------------- First  atom pair   
     |  |   +---------------- Second Residu pair   
     |  |   |  +------------- Second atom pair   
     |  |   |  |   +--------- Energy value for the first bin (0.25 to 0.75 A)    
     |  |   |  |   |     +--- Energy value for the second bin (0.75 to 1.25 A)    
     V  V   V  V   V     V    
    ALA CA ALA C 10.00 10.00 10.00 10.00 10.00 10.00 10.00 3.32 -1.11 -0.98 -0.52 -0.54 0.26 0.16 -0.36 0.16 0.20 -0.07 -0.16 -0.20 -0.14 -0.14 -0.10 -0.06 -0.03 -0.04 -0.07 -0.02 -0.02 -0.02  


For the study use only CA - CA pair energies  
