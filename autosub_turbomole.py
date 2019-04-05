#!/usr/bin/env python
import os
import time

# run from sub directory

root = os.getcwd() #put script in angles directory / directory that encompasses files
sub_dirs = os.listdir(root)

for subdir in sub_dirs:             #like doing ls in the directory I'm in
        if os.path.isdir(subdir):
                run_dir = os.path.join(root, subdir) # adding path and specific run dir 
                os.chdir(run_dir)                    #change into the specific run directory
                os.system("cp ../*.pbs .")           #cp premade pbs file that will run all jobs 
                os.system("cp ../define.inp .")      # a define file containing all inputs used in turbomole define 
                os.system("x2t *.xyz > coord")       #taking xyz coordinate in file and converting it to coord 
        #run commands
                time.sleep(1)
                os.system("define < define.inp > define.out")           # creating a control file 
                time.sleep(1)
                os.system("qsub *.pbs")                                 #submitting job 
                time.sleep(1)

                os.chdir(root) #go back to run directory

        else:
                continue
