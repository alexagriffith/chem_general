
#!/usr/bin/python
import numpy as np
import os
import math
import csv
import system

coords = []
vals = []
energy =[]


atom1 = sys.argv[1]
atom2 = sys.argv[2]
atom3 = sys.argv[3]


def main():

    global coords
    global vals

    root = os.getcwd()
    subdirs = os.listdir(root)


    for subdir in subdirs:
       if os.path.isdir(subdir):
          file_dir = os.path.join(root,subdir)
        #  print(file_dir)
          os.chdir(file_dir)
          xyz_name = (str(subdir) + 'out.xyz')
          os.system( "t2x -c > {}".format(xyz_name))
       else:
          continue


       f=open(xyz_name,"r")
       f.readline()
       for line in f:
          line =line.split()
          coords.append(line[1:4])
          #print(coords)
       f.close()


       input_atoms = [atom1,atom2,atom3]
       input_one = input_atoms[0]
       input_two = input_atoms[1]
       input_three = input_atoms[2]

       coords1 = coords[input_one]
       coords2 = coords[input_two]
       coords3 = coords[input_three]

       coords1 = coords[input_one]
       coords2 = coords[input_two]
       coords3 = coords[input_three]

      # print(coords1)
       coords1 = np.array(coords1)
       coords2 = np.array(coords2)
       coords3 = np.array(coords3)

       #print(coords1)
       coords1 = coords1.astype(np.float)
       coords2 = coords2.astype(np.float)
       coords3 = coords3.astype(np.float)

      # print(coords1)
      # print(coords2)
      # print(coords3)
       a_vector= np.transpose(np.subtract(coords2,coords1))
       b_vector = np.subtract(coords3,coords2)

      # print(b_vector)
      # print(a_vector)
       ab_dot = np.float(np.sum(a_vector * b_vector))
       a_mag = np.float(np.linalg.norm(a_vector))
       b_mag = np.float(np.linalg.norm(b_vector))

       #print(a_mag)
       #print(b_mag)
       #print(ab_dot)
       theta = -( (ab_dot) /  (a_mag * b_mag) )
       #print(theta)
       theta = math.acos(theta)
       theta = (theta) * float(180 / math.pi)

       row_output = [int(subdir), str(theta)]
       vals.append(row_output)


       coords = []
       os.chdir(root)

    myfile = open('angles.csv', 'w')


    with myfile:
        writer = csv.writer(myfile)
        writer.writerows(sorted(vals, key= lambda x:int(x[0])))             #sorted by ints 



if __name__ == '__main__':
    main()
