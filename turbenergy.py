#!/bin/python
f=open("energy","r") #to create a new file, use "w" instead of "r"
data= f.read() #& use write
f.close() #close the energy file being read
data= data.split("SCFPOT\n")[1].split("$end")[0].split("\n")#split the data xyz using string "y" giving [x,z[ then keep [1] or z, then split [z1,z2] and keep [0]
#print data
for i,x in enumerate(data):
        data[i] = x[len("     1 "):]
data= "\n".join(data)
data= data.replace("      ",",").replace("     ",",")
data= "SCF,SCFKIN,SCFPOT\n"+data

column= int(input("column 0-2 give it to meh!"))
data2 = "\n".join([x.split(',')[column] for x in data.split("\n")[:-1:]]) #[#] for SCF column (0,1, or 2)
#to create a new file, use "w" instead of "r" ----> variable = open(filemane,"w")
# variable.write(data2)
# variable.close()
print data2
