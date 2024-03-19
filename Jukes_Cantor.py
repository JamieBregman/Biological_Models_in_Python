#Assignment 5 Part 2

#Importing modules
import math

#Function that reads the file and loads the species names into a list of strings and loads the sequences into a list of strings
def getlistofsequences(filename):
    f=open(filename,"r")
    f_read=f.read()
    f.close()
    f_split=f_read.split()
    species_list=f_split[::2]
    seq_list=f_split[1::2]
    return seq_list, species_list

#Using the above function and separating data into species names and sequences
get_lists=getlistofsequences("Myotis_aligned.fa")
seq_list=get_lists[0]
species_list=get_lists[1]

#Abbreviating each species name to its first 9 letters so the results file is less clustered
species_list[:] = (elem[:9] for elem in species_list)

#Function that calculates the proportion of bases that are different and returns a value
def calculatep(seq1,seq2):
    S=0
    D=0
    #Making sure not to count "-" as matches or mismatches
    for i in range(len(seq1)):
        if seq1[i] == seq2[i] and seq1[i] != "-" and seq2[i] != "-":
            S +=1
        elif seq1[i] != seq2[i] and seq1[i] != "-" and seq2[i] != "-":
            D +=1
    p=D/(S+D)
    return p

#Function that calculates K using the Jukes-Cantor formula and a given value of p
def calculateKfromp(p):
    K=(-3/4)*(math.log(1-((4*p)/3)))
    return K

#Calculating p values for each sequence against one another
p_values=[]
for i in range(len(seq_list)):
    for j in range(len(seq_list)):
        p_values.append(calculatep(seq_list[i],seq_list[j]))

#Calculating K for each sequence against each sequence and creating a list
k_values=list(map(calculateKfromp,p_values))

#Separating the k values into a list of lists where each "list" represents a "row" in the final matrix
k_values_split = [k_values[i:i+75] for i in range(0, len(k_values), 75)]

#Creating a matrix with the species names as the column/row names and filling it with the k_values
#Fills each row of the array with each "list" from the list of lists in k_values_split
array = [k_values_split[i] for i in range(len(species_list))]
#Adding the species name to each row
array = [species_list] + array
#Adding the column names to the matrix
column_names = ["JC Distance Matrix"] + species_list
#Tying it all together by combining the column names with array containing row names
#Adding a length of "+1" to make room for the empty space of the title 
new_array = [[column_names[i]]+array[i] for i in range(len(species_list)+1)]
matrix=[]
#Final matrix that contains the sequences with the correct column/row names
for i in new_array:
       matrix.append(i)

#Writing the matrix to a results file
with open("matrix_results.txt","w") as testfile:
    for row in matrix:
          testfile.write(" ".join([str(i) for i in row]) + "\n")
