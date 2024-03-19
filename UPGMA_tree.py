#Assignment 6 Part 1

#First section of code is from assignment 5 (constructing a distance matrix)

#Importing modules
import math

#Function that reads the file and loads the species & the sequences into separate lists of strings
def getlistofsequences(filename):
    f=open(filename,"r")
    f_read=f.read()
    f.close()
    #splitting the file text
    f_split=f_read.split()
    #splitting the list into species & sequences
    species_list=f_split[::2]
    seq_list=f_split[1::2]
    return seq_list, species_list

#Using the above function and separating data into species names and sequences
get_lists=getlistofsequences("Myotis_aligned.fa")
seq_list=get_lists[0]
species_list=get_lists[1]

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

#Calculating p for every cell in the "matrix"
p_values=[]
for i in range(len(seq_list)):
    for j in range(len(seq_list)):
        p_values.append(calculatep(seq_list[i],seq_list[j]))

#Calculating K for every cell in the "matrix" using the p values
k_values=list(map(calculateKfromp,p_values))

#Separating the k values and creating a completed matrix where each "list" represents a "row"
k_values_split = [k_values[i:i+75] for i in range(0, len(k_values), 75)]

#Creating lower triangular matrix from the full matrix generated above
mylist=[]
for i, j in enumerate(k_values_split[1:]):
    mylist.append(j[:i+1])

table=[[]]+mylist

#Second section of code begins here; building a UPGMA tree w/ Newick format

#Function that locates the smallest cell in the table
#Function that locates the lowest value in the table
def minimum(table):
    # setting the initial minimum value to 1 
    minimum = 1
    x=0
    y=0

    #Iterates through every value until it finds the minimum
    for i in range(len(table)):
        for j in range(len(table[i])):
            #if value is lower than the minimum, set that value to be the new minimum
            if table[i][j] < minimum:
                minimum = table[i][j]
                #creating "coordinates" for the lowest value in the table
                x, y = i, j
    return x, y

#Function that creates a node using the species names (rows)
def create_node(species_list, x, y):
    if y < x:
        x,y = y,x
    #Formatting the node in newick format
    species_list[x] = "(" + species_list[x] + "," + species_list[y] + ")"
    # Removing the index per UPGMA algorithm (redundant)
    del species_list[y]

#Function that modifies the other values in the table by finding their average relative to the minimum
def modify_table(table, x, y):
    #Making sure the indexes are in the appropriate order
    if y < x:
        x,y = y,x
    #Iterate over the values in the row x where i < x
    row = []
    for i in range(0, x):
        #modifying the row so that the new values are updated per UPGMA algorithm
        row.append((table[x][i] + table[y][i])/2)
    #Replacing the row x with a new row of updates values
    table[x] = row
    
    #Iterating over the values in the row x where i > x
    for i in range(x+1, y):
        #modifying the row so that the new values are updated per UPGMA algorithm
        table[i][x] = (table[i][x]+table[y][i])/2
        
    #Iterating over values in row i
    for i in range(y+1, len(table)):
        #modifying the row so that the new values are updated per UPGMA algorithm
        table[i][x] = (table[i][x]+table[i][y])/2
        #Removing the index per UPGMA algorithm (now that new values are established)
        del table[i][y]

    # Removing the index row per UPGMA algorithm (redundant)
    del table[y]

#Function that combines the UPGMA steps and creates a node to output a labelled table in Newick format
def UPGMA(table, species_list):
    """Loop finds the minimum value, modifies the rest of the table in reference to that value and
        according to UPGMA process, then creates nodes in Newick format from that data"""
    #Loop runs until every species in species_list has been joined
    while len(species_list) > 1:
        x, y = minimum(table)
        modify_table(table, x, y)
        create_node(species_list, x, y)

    #Printing the final tree
    return species_list[0] + ";"

#Running the function
UPGMA_output = UPGMA(table,species_list)

#Writing the Newick string to a file
file=open("Assignment6_File.txt","w")
file.write(str(UPGMA_output))
file.close()
