#Assignment 3 Part 2
#Program to implement the Island Model of Migration
#Gathering information from the user
P=float(input("Enter P frequency on mainland: "))
p=float(input("Enter starting p frequency on the island: "))
m=float(input("Enter fraction of all genes on island that come from mainland each generation (must be less than or equal to 0.01): "))
s=int(input("Enter random number seed (must be a positive integer): "))

#Importing random module
import random

#Function for calculating allele frequency
def p_subt_calc(a,b,c,d):
    #Making sure that the parameters for m are met
    if c > 0.01 or c < 0:
        return print("Error: m value must be greater than zero and less than or equal to 0.01")
    #Making sure that the parameters for s are met
    if d < 0:
        return print("Error: s value must be a positive integer")
    #Starting the generation count from 0
    generation=0
    #Using the user data for the seed value
    random.seed(d)
    #While the value for allele frequency (b) is between 0 and 1 the function calculates frequency using the standard migration model.
    #It also uses random.unfirom() to add a new random value based on the seed for each generation
    while b<1 or b>0:
        b = ((c*(a - b))+b) + random.uniform(-0.05,0.05)
        #Ending the function once the allele is fixed and returning the corresponding statement
        if b >=1:
            return print("The allele frequency is fixed with value " + str(b) + " after " + str(generation) + " generations")
        #Ending the function once the allele is lost and returning the corresponding statement
        if b <=0:
            return print("The allele frequency is lost with value " + str(b) + " after " + str(generation) + " generations")
        #Ending the function if the process is taking longer than 10000 generations and returning the corresponding statement
        if generation == 10000:
            return print("Process is taking too long. Current allele frequency is " + str(b) + ".")
        #Instructing the function to increase by 1 "generation" for each calculation
        generation += 1

#Running the function
p_subt_calc(P,p,m,s)
