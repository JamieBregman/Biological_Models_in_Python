#Assignment 4 Part 2

#Importing modules
import random
import string

#Target string
target="methinks it is like a weasel"

#Asking user for mutation rate per character per generation
#Works best with low mutation rate
mutation_rate=int(input("Enter mutation rate per character per generation: "))

#Generating a random sequence using lowercase letters & space
seq="".join(random.choice(string.ascii_lowercase + " ") for i in range (28))

#Creating a child sequence with user-input number of changes
def create_child(sequence, mutation):
    sequence=list(sequence)
    #Telling the function where in the sequence it should replace
    m = random.sample(range(0,28),mutation)
    #Telling the function what it can replace with
    for index in m:
        sequence[index] = random.choice(string.ascii_lowercase + " ")
    new_sequence="" .join(sequence)
    return str(new_sequence)

#Creating 50 children to increase efficiency of program
def get_children(sequence,mutation):
    seq_list=[]
    #Iterating 50 times to create 50 children
    for i in range(50):
        seq_list.append(create_child(sequence,mutation))
    return seq_list

#Function for measuring the differences between new sequence and target sequence
def score(sequence, target):
	score = 0
	#For every character in the sequence that matches the target, score increases by 1
	for i in range(len(target)):
		if sequence[i] == target[i]:
			score += 1
	return score

#Selecting the best of the 50 children relative to the target sequence
def best_child(sequence, target, mutation):
        #Bringing in the 50 chlidren
	seq_list = get_children(sequence, mutation)
	#Assigning best_seq
	best_seq = seq_list[0]
	#Assigning the value of best_score
	best_score = score(best_seq,target)
	#For every sequence in the list of 50 children...
	for seq in seq_list:
		score1 = score(seq, target)
		#Choosing the best score and assigning it to best_seq
		if score1 > best_score:
			best_score = score1
			best_seq = seq
	return best_seq

#Final loop function    
def weasel(seq,target,mutation_rate):
    generation=0
    s=score(seq,target)
    #Operates while the score of the sequence is less than the final score
    while s < len(target):
        s=score(seq,target)
        #Using the best of the 50 children
        new_best=best_child(seq,target,mutation_rate)
        #Assigning that best child a new score
        snew=score(new_best,target)
        if snew > s:
            seq = new_best
        if s >= snew:
            seq = seq
        generation += 1
        print("generation: ",generation)
        print(seq)
        
#Calling function
weasel(seq,target,mutation_rate)
