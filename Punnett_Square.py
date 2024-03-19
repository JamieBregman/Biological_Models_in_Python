#Biological Models in Python Final Project

#Importing modules
import itertools
import pandas as pd

#FUNCTIONS OUTSIDE THE CLASS
#Function that finds unique parental allele combinations from a given genotype
def get_parent_pairs(parent_split):

        #Empty list
        pairs_list = []

        #Monohybrid cross
        if len(parent_split[0]) == 1:
           pairs_list = parent_split
           return pairs_list

        #Dihybrid cross
        elif len(parent_split[0]) == 2 and len(parent_split) == 2:
            for i in parent_split:
                if i == parent_split[0]:
                    for j in i:
                        #Pairing each allele in the first part of the parental genotype with the second part
                        pairs_list.append(j+parent_split[1][0])
                        pairs_list.append(j+parent_split[1][1])
            return pairs_list

        #Trihybrid cross
        elif len(parent_split[0]) == 2 and len(parent_split) == 3:
            for i in parent_split:
                if i == parent_split[0]:
                    for j in i:
                        #Pairing each allele in the first part of the parental genotype with the second part
                        pairs_list.append(j+parent_split[1][0]+parent_split[2][0])
                        pairs_list.append(j+parent_split[1][0]+parent_split[2][1])
                        pairs_list.append(j+parent_split[1][1]+parent_split[2][0])
                        pairs_list.append(j+parent_split[1][1]+parent_split[2][1])
            return pairs_list


#CALCULATING GENOTYPIC FREQUENCY
#Function that gets the frequency of each genotype in the offspring (but does not account for duplicates)
def get_frequency(genotype_list, genotype):
    counter_list = []
    counter = 0
    
    #Sorted() is case-sensitive so it can distinguish between dominant & recessive alleles
    for i in genotype_list:
            if sorted(i) == sorted(genotype):
                counter += 1
                count_value = [genotype,counter]
    counter_list.append(count_value)
    return counter_list

#Function that sorts a list of offspring, gets the counts of each genotype, and then removes duplicates to obtain a final %
def genotypic_frequencies(offspring_list):

    #Sorting the offspring using sorted() so that they can be compared for duplicates
    offspring_list_final = ["".join(i) for i in [sorted(i) for i in offspring_list]]

    #Creating a list that has each genotype's frequency in the square, but does not account for duplicates
    frequency_list = [item for sublist in [get_frequency(["".join(i) for i in [sorted(i) for i in offspring_list]],i) for i in ["".join(i) for i in [sorted(i) for i in offspring_list]]] for item in sublist]

    #Creating a dictionary with each genotype's frequency so that the appropriate frequencies can be referenced later after duplicates are removed 
    dictionary = {t[0]:t[1] for t in frequency_list}

    #Using itertools.groupby to remove duplicates
    #First sort() the list of offspring
    offspring_list_final.sort()
    #Remove duplicates
    non_duplicates = list(offspring_list_final for offspring_list_final,_ in itertools.groupby(offspring_list_final))

    #Accessing the counts of the non-duplicate genotypes after removing them in the step above using the dictionary created
    counts = [dictionary[key] for key in non_duplicates]

    
    #Creating list with non-duplicate genotypes and their corresponding counts in the punnett square
    frequency_list_final = [list(i) for i in list(zip(non_duplicates,counts))]

    #Replacing the "counts" with actual %'s
    for i in frequency_list_final:
        i[1] = str(((i[1]/len(offspring_list))*100)) + "%"

    #Printing the final list
    print("\n" + "Genotypic frequencies of your cross:" + "\n")
    for i in frequency_list_final:
        print(*i)


#CALCULATING PHENOTYPIC FREQUENCY
#Function to count the amount of uppercase & lowercase letters in an offspring
def count_case(genotype):
    lowercase = 0
    uppercase = 0
    for i in genotype:
        if i.islower():
            lowercase += 1
        else:
            uppercase += 1
    #Returning the counts in this specific format
    return [lowercase,uppercase]

#Function that calculates the percentage of dominant & recessive genotypes in the offspring  
def mono_phenotypic_frequencies(offspring_list):
    dominant_or_recessive_list = []
    dominant = 0
    recessive = 0

    #Using the count case function for every child in the offspring list
    for i in offspring_list:
        dominant_or_recessive_list.append(count_case(i))

    #i[0] == 2 indicates recessive because of the format in which we output the lowercase/uppercase counts in the count_case() function above
    for i in dominant_or_recessive_list:
        if i[0] == 2:
            recessive += 1
        else:
            dominant += 1

    #Printing the phenotypic frequencies as %'s
    print("\n" + "Phenotypic frequencies of your cross:" + "\n")
    print("Recessive: " + str(((recessive/4)*100)) + "%" + "\n" +
          "Dominant: " + str(((dominant/4)*100)) + "%")

def di_phenotypic_frequencies(offspring_list):
    #Creating counting variables for each possible phenotype
    dominant_dominant = 0
    dominant_recessive = 0
    recessive_dominant = 0
    recessive_recessive = 0

    #Creating a final allele list containing every offspring allele
    alleles = list(set(offspring_list[0].lower()))
    alleles_uppercase = []
    for i in alleles:
        alleles_uppercase.append(i.upper())
    alleles_final = sorted(alleles + alleles_uppercase)

    #Collection of if statements running in a for loop that count the occurrence of each phenotype
    for i in offspring_list:
        if i.count(alleles_final[0]) >= 1:
            if i.count(alleles_final[1]) >= 1:
                dominant_dominant += 1
            else:
                dominant_recessive += 1
        else:
            if i.count(alleles_final[1]) >= 1:
                recessive_dominant += 1
            else:
                recessive_recessive += 1

    #Printing the phenotypic frequencies as %'s
    print("\n" + "Phenotypic frequencies of your cross:" + "\n")
    print("Dominant/Dominant: " + str(((dominant_dominant/16)*100)) + "%" + "\n" +
          "Dominant/Recessive: " + str(((dominant_recessive/16)*100)) + "%" + "\n" +
          "Recessive/Dominant: " + str(((recessive_dominant/16)*100)) + "%" + "\n" +
          "Recessive/Recessive: " + str(((recessive_recessive/16)*100)) + "%")

def tri_phenotypic_frequencies(offspring_list):
    #Creating counting variables for each possible phenotype
    dom_dom_dom = 0
    dom_dom_rec = 0
    dom_rec_dom = 0
    dom_rec_rec = 0
    rec_dom_dom = 0
    rec_dom_rec = 0
    rec_rec_dom = 0
    rec_rec_rec = 0

    #Creating a final allele list containing every offspring allele
    alleles = list(set(offspring_list[0].lower()))
    alleles_uppercase = []
    for i in alleles:
        alleles_uppercase.append(i.upper())
    alleles_final = sorted(alleles + alleles_uppercase)

    #Collection of if statements running in a for loop that count the occurrence of each phenotype
    for i in offspring_list:
        if i.count(alleles_final[0]) >= 1:
            if i.count(alleles_final[1]) >= 1:
                if i.count(alleles_final[2]) >= 1:
                    dom_dom_dom += 1
                else:
                    dom_dom_rec += 1
            else:
                if i.count(alleles_final[2]) >= 1:
                    dom_rec_dom += 1
                else:
                    dom_rec_rec += 1
        else:
            if i.count(alleles_final[1]) >= 1:
                if i.count(alleles_final[2]) >= 1:
                    rec_dom_dom += 1
                else:
                    rec_dom_rec += 1
            else:
                if i.count(alleles_final[2]) >= 1:
                    rec_rec_dom += 1
                else:
                    rec_rec_rec += 1

    #Printing the phenotypic frequencies as %'s
    print("\n" + "Phenotypic frequencies of your cross:" + "\n")
    print("Dominant/Dominant/Dominant: " + str(((dom_dom_dom/64)*100)) + "%" + "\n" +
          "Dominant/Dominant/Recessive: " + str(((dom_dom_rec/64)*100)) + "%" + "\n" +
          "Dominant/Recessive/Dominant: " + str(((dom_rec_dom/64)*100)) + "%" + "\n" +
          "Dominant/Recessive/Recessive: " + str(((dom_rec_rec/64)*100)) + "%" + "\n" +
          "Recessive/Dominant/Dominant: " + str(((rec_dom_dom/64)*100)) + "%" + "\n" +
          "Recessive/Dominant/Recessive: " + str(((rec_dom_rec/64)*100)) + "%" + "\n" +
          "Recessive/Recessive/Dominant: " + str(((rec_rec_dom/64)*100)) + "%" + "\n" +
          "Recessive/Recessive/Recessive: " + str(((rec_rec_rec/64)*100)) + "%")


#Producing punnett square from parental genotypes
class PredictingOffspring:
    #Initializing the class with 2 variables (parent1 and parent 2)
    def __init__(self, parent1, parent2):
        self.parent1 = parent1
        self.parent2 = parent2

    #Calling a punnett square function using the two class variables
    def punnett(self):
        #Monohybrid cross
        if len(self.parent1) and len(self.parent2) == 2:

            #Creating lists with all possible PARENTAL allele combinations
            p1_gene_list = [self.parent1[i:i+1] for i in range(0, len(self.parent1), 1)]
            p2_gene_list = [self.parent2[i:i+1] for i in range(0, len(self.parent2), 1)]

            p1_final = column_list = p1_gene_list
            p2_final = row_list = p2_gene_list

            punnett_values = []

            #Using itertools to find all possible OFFSPRING allele combinations
            for pair in itertools.product(p1_final,p2_final):
                punnett_values.append("".join(pair))

            #Formatting offspring list/values so they can be read in a dataframe
            punnett_split = [punnett_values[i:i+2] for i in range(0, len(punnett_values), 2)]
            punnett_dict = {k:v for k, v in zip(column_list, punnett_split)}
            df = pd.DataFrame(punnett_dict, index = row_list)
            print("\n")
            print(df)

            #Printing genotypic frequencies
            genotypic_frequencies(punnett_values)

            #Printing phenotypic frequencies
            mono_phenotypic_frequencies(punnett_values)

        #Dihybrid or trihybrid cross
        elif len(self.parent1) and len(self.parent2) == 4 or 6:

            #Creating lists with parental alleles separated
            p1_gene_list = [self.parent1[i:i+2] for i in range(0, len(self.parent1), 2)]
            p2_gene_list = [self.parent2[i:i+2] for i in range(0, len(self.parent2), 2)]

            #Creating lists with all possible PARENTAL allele combinations
            p1_final= column_list = get_parent_pairs(p1_gene_list)
            p2_final = row_list = get_parent_pairs(p2_gene_list)

            punnett_values = []

            #Using itertools to find all possible OFFSPRING allele combinations
            for pair in itertools.product(p1_final,p2_final):
                punnett_values.append("".join(pair))

            #Dihybrid cross
            if len(self.parent1) and len(self.parent2) == 4:

                #Formatting offspring list/values so they can be read in a dataframe
                punnett_split = [punnett_values[i:i+4] for i in range(0, len(punnett_values), 4)]
                punnett_dict = {k:v for k, v in zip(column_list, punnett_split)}
                df = pd.DataFrame(punnett_dict, columns = column_list, index = row_list)
                print("\n")
                print(df)

                #Printing genotypic frequencies
                genotypic_frequencies(punnett_values)

                #Printing phenotypic frequencies
                di_phenotypic_frequencies(punnett_values)

            #Trihybrid cross
            elif len(self.parent1) and len(self.parent2) == 6:

                #Formatting offspring list/values so they can be read in a dataframe
                punnett_split = [punnett_values[i:i+8] for i in range(0, len(punnett_values), 8)]
                punnett_dict = {k:v for k, v in zip(column_list, punnett_split)}
                df = pd.DataFrame(punnett_dict, columns = column_list, index = row_list)
                print("\n")
                print(df)

                #Printing genotypic frequencies
                genotypic_frequencies(punnett_values)

                #Printing phenotypic frequencies
                tri_phenotypic_frequencies(punnett_values)

#FUNCTIONS FOR PREDICTING PARENTS
mono_parent_list = ['AA', 'Aa', 'aA', 'aa']
di_parent_list = ['ABAB', 'ABAb', 'ABaB', 'ABab', 'AbAB', 'AbAb', 'AbaB', 'Abab',
                  'aBAB', 'aBAb', 'aBaB','aBab', 'abAB', 'abAb', 'abaB', 'abab']
tri_parent_list = ['ABCABC', 'ABCABc', 'ABCAbC', 'ABCAbc', 'ABCaBC', 'ABCaBc', 'ABCabC', 'ABCabc',
                   'ABcABC', 'ABcABc', 'ABcAbC', 'ABcAbc', 'ABcaBC', 'ABcaBc', 'ABcabC', 'ABcabc',
                   'AbCABC', 'AbCABc', 'AbCAbC', 'AbCAbc', 'AbCaBC', 'AbCaBc', 'AbCabC', 'AbCabc',
                   'AbcABC', 'AbcABc', 'AbcAbC', 'AbcAbc', 'AbcaBC', 'AbcaBc', 'AbcabC', 'Abcabc',
                   'aBCABC', 'aBCABc', 'aBCAbC', 'aBCAbc', 'aBCaBC', 'aBCaBc', 'aBCabC', 'aBCabc',
                   'aBcABC', 'aBcABc', 'aBcAbC', 'aBcAbc', 'aBcaBC', 'aBcaBc', 'aBcabC', 'aBcabc',
                   'abCABC', 'abCABc', 'abCAbC', 'abCAbc', 'abCaBC', 'abCaBc', 'abCabC', 'abCabc',
                   'abcABC', 'abcABc', 'abcAbC', 'abcAbc', 'abcaBC', 'abcaBc', 'abcabC', 'abcabc'] 

def calculate_possible_genotypes(offspring_phenotype_list):
    if len(offspring_phenotype_list[0]) == 1:
        for j,i in enumerate(offspring_phenotype_list):
            if i == "d":
                offspring_phenotype_list[j] = ["AA", "Aa"]
            elif i == "r":
                offspring_phenotype_list[j] = ["aa"]

    if len(offspring_phenotype_list[0]) == 2:
        for j,i in enumerate(offspring_phenotype_list):
            if i == "dd":
                offspring_phenotype_list[j] = ["ABAB", "ABAb", "ABaB", "ABab", "AbAB", "AbaB", "abAB"]
            elif i == "dr":
                offspring_phenotype_list[j] = ["AbAb", "Abab", "abAb"]
            elif i == "rd":
                offspring_phenotype_list[j] = ["aBaB", "aBab", "abaB"]
            elif i == "rr":
                offspring_phenotype_list[j] = ["abab"]

    if len(offspring_phenotype_list[0]) == 3:
        for j,i in enumerate(offspring_phenotype_list):
            if i == "ddd":
                offspring_phenotype_list[j] = ['ABCABC', 'ABCABc', 'ABCAbC', 'ABCAbc', 'ABCaBC', 'ABCaBc',
                                               'ABCabC', 'ABCabc', 'ABcABC', 'AbCABC', 'AbcABC', 'aBCABC',
                                               'aBcABC', 'abCABC', 'abcABC', 'AbCABc', 'aBCABc', 'abCABc',
                                               'ABcAbC', 'aBCAbC', 'aBcAbC', 'aBCAbc', 'ABcaBC', 'AbCaBC',
                                               'AbcaBC', 'AbCaBc', 'ABcabC']
            if i == "ddr":
                offspring_phenotype_list[j] = ['ABcABc', 'AbcABc', 'aBcABc', 'abcABc', 'ABcAbc', 'aBcAbc',
                                               'ABcabc', 'ABcaBc', 'AbcaBc']
            if i == "drd":
                offspring_phenotype_list[j] = ['AbCAbC', 'AbCAbc', 'AbcAbC', 'abCAbC', 'abcAbC', 'abCAbC',
                                               'AbCabc', 'AbCabC', 'AbCabC']
            if i == "drr":
                offspring_phenotype_list[j] = ['AbcAbc', 'abcAbc', 'Abcabc']
            if i == "rdd":
                offspring_phenotype_list[j] = ['aBCaBC', 'aBCaBc', 'aBCabC', 'aBCabc', 'aBcaBC', 'abCaBC',
                                               'abcaBC', 'abCaBc', 'aBcabC']
            if i == "rdr":
                offspring_phenotype_list[j] = ['aBcaBc', 'abcaBc', 'aBcabc']
            if i == "rrd":
                offspring_phenotype_list[j] = ['abCabC', 'abcabC', 'abCabc']
            if i == "rrr":
                offspring_phenotype_list[j] = ['abcabc']

    return offspring_phenotype_list
        

def find_parents(offspring_list_of_lists):
    offspring_counter = 0
    
    for a in offspring_list_of_lists:
        offspring_counter += 1

        for i in a:
            
            print("offspring",offspring_counter,"(potential genotype: ",i + ")")
                  
            first_half = i[0:len(i)//2]
            second_half = i[len(i)//2:]

            poss_parent1 = []
            poss_parent2 = []

            if len(i) == 2:
                for j in mono_parent_list:
                    if first_half in j:
                        poss_parent1.append(j)
                for j in mono_parent_list:
                    if second_half in j:
                        poss_parent2.append(j)
                        
            elif len(i) == 4:
                for j in di_parent_list:
                    if first_half in j:
                        poss_parent1.append(j)
                            
                for j in di_parent_list:
                    if second_half in j:
                        poss_parent2.append(j)

            elif len(i) == 6:
                for j in tri_parent_list:
                    if first_half in j:
                        poss_parent1.append(j)
                            
                for j in tri_parent_list:
                    if second_half in j:
                        poss_parent2.append(j)

            print("\nparent1: ", poss_parent1, "\nparent2: ", poss_parent2,"\n")   
        
class PredictingParents:
    def __init__(self,offspring_list):
        self.offspring_list = offspring_list

    def get_possible_parents(self):
        possible_genotypes = calculate_possible_genotypes(self.offspring_list)

        find_parents(possible_genotypes)


#While loop that allows the program to be run continuously and allows the user to access the different modules
while True:
    module = input("\n" + "Please enter 'p' to access the punnett square module. \nPlease enter 'o' to access the parent prediction module. \nPlease enter any other key to quit.\n")

    #Punnett Square module
    if module == "p":
        #While loop that allows the program to create multiple punnett squares
        while True:

            #Welcome message
            print("\n" + "Welcome to the punnett square module. You can input a monohybrid (ex: 'Aa'), dihybrid (ex: 'AaBb'), or trihybrid ('AaBbCc') cross. \nYou can choose any letters to represent your alleles, but please follow the formatting examples above. \nAlso, please note that capital letters will be calculated as 'dominant' and lowercase letters will be calculated as 'recessive'. Please enjoy!")

            #Parent genotypes
            parent1 = input("\n" + "Please enter parent 1 genotype: ")
            parent2 = input("Please enter parent 2 genotype: ")

            #Assigning the parental genotypes to the class
            cross = PredictingOffspring(parent1,parent2)

            #Running the punnett method within the class
            cross.punnett()

            #Allowing the user to create multiple punnett squares
            action = input("\n" + "Enter 'C' to make another punnett square or any key to quit" + "\n")

            if action == "C":
                print("")

            else:
                quit()

    #Parent prediction module
    if module == "o":
        #While loop that allows the program to assess multiple offspring lists
        while True:

            #Welcome message
            print("\n" + "Welcome to the parent prediction module. You can input a list of offspring phenotypes and the program will provide all possible parental genotypes for each offspring.")
            print("You can enter as many (or as few) offspring as you would like, but please adhere to the following guidelines." + "\n")
            print("Please only enter offspring for one type of cross (i.e., mono, di, or trihybrid cross). Additionally, please denote dominant phenotype(s) as 'd' and recessive phenotype(s) as 'r' and separate all entries by a space" + "\n")
            print("For example: to denote 3 dominant and 1 recessive children in a monohybrid cross, please input:  d d d r  (order does not matter).")
            print("Similarly for a dihybrid cross, 1 dominant/recessive and one recessive/dominant child would be input as:  dr rd (order does not matter).")
            print("Finally, for a trihybrid cross, 1 dominant/dominant/recessive and 1 recessive/recessive/dominant child would be input as:  ddr rrd (order does not matter). Enjoy!")
            
            #Offspring list
            user_offspring_list = input("\nPlease enter offspring list (separated by spaces): ")
            user_offspring_list_final = user_offspring_list.split()
        
            #Assigning the offspring list to a class
            prediction = PredictingParents(user_offspring_list_final)

            #Running the prediction method within the class
            prediction.get_possible_parents()

            #Allowing the user to input multiple offspring lists
            action = input("\n" + "Enter 'C' to input another offspring list" + "\n")

            if action == "C":
                print("")

            else:
                quit()
    else:
        quit()
