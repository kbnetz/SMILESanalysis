#!/usr/bin/env python
# coding: utf-8

# In[1]:


# This code takes a SMILES string and returns a few pieces of useful information:
# Molecular formula
# Molecular mass
# Degrees of unsaturation
# Number of rings
# And a few bits of information about functional groups.
# Due to the complexities of SMILES strings and the fact that there are multiple ways to encode the same molecule,
# the scope of this code is limited. However, it could be expanded upon to give more information.
# Rdkit and pysmiles are other ways of dealing with SMILES strings but this is a from-scratch code to demonstrate 
# knowledge of dictionaries and modifications of strings / lists. 
# Here are some example SMILES strings to input:
# cyclohexane = C1CCCCC1
# ethanol = CCO
# water = O or [OH2]
# aspirin = O=C(C)Oc1ccccc1C(=O)O
# chlorobenzene = Clc1ccccc1
# bromochlorodifluoromethane = C(F)(F)(Cl)Br
# penicillin G = CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C
# biphenyl = c1ccc(cc1)-c1ccccc1
# ibuprofen = CC(C)CC1=CC=C(C=C1)C(C)C(=O)O
# cyclopropane = C1CC1

# Step 1. Ask the user for a SMILES string as input.
# Step 1.1 Make sure the SMILES string is made of characters we can recognize and won't break the program. 
# Step 1.2 Convert the SMILES string into a list and simplify some characters.

accepted_characters = ['C', 'c', 'H', 'N', 'n', 'O', 'o', 'B', 'S', 's', 'P', 'F', 'Cl', 'Br', '-', '.', '=', '#', '$', '@', '(', ')', '[', ']', '/', '\\', '1', '2', '3', '4', '5', '6', '7', '8', '9', '0']
SMILES = input('Please enter a SMILES string to analyze: ')

def check_SMILES(SMILES):
    SMILES_list = list(SMILES)
    SMILES_list.append('-')
    accepted_numbers = 0
    for index_i in range(len(SMILES_list) - 1):
        try:
            index_i_int = int(SMILES_list[index_i])
            index_i1_int = int(SMILES_list[index_i + 1])
            if SMILES_list[index_i].isdigit() and SMILES_list[index_i + 1].isdigit():
                accepted_numbers += 1
            else:
                pass
        except:
            pass 
    accepted_elements = 0
    for character in SMILES_list:
        if character not in accepted_characters:
            accepted_elements += 1
        else:
            pass 
    if accepted_numbers == 0 and accepted_elements == 0:
        SMILES_ok = True
    else:
        SMILES_ok = False
    return SMILES_ok

SMILES_ok = check_SMILES(SMILES)
    
while SMILES_ok == False:
    print('This program has limited functionality and can only accept SMILES strings with accepted characters.')
    print('Accepted elements are currently C, H, N, O, B, S, P, F, Cl, and Br. Please enter a SMILES string only consisting of those elements')
    print('Additionally, this program cannot handle SMILES strings containing double-digit numbers.')
    print('')
    SMILES = input('Please enter a SMILES string to analyze: ')
    SMILES_ok = check_SMILES(SMILES)

print('')
SMILES_ = SMILES.replace('-', '')
SMILES_ = SMILES_.replace('@', '')
SMILES_ = SMILES_.replace('\\', '')
SMILES_ = SMILES_.replace('/', '')
SMILES_ = SMILES_.replace('.', '')
SMILES_ = SMILES_.replace('[', '') # square brackets are used for stereochemistry and we are not building that in currently
SMILES_ = SMILES_.replace(']', '') # square brackets
SMILES_list = list(SMILES_)
#removes hyphens, @, stereochemistry symbols
'''print(f'Analyzing {''.join(SMILES_list)}')''' # written for checking logic

# Step 2. We need to write a function that will translate any aromatic rings written in lowercase to their explicit uppercase notations with alternating double bonds. For example, benzene can be written as (CC1=CC=CC=C1) or it can be written as c1ccccc1. 
# But to make things more complicated, pyridine can be c1ccncc1 or (C1=CC=NC=C1) or n1ccccc or (N1=CC=CC=C1) - the location of the nitrogen in the ring can differ. 
# For this purpose we will write a translating function that replaces a lowercase 'x' with its uppercase version 'X', and then insert alternating '=' signs.
# Step 2.1 We first want to extract where the lowercase letters are located and save that as a list of indexes. We will use the function get_aromatic_indexes to do this.
# Step 2.2 We then need to use that list of indexes to carefully replace those sequences with the "translated" version that has the explicit bonds written in. We will use the function translate_aromatics to do this.

def get_aromatic_indexes(SMILES_list): 
    aromatic_indexes = []
    for index, character in enumerate(SMILES_list):
        if character.islower() and character != 'l' and character != 'r':
            aromatic_indexes.append(index)
    return aromatic_indexes

def translate_aromatics(SMILES_list, aromatic_indexes): # this will travel to the index at which the first lowercase letter is found, then translate it to uppercase. There should be a number following this to signify the beginning of a ring. 
    # So we don't want to change that number, we need to insert an '=' sign after that number, then translate the other lowercase letters to uppercase while inserting an '=' sign in alternating places.
    translated_SMILES_aromatics = []
    for index in range(len(SMILES_list)):
        if index not in aromatic_indexes:
            translated_SMILES_aromatics.append(SMILES_list[index]) # Copies the regular symbol if it's not part of an aromatic group. This will preserve Cl's, etc. 
        else:
            translated_SMILES_aromatics.append(SMILES_list[index].upper()) # Copies the uppercase symbol if it's part of the aromatic group. This means we're only changing the case of aromatic carbons/nitrogens/oxygens, not Cl, Br, etc. 
    # Now we need to split the list of aromatic indexes into sets of 6. This will handle if we have multiple aromatic rings. We will deal with each one separately.
    aromatic_index_sets = []
    for i in range(0, len(aromatic_indexes), 6):
        aromatic_index_sets.append(aromatic_indexes[i:i+6])
    '''print(f'The groups of 6 for the aromatic sets are {aromatic_index_sets}')''' # written for checking logic
    # for aromatic_set in aromatic_index_sets: # For each set, we want to change the translated_SMILES_aromatics to include the appropriate equals signs.
    # Each "set" in our aromatic_index_sets list of list corresponds to a string x#xxxxx# will translate to x#=xx=xx=x#. Our index_set for this will be [0,2,3,4,5,6] and we want to add an = at 2, then at 4, then at 6. 
    # But sometimes there are parentheses, other characters, etc so we need to follow the directions given to us by the set of aromatic indexes.
    # For example, biphenyl is encoded by the SMILES string "c1ccc(cc1)-c1ccccc1", which gives you the following groups for the aromatic sets: [[0, 2, 3, 4, 6, 7], [11, 13, 14, 15, 16, 17]]"
    # We want to add an equals sign right before index 2, 4, and 7, and then right before 13, 15, and 17. You can see that the spacing here isn't constant so we can't just iterate by 2. 
    # get values to add an '=' sign at this spot.
    add_bond_spots = []
    for aromatic_set in aromatic_index_sets:
        add_bond_spots.append(aromatic_set[1])
        add_bond_spots.append(aromatic_set[3])
        add_bond_spots.append(aromatic_set[5])
    '''print(f'We need to add the = sign at these positions: {add_bond_spots}')''' # written for checking logic
    # But, each time we add an '=' sign, the string gets longer by 1. So we need to add '=' at index 2, index 4 + 1, index 7 + 2, index 13 + 3, index 15 + 4, index 17 + 5, etc. We will create a new list with the additions.
    add_bond_spots_adjusted = add_bond_spots
    for i in range(len(add_bond_spots_adjusted)):
        add_bond_spots_adjusted[i] += i
    '''print(f'The adjusted places to add the = sign are {add_bond_spots_adjusted}')''' # written for checking logic
    # Now, we can add the '=' sign at the indicated spots from our list.
    for index in add_bond_spots_adjusted:
        translated_SMILES_aromatics.insert(index, '=')
    return translated_SMILES_aromatics
    
aromatic_indexes = get_aromatic_indexes(SMILES_list)
'''print(f'The indexes corresponding to aromatic groups indicated by lowercase are {aromatic_indexes}')''' # written for checking logic
translated_SMILES_aromatics = translate_aromatics(SMILES_list, aromatic_indexes)
'''print(f'The translated SMILES string where lowercase aromatics are denoted by uppercase letters with explicit bonds is {''.join(translated_SMILES_aromatics)}')''' # written for checking logic

# Step 3. Now that we don't have to deal with any lowercase aromatic rings, we can further convert the SMILES string into a list with a few changes in order to enable processing of this list:
# Step 3.1 Convert multi-lettered atom symbols into placeholders. This will make sure that "Cl" isn't interpreted as "C" + "l". "Cl" will map to "L" and "Br" will map to "R". Save this as "abbreviated SMILES"
# Step 2.3 Remove any explicitly defined hydrogens. Most SMILES strings omit implicit hydrogens, but sometimes they can be written in. We will have to later add in the correct number of H's based on the degree of unsaturation present in the molecule. 

def abbreviate_SMILES(translated_SMILES_aromatics):
    translated_SMILES_aromatics.append(' ')
    abbrev_SMILES = []
    for index in range(len(translated_SMILES_aromatics)):
        if translated_SMILES_aromatics[index] == 'C':
            if translated_SMILES_aromatics[index + 1] == 'l':
                abbrev_SMILES.append('L')
            else:
                abbrev_SMILES.append('C')
        elif translated_SMILES_aromatics[index] == 'B':
            if translated_SMILES_aromatics[index + 1] == 'r':
                abbrev_SMILES.append('R')
            else:
                abbrev_SMILES.append('B')
        elif translated_SMILES_aromatics[index] != 'l' and translated_SMILES_aromatics[index] != 'r' and translated_SMILES_aromatics[index] != 'H' and translated_SMILES_aromatics[index] != ' ':
            abbrev_SMILES.append(translated_SMILES_aromatics[index])
    return abbrev_SMILES

abbrev_SMILES = abbreviate_SMILES(translated_SMILES_aromatics)

print(f'The SMILES string with standardized aromatic rings and stereochemistry removed is {''.join(abbrev_SMILES)}')
print('Calculating some properties about this molecule...')
print('')

# Now we have a list of all the elements, besides hydrogen, where lowercase and uppercase letters have been neatened up, and multi-letter atom symbols have been converted into single letter atom symbols.
# This abbreviated list still contains symbols that indicate order of bonding, stereochemistry, or aromatic bonds, if written out. However, it has lost the lowercase signifier of aromatic rings, so we will have to pay attention to that later.

# Step 3. Calculate the degrees of unsaturation present in the molecule.
# Step 3.1 Calculate the number of rings. Each ring adds 1 to our unsaturation count. This requires us to look for pairs of numbers that indicate connectivity.
# Since we already converted any lowercase-encoded aromatic rings into uppercase and explicit unsaturation, we won't have to deal with that exception!
# First, we will have to look through our list of characters, which are all in string form, and see if we can convert them to integers.
# Step 3.2 Calculate the number of double, triple, and quadruple bonds. Each increased order of bond adds 1 to our unsaturation count. 

def make_ring_list(abbrev_SMILES):
    SMILES_ring_list = []
    for index in range(len(abbrev_SMILES)):
        try:
            number_index = int(abbrev_SMILES[index])
            SMILES_ring_list.append(number_index)
        except:
            pass
    return SMILES_ring_list
    
def count_rings(SMILES_ring_list):
    ring_count = 0
    SMILES_ring_list.sort() #the order of the ring connectivity doesn't matter for this calculation, all that matters is that the numbers are found in pairs.
    for index in range(0, len(SMILES_ring_list) - 1):
        if SMILES_ring_list[index] == SMILES_ring_list[index + 1]:
            ring_count +=1
    return ring_count

SMILES_ring_list = make_ring_list(abbrev_SMILES)
'''print(f'The numbers indicating rings in your SMILES string are {SMILES_ring_list}')''' # written for checking logic

ring_count = count_rings(SMILES_ring_list)
'''print(f'The total number of rings in your SMILES string is {ring_count}')''' # written for checking logic

def count_bond_unsaturation(abbrev_SMILES):
    bond_unsaturation = 0
    for index in range(len(abbrev_SMILES)):
        if abbrev_SMILES[index] == '=':
            bond_unsaturation +=1
        elif abbrev_SMILES[index] == '#':
            bond_unsaturation +=2
        elif abbrev_SMILES[index] == '$':
            bond_unsaturation +=3
    return bond_unsaturation

bond_unsaturation = count_bond_unsaturation(abbrev_SMILES)
'''print(f'The degree of unsaturation corresponding to bonds of multiple order is {bond_unsaturation}')''' # written for checking logic

degrees_of_unsaturation = ring_count + bond_unsaturation 
print(f'The total degree of unsaturation in the molecule is {degrees_of_unsaturation}')

# Step 4. Now we need to set up to calculate the molecular formula. 
# Step 4.1 We will first have to count how many of each atom appear, excepting hydrogen.
# Step 4.2 Then, we will use the degrees of unsaturation to calculate how many hydrogens should be on the molecule.
# Step 4.3 Then, we will display the molecular formula.

def create_atom_dictionary(abbrev_SMILES):
    atom_counts = {'C': 0, 'N': 0, 'O': 0, 'S': 0, 'P': 0, 'B': 0, 'F': 0, 'L': 0, 'R': 0}
    for character in abbrev_SMILES:
        if character.isupper(): # this will only choose the atoms, not any parentheses or symbols for bond order
            atom_counts[character] +=1
    return atom_counts

atom_counts = create_atom_dictionary(abbrev_SMILES)
'''print(f'The total counts of each atom in your SMILES string is {atom_counts}')''' # written for checking logic

def calculate_hydrogens(atom_counts, degrees_of_unsaturation):
    H_count = -(2 * degrees_of_unsaturation) + (2 * int(atom_counts['C'])) + 2 + int(atom_counts['N']) + int(atom_counts['B']) + int(atom_counts['P']) - int(atom_counts['F']) - int(atom_counts['L']) - int(atom_counts['R'])
    return H_count

H_count = calculate_hydrogens(atom_counts, degrees_of_unsaturation)
'''print(f'The total number of hydrogens in your molecule should be {H_count}')''' # written for checking logic

def get_molecular_formula(atom_counts, H_count):
    formula_list = []
    for character in atom_counts:
        if atom_counts[character] != 0:
            if character != 'R' and character != 'L':
                formula_list.append(str(character) + str(atom_counts[character]))
            if character == 'R':
                formula_list.append('Br' + str(atom_counts[character]))
            if character == 'L':
                formula_list.append('Cl' + str(atom_counts[character]))
    formula_list.append('H' + str(H_count))
    molecular_formula = ''.join(formula_list)
    return molecular_formula
            
molecular_formula = get_molecular_formula(atom_counts, H_count)
print("The molecular formula of your molecule is %s" % molecular_formula)

# Step 5. Now we will calculate the molecular weight!
# Step 5.1 Create a dictionary with the atomic weights of all of the elements included in our SMILES parser. This will be constant
# Step 5.2 Create a function that multiplies atom_counts by ATOMIC_WEIGHTS

ATOMIC_WEIGHTS = {'C': 12.011, 'N': 14.007, 'O': 15.999, 'S': 32.066, 'P': 30.974, 'B': 10.811, 'F': 18.998, 'L': 35.453, 'R': 79.904, 'H': 1.008}

def get_molecular_weight(atom_counts, ATOMIC_WEIGHTS):
    atom_counts['H'] = H_count
    MW = 0
    for atom in atom_counts:
        MW += atom_counts[atom] * ATOMIC_WEIGHTS[atom]
    MW = round(MW, 3)
    return MW

MW = get_molecular_weight(atom_counts, ATOMIC_WEIGHTS)
print(f'The molecular weight of your molecule is {MW}')

# Step 6. We want to present some information about the functional groups present in the molecule. This can be difficult given the different orders, parentheses, numbers, etc present in SMILES strings, and how the same molecule can be written in different ways. 
# Step 6.1 We need to create a list that represents our SMILES string that removes any numbers. This will clear up the formatting so we can find functional groups. It removes any information regarding rings. 
# We have already removed information regarding stereochemistry including square brackets.
# Step 6.2 We will now look through this cleaned SMILES for functional groups that the program can recognize.
# Here is a list of a number of functional groups and their SMILES maps: Note that there could be numbers interspersed between these characters if these are part of a ring. So we will have to remove numbers if we want to identify functional groups.
# This is not comprehensive and due to the limitations of SMILES and not converting them into a graph representation, we cannot handle every case especially if things are written in different orders
# Note that if a functional group bridges a ring junction this information is lost and this functional group will not be captured. For example N2C(S)(C)(C2=O) has an amide where N2 is linked to C2. But we lose that information when we remove numbers. We will have to figure this out later
# carboxylic acid = C(=O)O or OC(=O)
# ester = CC(OC)=O
# amide = CC(N)=O
# acetone (ketone) = CC(C)=O
# acetaldehyde = CC([H])=O or C=O without the other qualifiers that make things an ester, amide, ketone, or carboxylic acid. 
# cyano = C#N
# trifluoromethyl = C(F)(F)F
# difluoromethyl = C(F)F
# sulfonyl = S(=O)(=O)
# nitromethane = C[N+]([O-])=O
# phosphate = CP(O)(O)=O
# ether = COC

def cleaned_SMILES(abbrev_SMILES):
    cleaned_SMILES = []
    for character in abbrev_SMILES:
        if '0' <= character <= '9':
            pass
        else:
            cleaned_SMILES.append(character)
    return cleaned_SMILES

def find_functional_groups(cleaned_SMILES):
    for i in range(13):
        cleaned_SMILES.append(' ') # this just adds empty space on to the end of our cleaned_SMILES so that we can check out to the end of our nitromethane, which is the longest specified functional group.    
    if 'C#N' in ''.join(cleaned_SMILES):
        print('Your compound has a nitrile / cyano group!')
    if 'C(=O)O' in ''.join(cleaned_SMILES) or 'OC(=O)' in ''.join(cleaned_SMILES):
        print('Your compound has a carboxylic acid or ester group. Sorry I am not able to distinguish these yet! If I do not print that there is an ester after this, then it would represent a carboxylic acid.')
    if 'C(OC)=O' in ''.join(cleaned_SMILES):
        print('Your compound has an ester in it!')
    if 'C(N)=O' in ''.join(cleaned_SMILES) or 'NC(=O)' in ''.join(cleaned_SMILES):
        print('Your compound has an amide in it!')
    else:
        print('There are no other functional groups that I can identify right now')

cleaned_SMILES = cleaned_SMILES(abbrev_SMILES)
'''print(f'Your SMILES string with the numbers removed is {''.join(cleaned_SMILES)}')''' # written for checking logic
print('')
print('I can identify some basic functional groups, but this code is not comprehensive and I will not be able to find everything.')
print('Here are the ones that this program recognizes currently:')
find_functional_groups(cleaned_SMILES)

# Step 7. Representing molecules as lists loses / obfuscates the 2- and 3-dimensional connectivity in a molecule. 
# It is better to represent a molecule as an adjacency matrix. Maybe later I could make atom classes and save information, but for now we will try to make adjacency matrices.
# 7.1 Take the SMILES list where we have converted it to uppercase, explicitly defined aromatic rings, and replaced Cl and Br with L and R, respectively. This is the list abbrev_SMILES
# Use ibuprofen to test. SMILES = CC(C)CC1=CC=C(C=C1)C(C)C(=O)O

# 7.2 Write a function that makes an empty adjacency matrix.
# First, we want the indexes that correspond to letters, since this will be the atoms. All other characters determine bonding and connectivity.
# Then, we will create a matrix that has '0' in each spot with length and width equal to the length of our SMILES string. This will include special symbols such as bonds, etc, which will always end up being rows or columns of 0s.
# We will remove the all-zero rows and columns later, but this will help keep track of indexes.

def get_atom_indexes(abbrev_SMILES):
    atom_indexes = []
    for index, character in enumerate(abbrev_SMILES):
        if character.isupper():
            atom_indexes.append(index)
    return atom_indexes

def create_empty_matrix(abbrev_SMILES):
    empty_matrix = [[0 for _ in range(len(abbrev_SMILES))] for _ in range(len(abbrev_SMILES))]
    return empty_matrix

atom_indexes = get_atom_indexes(abbrev_SMILES)
# print(f'There are atoms at the following locations of the SMILES string: {atom_indexes}')

adjacency_matrix = create_empty_matrix(abbrev_SMILES)

def modify_empty_matrix(adjacency_matrix, abbrev_SMILES, atom_indexes):
    modified_matrix = [[0 for _ in range(len(abbrev_SMILES))] for _ in range(len(abbrev_SMILES))]
    for i in range(len(abbrev_SMILES)):
        for j in range(len(abbrev_SMILES)):
            if i in atom_indexes and j in atom_indexes:
                modified_matrix[i][j] = adjacency_matrix[i][j]
            else:
                modified_matrix[i][j] = '.' # if the character at that spot doesn't correspond to an atom, it will be denoted by '.'
    return modified_matrix
    
adjacency_matrix = modify_empty_matrix(adjacency_matrix, abbrev_SMILES, atom_indexes)

'''print('Here is your empty matrix:')
for row in adjacency_matrix:
    print(row)''' # added for testing purposes
'''adjacency_matrix[0][1] = 1
print('Testing:')
for row in adjacency_matrix:
    print(row)''' # added for testing purposes

# 7.3 Now we need to determine criteria for how things are bonded
# Function for determining ring connectivity: we want to know at which indexes do numbers indicating rings appear. So we will look along the abbrev_SMILES list and pull out the positions at which numbers occur. 

def get_ring_indexes(abbrev_SMILES):
    ring_indexes = []
    for index, character in enumerate(abbrev_SMILES):
        try: 
            intcheck = int(character)
            ring_indexes.append(index)
        except:
            pass
    return ring_indexes

ring_indexes = get_ring_indexes(abbrev_SMILES)
# print(f'There are rings denoted at the following locations of the SMILES string: {ring_indexes}')

# For molecules with multiple rings, we need to match which position is connected to which. We have a list of indexes that correspond to where numbers appear. These numbers refer to the character directly preceding them.
# Let's create a dictionary of ring pairs. This will show which indexes are matched together.
def get_ring_pairs(abbrev_SMILES):
    ring_indexes = get_ring_indexes(abbrev_SMILES)
    ring_pairs = []
    for k in range(len(abbrev_SMILES)):
        for l in range(len(abbrev_SMILES)):
            if k in ring_indexes and l in ring_indexes and k != l:
                if abbrev_SMILES[k] == abbrev_SMILES[l]:
                    ring_pairs.append([k,l])
                    ring_indexes.remove(k) # removing k and l prevents double-counting of ring pairs in this function
                    ring_indexes.remove(l)
    return ring_pairs

ring_pairs = get_ring_pairs(abbrev_SMILES)
'''print(f'The pairs of indexes corresponding to rings are: {ring_pairs}')''' # for tracking purposes

# 7.4 Now, let's add these rings to the adjacency matrix. Rings denoted by numbers will always be joined by single bonds, so we don't have to worry about bond order here. We will add a '1' at the spots indicated by the ring_pairs.
# Complicating this is the fact that the ring pairs indicate where the number is written, not the atom that is actually linked in the SMILES string. We will add in a "-1" to add the adjacency to the atom itself.
# So for a ring pair of [a, b], row_a = empty_matrix[a] and column_b = row_a[b]. This will be set equal to 1. We will also need to set row_b = empty_matrix[b] and column_a = row_b[a] = 1
# A ring pair is denoted by ring_pairs[i]
# troubleshoot with cyclopropane = C1CC1

def adjacency_matrix_add_rings(adjacency_matrix, ring_pairs):
    for i in range(len(ring_pairs)): # iterates on each ring pair item in the list of ring pairs
        ring_pair = ring_pairs[i]
        ring_index_a = ring_pair[0]
        ring_index_b = ring_pair[1]
        adjacency_matrix[ring_index_a - 1][ring_index_b - 1] += 1
        adjacency_matrix[ring_index_b - 1][ring_index_a - 1] += 1
    return adjacency_matrix

adjacency_matrix = adjacency_matrix_add_rings(adjacency_matrix, ring_pairs)
'''print('adjacency matrix with rings added:')
for row in adjacency_matrix:
    print(row)''' # for testing purposes

# 7.5 Now we need to address the complicated system of bonding. There are a number of cases:
# Case 1: The simplest: no bonds specified. Examples are CC, CO, CN, etc. There are no parentheses or special characters. In this case, index(atom 1) + 1 = index(atom 2). We will add an adjacency of "1" to the spots in our matrix because this is a single bond.

def add_simple_bonds(adjacency_matrix, abbrev_SMILES):
    for i in range(len(abbrev_SMILES) - 1): #atom_indexes:
        if abbrev_SMILES[i].isupper():
            if abbrev_SMILES[i + 1].isupper(): # indicates a single bond to atom(i+1). 
                adjacency_matrix[i][i + 1] += 1
                adjacency_matrix[i + 1][i] += 1
            if abbrev_SMILES[i + 1] == '=': # indicates a double bond to the next atom at (i+2)
                adjacency_matrix[i][i + 2] += 2
                adjacency_matrix[i + 2][i] += 2
            if abbrev_SMILES[i + 1] == '#': # indicates a triple bond to the next atom at (i+2)
                adjacency_matrix[i][i + 2] += 3
                adjacency_matrix[i + 2][i] += 3
            if abbrev_SMILES[i + 1] == '$': # indicates a quadruple bond to the next atom at (i+2)
                adjacency_matrix[i][i + 2] += 4
                adjacency_matrix[i + 2][i] += 4
            try:
                if int(abbrev_SMILES[i + 1]) >= 0: # if the next character is a number, then we need to repeat the above if statements just adding 1 to the index
                    # this covers cases if there is a ring number denoted right after our atom at index i
                    # note that this code will break if there are double-digit ring numbers. That's not usually an issue for small organic molecules though
                    if abbrev_SMILES[i + 2].isupper(): # indicates a single bond to atom(i+2). - Note how indexes increased again by 1
                        adjacency_matrix[i][i + 2] += 1
                        adjacency_matrix[i + 2][i] += 1
                    if abbrev_SMILES[i + 2] == '=': # indicates a double bond to the next atom at (i+3)
                        adjacency_matrix[i][i + 3] += 2
                        adjacency_matrix[i + 3][i] += 2
                    if abbrev_SMILES[i + 2] == '#': # indicates a triple bond to the next atom at (i+3)
                        adjacency_matrix[i][i + 3] += 3
                        adjacency_matrix[i + 3][i] += 3
                    if abbrev_SMILES[i + 2] == '$': # indicates a quadruple bond to the next atom at (i+3)
                        adjacency_matrix[i][i + 3] += 4
                        adjacency_matrix[i + 3][i] += 4 
            except:
                pass
    return adjacency_matrix

adjacency_matrix = add_simple_bonds(adjacency_matrix, abbrev_SMILES)
'''for row in adjacency_matrix:
    print(row)''' # for testing purposes

# 7.6 Case 2: Branching. This is specified by a parentheses. 
# Use this SMILES to troubleshoot: C(F)(F)(Cl)Br
# This means that the first atom is bonded to all of the following atoms that are specified by parentheses, and it is also bonded to the last atom, which is not specified in parentheses.
# Atom at index i: If abbrev_SMILES[i + 1] == (, then atom at index i is bonded to atom at index i + 2, if there is an atom. However, there might be more parentheses specifying more branching. Example CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C
# There are also brackets. Brackets denote stereochemistry so we will replace them for now with round parentheses. 
# Define a function called find_parentheses: it will create a list of indexes where parentheses are found.

def translate_brackets(abbrev_SMILES):
    for index, character in enumerate(abbrev_SMILES):
        if character == '[':
            abbrev_SMILES[index] = '('
        elif character == ']':
            abbrev_SMILES[index] = ')'
    return abbrev_SMILES

abbrev_SMILES = translate_brackets(abbrev_SMILES)
# print(abbrev_SMILES)

def find_nested_parentheses(abbrev_SMILES):
    nested_parentheses_map = {}
    temporary_stack = []
    for index, character in enumerate(abbrev_SMILES):
        if character == '(':
            temporary_stack.append(index) # adds this index to the stack as we look for the next closing parentheses
        elif character == ')':
            try:
                matched_parentheses_index = temporary_stack.pop() # removes the previous index of '(' from the temproary stack
                nested_parentheses_map[matched_parentheses_index] = index # maps the opening parentheses index to the closed parentheses index we just found
            except:
                pass
    return nested_parentheses_map

nested_parentheses_map = find_nested_parentheses(abbrev_SMILES)
# print(nested_parentheses_map)

# 7.7 Now that we have a map of which parentheses are matched, then we can start tracing them and pulling out the atoms that are specified in each parentheses block. 
# An atom is bonded to all following atoms that start with a parentheses. So any atom that is '(X' after atom Y is bonded to Y. For example, triethylamine would be written as N(CC)(CC)CC and the nitrogen is bonded to the three carbons right after the open parentheses. 
# The exception is the last atom - it is not present in the parentheses, but the nitrogen is still bonded to it. 

def find_branched_atoms(abbrev_SMILES, atom_indexes, nested_parentheses_map, ring_indexes):
    branched_atoms = []
    for index in range(len(abbrev_SMILES)):
        if index in atom_indexes: # selects indexes that indicate an atom
            if index + 1 in nested_parentheses_map.keys(): # this means that the next character corresponds to an open bracket
                branched_atoms.append(index)
            elif index + 1 in ring_indexes:
                if index + 2 in nested_parentheses_map.keys(): # this is the case when branching directly follows an atom that has a number indicating it is part of a ring junction
                    branched_atoms.append(index)
    return branched_atoms

branched_atoms = find_branched_atoms(abbrev_SMILES, atom_indexes, nested_parentheses_map, ring_indexes)
# print(f'The atoms at the following position have been flagged as branched atoms: {branched_atoms}')

def simplify_nested_parentheses(abbrev_SMILES, nested_parentheses_map, atom_index): # We will give the function an atom index. Then, any nested parentheses after that will be collapsed.
# Example: CC(CC(CC))C should be converted to CC(C)C and we record that we removed 5 characters: we removed between the *'s in CC(C*C(CC)*)C
# We are given an atom index. In this case, we have index 1. This will come from our branched_atoms list later. 
# Our nested set of parentheses in this case is {5: 8, 2: 9}
# Look at nested_parentheses_map.keys() and select keys for which key > atom_index. Delete any other parentheses. 
# This will prevent us from accidentally collapsing what we are looking at if our atom index is part of a nested set of parentheses. 
# Then, look at the remaining parentheses. If one is contained within the other, we will delete that. Use CC(CC)CC(CC(CC))C with atom index = 7 to test.
    characters_removed = 0
    truncated_SMILES = abbrev_SMILES[atom_index:]
    nested_parentheses_modified = find_nested_parentheses(truncated_SMILES)
    characters_removed += len(abbrev_SMILES[:atom_index])
    #print(f'The set of open and closed parentheses remaining in the truncated list is {nested_parentheses_modified}')
    #print(f'The truncated SMILES string is {truncated_SMILES}')
    #print(f'Removed {characters_removed} characters from the beginning of the SMILES string') # This index has to be added to ALL following numbers when they are pulled out. 
    #return truncated_SMILES Test for debugging
    ''' def collapse_parentheses(truncated_SMILES, nested_parentheses_modified): Test for debugging individual functions'''
    # Collapses nested parentheses so that only the outer parentheses is kept
    # For key1, value1 in nested_parentheses_modified: if there is key2 < key1 and value2 > value1, we will delete the string contained within key1, value1
    # Look at all nested_parentheses_modified.keys() and find a key2 < key1
    parentheses_listed = []
    for key in nested_parentheses_modified.keys():
        key_listed = [key, nested_parentheses_modified[key]]
        parentheses_listed.append(key_listed)
    # print(parentheses_listed)
    removed_parentheses = []
    for i in range(len(parentheses_listed)):
        pair1 = parentheses_listed[i]
        for n in range(len(parentheses_listed)):
            pair2 = parentheses_listed[n]
            if pair1[0] > pair2[0] and pair1[1] < pair2[1]:
                removed_parentheses.append(pair1) # This is a list of the parentheses that are part of nested groups. So now we will remove them from listed parentheses.
    parentheses_removed_dict = {}
    for pair in removed_parentheses:
        parentheses_removed_dict[pair[0]] = pair[1]
    # print(parentheses_removed_dict)
    cleaned_parentheses = {}
    for key in nested_parentheses_modified:
        if key not in parentheses_removed_dict:
            cleaned_parentheses[key] = nested_parentheses_modified[key]
        else:
            pass
    # print(cleaned_parentheses)
    # Now we just need to correct back for the indexes we removed from the beginning of the string before atom i. We can add them back in now, characters_removed
    almost_simplified_nested_parentheses = {}
    for key in cleaned_parentheses:
        corrected_key = key + characters_removed
        almost_simplified_nested_parentheses[corrected_key] = cleaned_parentheses[key] + characters_removed
    # print(f'The almost simplified parentheses are {almost_simplified_nested_parentheses}') 
    # Another correction: if there are multiple nested parentheses following an atom, but aren't connected to that atom, they will still show up.
    # This is evident when the simplified nested parentheses are not sequential: ex. {2: 5, 8: 15}. The pair encoded by {8, 15} is not connected to atom at index 1.
    # Remove any pair2 if key2 - value1 > 0. 
    '''This needs to be fixed because it will remove multiple branching! Ex. C(C)(C)(C). We will add this back in later.'''
    cleaned_parentheses_list = []
    for key in almost_simplified_nested_parentheses:
        key_listed = [key, almost_simplified_nested_parentheses[key]]
        cleaned_parentheses_list.append(key_listed)
    #nonsequential_parentheses = {}
    linked_parentheses = [cleaned_parentheses_list[0]]
    for i in range(len(cleaned_parentheses_list) - 1):
        selected_pair = cleaned_parentheses_list[i]
        next_pair = cleaned_parentheses_list[i + 1]
        if selected_pair[1] + 1 == next_pair[0]:
            linked_parentheses.append(next_pair)
        else:
            pass
    linked_parentheses_dict = {}
    for i in range(len(linked_parentheses)):
        pair = linked_parentheses[i]
        key = pair[0]
        value = pair[1]
        linked_parentheses_dict[key] = value
    # print(f'The finalized simplified parentheses corresponding to the atom at {atom_index} are {linked_parentheses_dict}')
    return linked_parentheses_dict

# Step 7.8 We can add these identified branching atoms as single bonds. Later, we will figure out how to handle the case of double and triple bonds. 

def matrix_add_branched_atoms(abbrev_SMILES, nested_parentheses_map, branched_atoms):
    for atom in branched_atoms:
        # print(f'Analyzing atom at position {atom}')
        linked_parentheses_dict = simplify_nested_parentheses(abbrev_SMILES, nested_parentheses_map, atom)
        # Now we have a final dictionary removing nested parentheses given a certain atom index that gives the sets of parentheses indicating bonding to an atom at that position.
        # At this point we can retrieve the atoms directly bonded to an atom at a particular index.
        # An atom at index i with simplified nested parentheses following i will be bonded to the first letter following the parentheses
        # So if the nested parentheses cleaned up dict = {k : l, i : j}, then atom i will be bonded to the atoms at index k + 1 and index i + 1. 
        # If the value of the final pair, j, has j + 1 = an atom, then atom i is bonded to the atom at l + 1.
        atom_index_a = atom
        linked_parentheses_keys = linked_parentheses_dict.keys()
        for key in linked_parentheses_keys:
            if 'A' < abbrev_SMILES[key + 1] < 'Z': # If the character at position key + 1 in the string is an atom, denoted by a letter:
                atom_index_b = key + 1
                try:
                    adjacency_matrix[atom_index_a][atom_index_b] += 1
                    adjacency_matrix[atom_index_b][atom_index_a] += 1
                except:
                    pass
            elif abbrev_SMILES[key + 1] == '=':
                atom_index_b = key + 2 # Indicates a double bond to the following atom
                try:
                    adjacency_matrix[atom_index_a][atom_index_b] += 2
                    adjacency_matrix[atom_index_b][atom_index_a] += 2
                except:
                    pass
            elif abbrev_SMILES[key + 1] == '#':
                atom_index_b = key + 2 # Indicates a triple bond to the following atom. 
                # Not actually possible with branching to main organic elements but if I expand to further periodic table elements, it will be necessary
                try:
                    adjacency_matrix[atom_index_a][atom_index_b] += 3
                    adjacency_matrix[atom_index_b][atom_index_a] += 3
                except:
                    pass
            else:
                pass
        linked_parentheses_values = linked_parentheses_dict.values()
        end_of_parentheses = max(linked_parentheses_values)
        try:
            if 'A' < abbrev_SMILES[end_of_parentheses + 1] < 'Z': # If the character after the last parentheses is a letter, that will be the final thing bonded to atom 1.
                atom_index_b = end_of_parentheses + 1
                adjacency_matrix[atom_index_a][atom_index_b] += 1
                adjacency_matrix[atom_index_b][atom_index_a] += 1  
        except:
            pass
    return adjacency_matrix     

adjacency_matrix = matrix_add_branched_atoms(abbrev_SMILES, nested_parentheses_map, branched_atoms)
'''for row in adjacency_matrix:
    print(row)'''

# Step 8. The last task is to print our adjacency matrix and clean it up so that the numbers correspond to bonding orders and coordinates between atoms.
# We don't need any special characters in the matrix anymore, we will just shorten it to atoms. 

def atoms_only_adjacency_matrix(adjacency_matrix): # this is only for printing because it turns 0's into spaces: ' '
    atoms_adjacency_matrix = []
    for row in adjacency_matrix:
        atoms_row = []
        for char in row:
            try:
                if 0 < int(char):
                    atoms_row.append(str(char))
                if 0 == int(char):
                    atoms_row.append(' ') # If you want to maintain the information of 0-order bonds, change this to '0' instead of '.'
            except:
                pass
        if atoms_row != []:
            atoms_adjacency_matrix.append(atoms_row)
    return atoms_adjacency_matrix
    
atoms_adjacency_matrix = atoms_only_adjacency_matrix(adjacency_matrix)

# The following functions will print the adjacency matrix in a way that has atom labels on the top and on the rows for easier viewing.

def create_atom_list(abbrev_SMILES, atom_indexes):
    atom_list = []
    for index in atom_indexes:
        atom = abbrev_SMILES[index]
        atom_list.append(atom)
    return atom_list

def SMILES_atom_labels(atom_list): # This is a function that creates a row that looks like a list printed out
    # It will align with our adjacency matrix and make it nice to look at so we can correspond atoms
    row_string = ', '.join(atom_list)
    final_label = '      [' + row_string + ']'
    return final_label

def create_row_label(atom_list, index):
    atom_index_label = str(atom_list[index])
    row_label = '[' + atom_index_label + ']   '
    return row_label

atom_list = create_atom_list(abbrev_SMILES, atom_indexes)
top_labels_aligned = SMILES_atom_labels(atom_list)

print('')
print('The bonding order and adjacency matrix corresponding to this SMILES string is:')
print('If there is no bond present between two atoms, this is indicated by an empty space in the matrix rather than a 0.')
print('')
print(top_labels_aligned)
print('')

for i in range(len(atom_list)):
    row_atom_label = create_row_label(atom_list, i)
    row = atoms_adjacency_matrix[i]
    row_string = ', '.join(row)
    final_row = '[' + row_string + ']'
    final_row_labeled = row_atom_label + final_row
    print(final_row_labeled)


# In[ ]:




