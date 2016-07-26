import numpy as np
import matplotlib.pyplot as plt
import time      # TEST
import sys   

position = "F:\\BaiduYunDownload\\test\\"
filename = position + "SO2.outmol"

############ User Options ########
atomSelection = None                     # None: analyse all atoms
#atomSelection = [1, 2]                  # Include the atoms to be analysed
valanceElectronsSelection = True         # True: only consider valance electrons NOT FULLY TESTED!
ms_mergeSelection = True                 # True: ignore the differences in ms
atom_mergeSelection = False              # True: ignore the differences in the same atom
                                         # If True, valance/ms switches are TURN OFF automatically
color = "hot"                            # "bone"/"jet"/"spectral"/...

############ Built-in Functions ###########

def sequenceGenerator():                # Generate number sequence
    output = 0
    while True:
        yield output
        output += 1
naturalNumber = sequenceGenerator()

class Electron:
    """
    Contain main character of an electron, including atom number, n, l, ms, and occupation
    """
    def __init__(self, atomNum, electron_n, electron_l, electron_ms, occupation):
        self.atomNum = atomNum
        self.electron_n = electron_n
        self.electron_l = electron_l
        self.electron_ms = electron_ms
        self.occupation = occupation
        self.num = naturalNumber.next()
    def __lt__(self, other):
        if self.get_atomNum() < other.get_atomNum():
            return True
        elif self.get_atomNum() == other.get_atomNum():
            if self.get_electron_n() < other.get_electron_n():
                return True
            elif self.get_electron_n() == other.get_electron_n():
                if self.get_electron_l() < other.get_electron_l():
                    return True
                elif self.get_electron_l() == other.get_electron_l():
                    if self.get_electron_ms() < other.get_electron_ms():
                        return True
                    elif self.get_electron_ms() == other.get_electron_ms():
                        if self.get_num() < other.get_num():
                            return True
        return False
    def __repr__(self):
        if self.electron_n == None:
            return "A-%i" %self.atomNum
        elif self.electron_l == None:
            return "A-%i: %i" %(self.atomNum, self.electron_n)
        elif self.electron_ms == None:
            return "A-%i: %i %i" %(self.atomNum, self.electron_n, self.electron_l)
        else:
            return "A-%i: %i %i %i" %(self.atomNum, self.electron_n, self.electron_l, self.electron_ms)
    def __str__(self):
        if self.electron_n == None:
            return "Atom %i" %self.atomNum
        elif self.electron_l == None:
            return "Atom %i: n = %i" %(self.atomNum, self.electron_n)
        elif self.electron_ms == None:
            return "Atom %i: n = %i, l = %i" %(self.atomNum, self.electron_n, self.electron_l)
        else:
            return "Atom %i: n = %i, l = %i, ms = %i" %(self.atomNum, self.electron_n, self.electron_l, self.electron_ms)
    
    def get_atomNum(self):
        return self.atomNum
    def get_electron_n(self):
        return self.electron_n
    def get_electron_l(self):
        return self.electron_l
    def get_electron_ms(self):
        return self.electron_ms
    def get_electron_occupation(self):
        return self.occupation
    def get_num(self):
        return self.num

    def change_atomNum(self, newnum):
        self.atomNum = newnum
    def change_electron_n(self, newnum):
        self.electron_n = newnum
    def change_electron_l(self, newnum):
        self.electron_l = newnum
    def change_electron_ms(self, newnum):
        self.electron_ms = newnum
    def change_electron_occupation(self, newnum):
        self.occupation = newnum
    def change_num(self, newnum):
        self.num = newnum

def exchangeSymMatrix(matrix, i, j):
    """
    Exchange the i'th and j'th row and column of a symmatrix matrix
    """
    matrix[[i, j], :] = matrix[[j, i], :]
    matrix[:, [i, j]] = matrix[:, [j, i]]

def sub_sort(electrons,matrix,low,high):
    key = electrons[low]
    while low < high:
        while low < high and electrons[high] > key:
            high -= 1
        while low < high and electrons[high] < key:
            electrons[low] = electrons[high]
            exchangeSymMatrix(matrix, low, high)
            low += 1
            electrons[high] = electrons[low]
            exchangeSymMatrix(matrix, low, high)
    electrons[low] = key
    return low

def quick_sort(electrons,matrix,low = None,high = None):
    sys.setrecursionlimit(100000)
    if low == None:
        low = 0
    if high == None:
        high = len(electrons) - 1
    if low < high:
        key_index = sub_sort(electrons,matrix,low,high)
        quick_sort(electrons,matrix,low,key_index)
        quick_sort(electrons,matrix,key_index+1,high)
    
def reIndex(electrons):
    """
    Re-index the electrons after sorting
    """
    naturalNumber = sequenceGenerator()
    for electron in electrons:
        electron.change_num(naturalNumber.next())

def atomFilter(electrons, matrix, selection):
    """
    Leave the electrons belong to the selected atoms
    selection: a list, containing the atom numbers
    return: new electrons and matrix
    """
    del_list = []
    for i in xrange(len(electrons)):
        if electrons[i].get_atomNum() not in selection:
            del_list.append(i)
    electrons = np.delete(electrons, del_list).tolist()
    matrix = np.delete(matrix, del_list, 0)
    matrix = np.delete(matrix, del_list, 1)
    return [electrons, matrix]

def valenceElectronsFilter(electrons, matrix):
    """
    Leave only valence electrons
    """
    atom_n_l_index = {}                    # Create index
    for i in xrange(len(electrons)):
        atom, n, l = electrons[i].get_atomNum(), electrons[i].get_electron_n(), electrons[i].get_electron_l()
        if atom not in atom_n_l_index.keys():
            atom_n_l_index[atom] = [(n, l)]
        else:
            atom_n_l_index[atom].append((n,l))
    for atom in atom_n_l_index.keys():
        atom_n_l_index[atom] = sorted(atom_n_l_index[atom])

    valenceDictionary = {}                 # A mark of valence electrons
    for atom in atom_n_l_index.keys():
        if (6, 0) in atom_n_l_index[atom]:                  # 6s, 6p, 6d, 5d, 4f
            valenceDictionary[atom] = [(6, 0), (6, 1), (6, 2), (5, 2), (4, 3)]
        elif (5, 0) in atom_n_l_index[atom]:                # 5s, 5p, 5d, 4d
            valenceDictionary[atom] = [(5, 0), (5, 1), (5, 2), (4, 2)]
        elif (4, 0) in atom_n_l_index[atom]:                # 4s, 4p, 4d, 3d
            valenceDictionary[atom] = [(4, 0), (4, 1), (4, 2), (3, 2)]
        elif (3, 0) in atom_n_l_index[atom]:                # 3s, 3p, 3d
            valenceDictionary[atom] = [(3, 0), (3, 1), (3, 2)]
        elif (2, 0) in atom_n_l_index[atom]:                # 2s, 2p
            valenceDictionary[atom] = [(2, 0), (2, 1)]
        elif (1, 0) in atom_n_l_index[atom]:                # 1s
            valenceDictionary[atom] = [(1, 0)]        

    delete_list = []                       # A mark for deleting electrons
    for i in xrange(len(electrons)):
        atom, n, l = electrons[i].get_atomNum(), electrons[i].get_electron_n(), electrons[i].get_electron_l()
        if (n, l) not in valenceDictionary[atom]:
            delete_list.append(i)

    electrons = np.delete(electrons, delete_list).tolist()
    matrix = np.delete(matrix, delete_list, 0)
    matrix = np.delete(matrix, delete_list, 1)

    return electrons, matrix
    
def merge_electrons(electrons, index, deleteFlag = True):
    """
    Merge electrons to index's first number
    index: a list from i to j
    return: new list
    """
    index.sort()
    for i in index[1:]:
        new_occupation = electrons[index[0]].get_electron_occupation() + electrons[i].get_electron_occupation()
        electrons[index[0]].change_electron_occupation(new_occupation)
   
    if deleteFlag == True:
        electrons = np.delete(electrons, index[1:]).tolist()

        for electron in electrons:                             # Clear the ms attribute
            electron.change_electron_ms(None)

    return electrons

def merge_sym_matrix(matrix, index, deleteFlag = True):
    """
    Merge the rows and columns of a symmetric matrix to index's first row and column
    index: a list from i to j
    return: new matrix
    """
    index.sort()
    for i in index[1:]:
        matrix[index[0],:] += matrix[i,:]
        matrix[:,index[0]] += matrix[i,:]
    if deleteFlag == True:
        for i in index[:0:-1]:
            matrix = np.delete(matrix,i,0)
            matrix = np.delete(matrix,i,1)
    return matrix

def merge_ms(electrons, matrix):
    """
    Merge electrons which only differ in ms
    """
    atom_n_l_index = {}                  # Create index
    for i in xrange(len(electrons)):
        (atom, n, l) = electrons[i].get_atomNum(), electrons[i].get_electron_n(), electrons[i].get_electron_l()
        if (atom, n, l) not in atom_n_l_index.keys():
            atom_n_l_index[(atom, n, l)] = [i]
        else:
            atom_n_l_index[(atom, n, l)].append(i)

    for key in atom_n_l_index.keys():    # Merge
        electrons = merge_electrons(electrons, atom_n_l_index[key], deleteFlag = False)
        matrix = merge_sym_matrix(matrix, atom_n_l_index[key], deleteFlag = False)

    # Delete the extra information
    delete_list = []
    for value in atom_n_l_index.values():
        delete_list += value[1:]
    electrons = np.delete(electrons, delete_list).tolist()
    matrix = np.delete(matrix, delete_list, 0)
    matrix = np.delete(matrix, delete_list, 1)
    for electron in electrons:                             # Clear the ms attribute
        electron.change_electron_ms(None)

    return electrons, matrix

def merge_atom(electrons, matrix):
    """
    Merge electrons which only differ in ms
    """
    atom_index = {}                    # Create index
    for i in xrange(len(electrons)):
        atom = electrons[i].get_atomNum()
        if atom not in atom_index.keys():
            atom_index[atom] = [i]
        else:
            atom_index[atom].append(i)
    
    for atom in atom_index.keys():    # Merge
        electrons = merge_electrons(electrons, atom_index[atom], deleteFlag = False)
        matrix = merge_sym_matrix(matrix, atom_index[atom], deleteFlag = False)

    # Delete the extra information
    delete_list = []
    for value in atom_index.values():
        delete_list += value[1:]
    electrons = np.delete(electrons, delete_list).tolist()
    matrix = np.delete(matrix, delete_list, 0)
    matrix = np.delete(matrix, delete_list, 1)
    for electron in electrons:                             # Clear the n, l, ms attribute
        electron.change_electron_ms(None)
        electron.change_electron_l(None)
        electron.change_electron_n(None)

    return electrons, matrix



############ Main Code ###########

clock0 = time.time()           # TEST

# Step1: generate an array
MullikenFlag = False
num = 0   # TEST
listLikeArray = []
listElectrons = []
for line in open(filename):
    num += 1   # TEST
    #print(str(num))   # TEST
    if MullikenFlag == False:           # Start flag
        if len(line) == 52 and ("Population analysis for representation" in line):
            MullikenFlag = True
            continue
    else:
        try:
            assert '.' in line
            if line[3] != ' ':                                   # Avoid line break
                atomNum = int(line[0:4])
                electron_n = int(line[5])
                electron_l = int(line[6])
                electron_ms = int(line[7:9])
                occupation = float(line[9:16])

                info = line[16:].strip()
                splitting = info.split('.')
            else:
                splitting += line.strip().split('.')
                del listLikeArray[-1]                             # Delete the incomplete array unit
                del listElectrons[-1]

            while '' in splitting:
                splitting.remove('')

            splitting = map(lambda x: int(x), splitting)
            listLikeArray.append(splitting)
            listElectrons.append(Electron(atomNum, electron_n, electron_l, electron_ms, occupation))

        except AssertionError:                                     # End of the Mulliken array
            break

assert len(listLikeArray) == len(listLikeArray[-1])
MullikenPopulationArray = np.zeros([len(listLikeArray), len(listLikeArray[-1])])
for i in xrange(len(listLikeArray)):
    for j in xrange(len(listLikeArray[i])):
        MullikenPopulationArray[i,j] = listLikeArray[i][j]
MullikenPopulationArray /= 1000.0

mirror = MullikenPopulationArray.transpose()                      # Generate the other half of the triangle
MullikenPopulationArray += mirror
for i in xrange(len(listLikeArray)):
    MullikenPopulationArray[i,i] /= 2

del listLikeArray

clock1 = time.time()    # TEST
print(clock1-clock0)

# Step 1.x: Simplification
if atom_mergeSelection == True:
    valanceElectronsSelection, ms_mergeSelection = False, False
if atomSelection != None:                                        # Atom filter
    listElectrons, MullikenPopulationArray = atomFilter(listElectrons, MullikenPopulationArray, atomSelection)

if valanceElectronsSelection == True:                            # Valence electrons filter
    listElectrons, MullikenPopulationArray = valenceElectronsFilter(listElectrons, MullikenPopulationArray)

if ms_mergeSelection == True:                                    # ms merging
    listElectrons, MullikenPopulationArray = merge_ms(listElectrons, MullikenPopulationArray)

if atom_mergeSelection == True:                                  # Atom merging
    listElectrons, MullikenPopulationArray = merge_atom(listElectrons, MullikenPopulationArray)

quick_sort(listElectrons, MullikenPopulationArray)               # Sort the result
reIndex(listElectrons)
clock2 = time.time()        # TEST
print(clock2-clock1)

# Step 2: plotting
plt.imshow(MullikenPopulationArray, interpolation = "nearest", cmap = color)
plt.colorbar()
plt.show()

# Step 3: print electrons list
outfile = open(position + "ElectronsList.txt", "w")
num = 0
for electron in listElectrons:
    outfile.write("%7i  --  " % num + str(electron) + "\n")
    num += 1
outfile.close()
