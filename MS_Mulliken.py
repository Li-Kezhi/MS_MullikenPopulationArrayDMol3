import numpy as np
import matplotlib.pyplot as plt
import time      # TEST
import sys   

position = "F:\\BaiduYunDownload\\test\\"
filename = position + "SO2.outmol"

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
        return "Atom %i: n = %i, l = %i, ms = %i" %(self.atomNum, self.electron_n, self.electron_l, self.electron_ms)
    def __str__(self):
        return "Atom %i: n = %i, l = %i, ms = %i" %(self.atomNum, self.electron_n, self.electron_l, self.electron_ms)
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
MullikenPopulationArray /= 2000.0

mirror = MullikenPopulationArray.transpose()                      # Generate the other half of the triangle
MullikenPopulationArray += mirror
for i in xrange(len(listLikeArray)):
    MullikenPopulationArray[i,i] /= 2

clock1 = time.time()    # TEST
print(clock1-clock0)

print(listElectrons)        # TEST
quick_sort(listElectrons, MullikenPopulationArray)                                # Sort the result
reIndex(listElectrons)
clock2 = time.time()        # TEST
print(clock2-clock1)
print(listElectrons)        # TEST

# Step 2: plotting
plt.imshow(MullikenPopulationArray, interpolation = "nearest", cmap = "bone")
plt.colorbar()
plt.show()

# Step 3: print electrons list
outfile = open(position + "ElectronsList.txt", "w")
num = 0
for electron in listElectrons:
    outfile.write("%7i  --  " % num + str(electron) + "\n")
    num += 1
outfile.close()
