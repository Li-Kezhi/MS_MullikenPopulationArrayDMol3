import numpy as np
import matplotlib.pyplot as plt

filename = "E:\\temp\\For python\\Mulliken\\SO2.outmol"

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


# Step1: generate an array
MullikenFlag = False
num = 0   # TEST
listLikeArray = []
listElectrons = []
for line in open(filename):
    num += 1   # TEST
    #print(str(num))   # TEST
    if MullikenFlag == False:           # Start flag
        if "Population analysis for representation" in line:
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

# Step 2: plotting
plt.imshow(MullikenPopulationArray, interpolation = "nearest", cmap = "bone")
plt.colorbar()
plt.show()
