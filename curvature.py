# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 2017

@author: snfrbert 
based on work of schustmc (Sept 01 2016)
"""
import sys
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, colors

class Bilayer:

    def __init__(self):
        self.positions = []
        self.angstrom = True
        # r and d are given in nanometers
        # if a is True, they have to be converted to 

    def readIn(self, filename, angstrom = True):
        self.angstrom = angstrom
        try:
            f = open(filename)
        except IOError as info: 
            print(info)
            sys.exit(1)
        else:
            x = 0
            y = 0
            z = 0       
            for line in f:
                if (not line.startswith("ATOM")) and (not "PO4" in line):
                    continue
                ls = line.split()
                if self.angstrom == True:
                    x = float(ls[5]) / 10
                    y = float(ls[6]) / 10
                    z = float(ls[7]) / 10
                else:
                    x = float(ls[5])
                    y = float(ls[6])
                    z = float(ls[7])
                self.positions.append((x,y,z))
            f.close()   
    """ 
    reads a pdb file, identifies all the lines with atom postitions and stores 
    the values in a list of tuples
    if angstrom is True, they are given in angstrom and have to be converted to nanometer 
    """

    @staticmethod
    def euclidean(p1, p2):
        if len(p1) != len(p2):
            print("Vectors not of same size in euclidean!")
            sys.exit(-1)
        s = 0
        for i in range(len(p1)):
            s += (p1[i] - p2[i])**2
        return np.sqrt(s)
    """
    return the euclidean distance for two vectors of
    arbitrary size
    """

    @staticmethod
    def euclideanDistanceMatrix(positions):
        result = np.zeros((len(positions), len(positions)))
        for i in range(len(positions)):
            for j in range(len(positions)):
                result[i,j] = Bilayer.euclidean(positions[i], positions[j])
        return result
    """ 
    calculates all pairwise distances between all atom positions and stores those 
    distances in a matrix. 
    alternativ approach:
    """
    @staticmethod
    def computeCom(liste):
        x = y = z = 0
        for pos in liste:
            x += pos[0]
            y += pos[1]
            z += pos[2] 
        try:
            x /= len(liste)
            y /= len(liste)
            z /= len(liste)
        except ZeroDivisionError:
            print("Warning: List of positions is empty")
        return (x,y,z)
    """ 
    calculates the center of mass from a list of atom positions
    """

    def computeComIndices(self,liste):
        x = y = z = 0
        for i in liste:
            x += bl.positions[i][0]
            y += bl.positions[i][1]
            z += bl.positions[i][2] 
        try:
            x /= len(liste)
            y /= len(liste)
            z /= len(liste)
        except ZeroDivisionError:
            print("Warning: List of positions is empty")
        return (x,y,z)
    """ 
    calculates the center of mass from a list of atom positions
    """
    
    @staticmethod
    def maxDist(liste):
        maxDist = 0
        com = Bilayer.computeCom(liste)
        for coordinates in liste:
            dist_to_check = Bilayer.euclidean(coordinates, com)
            if dist_to_check > maxDist:
                maxDist = dist_to_check
        if maxDist == 0:
            print("Error: could not calculate maximum distance!")
            sys.exit()
        return maxDist
    """ 
    calculates the maximal distance of an atom to the center of mass in a given
    list of positions
    """

    
    def getRelevant(self, distanceFromRim):
        indices = []
        indices2 = []
        c = Bilayer.computeCom(self.positions)
        mdist = Bilayer.maxDist(self.positions)
        if distanceFromRim > mdist:
            print("Error: Distance from rim is larger than calculated maximum distance!")
            print(distanceFromRim, mdist)
            sys.exit()
        for index, coordinate in enumerate(self.positions):
            d = Bilayer.euclidean(c, coordinate)
            if (mdist - distanceFromRim - d) > 0:
                indices.append(index)
            else:
                indices2.append(index)
        return indices,indices2
    """ 
    Determines all atoms, which are in the range maxDist-X from the COM. 
    Returns list of the indices which refer to the corresponding indices in the pdb file
    distanceFromRim (=X) is in nanometer
    """

    
    def bilayerSeperation(self, distance, searchRadius, warning=True, DOUBLECHECK=True, START=0):
        indices = self.getRelevant(distance)[0]
        dist_matrix = Bilayer.euclideanDistanceMatrix(self.positions)
        lst_1 = [indices[START]]
        lst_2 = []
        for i in lst_1:
            for j in indices:
                if dist_matrix[i,j] != 0. and dist_matrix[i,j] < searchRadius:
                    if j not in lst_1:
                        lst_1.append(j) 
        for i in indices:
            if i not in lst_1:
                lst_2.append(i)

        if len(lst_2) == 0:
            DOUBLECHECK = False
            if warning:
                print("Warning: ")
                print("List 2 is empty. DistanceToRim not appropriate or chosen searchradius is too large!")
                sys.exit()

        
        if DOUBLECHECK == True:            
            lst_3 = [lst_2[0]]
            lst_4 = []
            for i in lst_3:
                for j in indices:
                    if dist_matrix[i,j] != 0. and dist_matrix[i,j] < searchRadius:
                        if j not in lst_3:
                            lst_3.append(j) 
            for i in indices:
                if i not in lst_3:
                    lst_4.append(i)
            
            if len(lst_1) != len(lst_4):
                print("Warning: Found two different sets of atoms after double check!\nCheck parameter combinations!")
                sys.exit()
            
        return np.sort(lst_1), np.sort(lst_2)
    """ 
    first, use the EuclideanDistance function to calculate all pairwise 
    distances of the atoms. Then calculate the euclidean distance for a given Radius.
    Initialize the two lists and choose a starting atom.
    Find all neighbors in between a given RADIUS. Do this for all found neighbors.
    All other atoms in the indices list are stored in a different list. The
    output are two list of indices (corresponding to atom postitions in pdb file)
    referring to the two layers in a lipid bilayer.
    Note: A double check routine is possible, if True, starts with an atom of the second
    list and checks, if same atom positions are found on second run.
    """

    def getPositions(self,indices):
        x = []
        y = []
        z = []
        for i in indices:
            xi,yi,zi = self.positions[i]
            x.append(xi)
            y.append(yi)
            z.append(zi)

        return (x,y,z)
    """
    returns a tuple with x,y,z values (a list each) for the given indices
    for plotting e.g
    """

    def fitSphere(self,indices):
        # first separate the x, y and z coordinates from each other and store them in individual lists
        x,y,z = self.getPositions(indices)

        # set up least squares system and solve with numpy
        spX = np.array(x)
        spY = np.array(y)
        spZ = np.array(z)
        A = np.zeros((len(spX),4))
        A[:,0] = spX
        A[:,1] = spY
        A[:,2] = spZ
        A[:,3] = 1
        ## assemble the f matrix
        f = np.zeros((len(spX),1))
        f[:,0] = -(spX**2+spY**2+spZ**2)    
        C, residues, rank, singval = np.linalg.lstsq(A,f)
        xc = -C[0][0] / 2
        yc = -C[1][0] / 2
        zc = -C[2][0] / 2

        radius = np.sqrt(xc**2 + yc**2 + zc**2 - C[3][0])
        center = [xc,yc,zc]
        return center, radius


    @staticmethod
    def sphere(center, radius):
        x = []
        y = []
        z = []
        for i in range(0,360,2):      #Kugel mit Radius R um Mittelpunkt (xc,yc,zc) mit allen Winkelkombinationen
            for j in range(0,360,2):
                x.append(center[0]+radius*np.cos(i)*np.cos(j))
                y.append(center[1]+radius*np.cos(i)*np.sin(j))
                z.append(center[2]+radius*np.sin(i))
        return (x,y,z)

    def error_fit(self,indices,center,radius):
        if len(indices) == 0:
            return
        error = 0
        for i in indices:       
            error += (self.euclidean(self.positions[i], center) - radius)**2
        error = np.sqrt(error/len(indices))   #rmsd
        return error

    def checkVesicle(self):
        up,low = self.bilayerSeperation(0,2,False)
        if (len(low) > 2) and (len(up) > 2):
            com1 = self.computeComIndices(up)
            com2 = self.computeComIndices(low)
            if Bilayer.euclidean(com1,com2) < 1:
                print("The patch formed a vesicle")
                return 1
            else:
                print("Probably a flat patch, not connected at the rim!")
                return 2

        return 3



def plotCut(bl, upp, lowp, cut, interactive):
    com_x,com_y,com_z = Bilayer.computeCom(bl.positions)
    relevant,irrelevant = bl.getRelevant(cut)
    #relp = bl.getPositions(relevant)
    irrp = bl.getPositions(irrelevant)
    fig = plt.figure()
    fig.suptitle('Layer1, Layer2 and CutAway', fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(com_x,com_y,com_z, c="b", label = "COM")
    ax.plot(irrp[0],irrp[1], irrp[2],"xr", label = "Cut")
    ax.plot(upp[0],upp[1], upp[2],"xg", label = "Layer1")
    ax.plot(lowp[0], lowp[1], lowp[2],"xb", label = "Layer2")

    ax.set_xlim([-10,40])
    ax.set_ylim([-10,40])
    ax.set_zlim([-10,40])
    ax.set_aspect("equal")
    plt.tight_layout()

    if interactive:
        plt.show()  
    else:
        plt.savefig("cut-plot.png") 

def plotFit(bl,c,r,p,i,interactive):
    x, y, z = Bilayer.sphere(c, r)
    fig = plt.figure()
    fig.suptitle('Sphere Fit for Layer {}'.format(i), fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(p[0],p[1],p[2],"xg", label = "Layer1")
    ax.plot(x, y, z,color="c", alpha=0.3)
    ax.set_xlim([-10,40])
    ax.set_ylim([-10,40])
    ax.set_zlim([-10,40])
    ax.set_aspect("equal")
    plt.tight_layout()

    if interactive:
        plt.show()  
    else:
        plt.savefig("fit-plot-{}.png".format(i)) 

def printHelp():
    print("Usage:")
    print("python3 curvature.py filename [distanceFromRim searchRadius plots interactive]")
    print()
    print("filename is the path to your .pdb file")
    print("distanceToRim is a float in nm")
    print("     default: 2.0")
    print("searchRadius is a float in nm")
    print("     default: 2.0")
    print("plots is True or False (whether you want plots to check the fitting results")
    print("     default: True")
    print("interactive is True or False (whether you want to be asked to proceed for a presumeably flat patch and have interactive plots)")
    print("     default: True")

def readInput():
    l = len(sys.argv)
    distanceFromRim = 2  #2,2 funktioniert gut fuer alles außer vesikel
    searchRadius = 2
    plots = True
    interactive = True
    if l < 1 or l > 6:
        print("Wrong number of parameters")
        printHelp()
        sys.exit()
    if sys.argv[1] == "-h":
        printHelp()
        sys.exit()
    filename = sys.argv[1]
    try:
        if l == 6:
            interactive = eval(sys.argv[5])
            plots = eval(sys.argv[4])
            searchRadius = float(sys.argv[3])
            distanceFromRim = float(sys.argv[2])
        if l == 5:
            plots = eval(sys.argv[4])
            searchRadius = float(sys.argv[3])
            distanceFromRim = float(sys.argv[2])
        if l == 4:
            searchRadius = float(sys.argv[3])
            distanceFromRim = float(sys.argv[2])
        if l == 3:
            distanceFromRim = float(sys.argv[2])
    except:
        print("Make sure your parameters are in the right order and of the right type!")
        printHelp()
        sys.exit()
    return(filename,distanceFromRim,searchRadius,plots,interactive)




filename, distanceFromRim, searchRadius, plots, interactive = readInput()



bl = Bilayer()       
bl.readIn(filename)
v = bl.checkVesicle()

if v == 1:
    distanceFromRim = 0
if v == 2:
    if interactive:
        ci = input("Try to fit a sphere anyway (please doublecheck result!!) y/[n]")
        if ci.lower() != "y":
            print("Finished")
            sys.exit()

up,low = bl.bilayerSeperation(distanceFromRim, searchRadius) 

c_u, r_u = bl.fitSphere(up)
print("Fitted sphere for first layer with center ({:.2f},{:.2f},{:.2f}) and radius {:.2f}.".format(*(c_u+[r_u])))
print("The rmsd is: {:.4f} nm".format(bl.error_fit(up, c_u, r_u)))
print()

c_l, r_l = bl.fitSphere(low)
print("Fitted sphere for second layer with center ({:.2f},{:.2f},{:.2f}) and radius {:.2f}.".format(*(c_l+[r_l])))
print("The rmsd is: {:.4f} nm".format(bl.error_fit(low, c_l, r_l)))
print()

print("Curvature of layer1: {:.4f} = 1/nm".format(1/r_u))
print("Curvature of layer2: {:.4f} = 1/nm".format(1/r_l))
print("Mean curvature of layer1 and layer2: {:.4f} = 1/nm".format((1/r_l+1/r_u)/2))
print("Result curvature (2/(r1+r2)) is: {:.4f} = 1/nm".format(2/(r_l+r_u)))

print()
print()

if plots:
    print("Plotting for validation")

    upp = bl.getPositions(up)  
    lowp = bl.getPositions(low)  
     
    plotCut(bl, upp, lowp, distanceFromRim, interactive)
    plotFit(bl, c_u,r_u,upp,1,interactive)
    plotFit(bl, c_l,r_l,lowp,2,interactive)



print("Finished")
