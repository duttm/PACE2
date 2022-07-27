import sys
import random
import math
import numpy as np
import subprocess

# References: https://www.geeksforgeeks.org/k-means-clustering-introduction/

# Purpose of this code: To determine if simulation should be extended or not

# How is this determined ?: 

# 1. We take the input structure as a gro file (coordinate file). This file could contain: 

#       a. A closed vesicle (desirable structure)
#       b. A flat membrane (undesirable structure)
#       c. Something between a and b (undesirable structure)

# 2. This code extracts the lipid coordinates from the gro file. These are the main coordinates that will determine the status of the simulation. 

# 3. This code will evaluate the characteristic length of the structure using K Means clustering. 

# 4. A predetermined cut-off is used to do the make the final decision (whether the simulation should be extended or not). The logic

# is :    If (characteristic length of the structure) > is cut-off: 
#                   then -> the simulation is terminated.
#         Else:
#                     the simulation is extended. 


# execution : python < name of input coordinate file > <cut-off>

# < name of input coordinate file >: Please enter the gromacs coordinates file (should have coordinates of all molecules, i.e. dendrons, lipids, water and ions)

# <cut-off>: Value should be in nanometers. For this class of systems, its between 10 to 11 nanometers. 


def ReadData(filename,a,NOF):
    # Read the gromacs input file

    y = 0
    f = open(filename, "r")
    lines = f.readlines()

    result = []
    for x in lines:
        result.append(x)
        y = y + 1

    f.close()


    counter = 0
    refine = []
    x = 0

    while x <= a * (NOF + 1):
        if x == 0 or x == 1:  # skip first 2 lines
            counter = counter + 1
            x = x + 1
           

        elif 2 <= counter <= a + 1:
            refine.append(result[x])
            x = x + 1
            counter = counter + 1
           

        elif counter == a + 2:
            counter = 2
            x = x + 3
            if x > y:
                break
       

        else:
            print(0)

    file_out = open("array.txt", "w")

    # admit values into a large 2D array

    rows, cols = (a * NOF, 4)
    array = [[0 for i in range(cols)] for j in range(rows)]

    for i in range(a * NOF):
        t = refine[i][0:8]
        id = refine[i][15:20]  # Column values may vary on the basis of gromacs version
        x = refine[i][21:28]
        y = refine[i][29:36]
        z = refine[i][37:44]
        array[i][0] = id
        array[i][1] = x
        array[i][2] = y
        array[i][3] = z

    for i in range(rows):
        file_out.writelines(str(array[i]) + '\n')
        # view array.txt to check
    file_out.close()

    # 2d big to 3D small conversion

    rows, cols, pages = (a, 4, NOF)
    array3d = [[[0 for k in range(pages)] for i in range(cols)] for j in range(rows)]

    h = 0

    for k in range(pages):
        for i in range(rows):
            array3d[i][0][k] = float(array[h][0])
            array3d[i][1][k] = float(array[h][1])
            array3d[i][2][k] = float(array[h][2])
            array3d[i][3][k] = float(array[h][3])
            h = h + 1

    file_out1 = open("array3d.txt", "w")

    # view array3d.txt to check

    for k in range(pages):
        for i in range(rows):
            file_out1.writelines(
                str(array3d[i][1][k]) + ',' + str(array3d[i][2][k]) + ',' + str(array3d[i][3][k]) + '\n')

    file_out1.close()

    fileName = "array3d.txt"
    f = open(fileName, 'r')
    lines = f.read().splitlines()
    f.close()

    items = []

    for i in range(1, len(lines)):
        line = lines[i].split(',')
        itemFeatures = []

        for j in range(len(line)):
            v = float(line[j])  # Convert feature value to float
            itemFeatures.append(v)  # Add feature value to dict

        items.append(itemFeatures)

    random.shuffle(items)

    ##fig = plt.figure()
    ##ax = fig.add_subplot(111, projection='3d')

    return items;

# References: https://www.geeksforgeeks.org/k-means-clustering-introduction/

# A kmeans function
def FindColMinMax(items):
    n = len(items[0]);

    minima = [sys.maxsize for i in range(n)];
    maxima = [-sys.maxsize - 1 for i in range(n)];

    for item in items:
        for f in range(len(item)):

            if (item[f] < minima[f]):
                minima[f] = item[f];

            if (item[f] > maxima[f]):
                maxima[f] = item[f];

    # print(maxima, minima)
    return minima, maxima;
# References: https://www.geeksforgeeks.org/k-means-clustering-introduction/

# A kmeans function

def InitializeMeans(items, k, cMin, cMax):
    # Initialize means to random numbers between
    # the min and max of each column/feature
    f = len(items[0]);  # number of features
    means = [[0 for i in range(f)] for j in range(k)];

    for mean in means:
        for i in range(len(mean)):
            # Set value to a random float
            # (adding +-1 to avoid a wide placement of a mean)
            mean[i] = random.uniform(cMin[i] + 1, cMax[i] - 1);

    # print(means)
    return means;
# References: https://www.geeksforgeeks.org/k-means-clustering-introduction/

# A kmeans function

def EuclideanDistance(x, y):
    S = 0  # The sum of the squared differences of the elements
    for i in range(len(x)):
        S += math.pow(x[i] - y[i], 2)

    return math.sqrt(S)  # The square root of the sum

# References: https://www.geeksforgeeks.org/k-means-clustering-introduction/

# A kmeans function
def UpdateMean(n, mean, item):
    for i in range(len(mean)):
        m = mean[i];
        m = (m * (n - 1) + item[i]) / float(n);
        mean[i] = round(m, 3);

    return mean;

# References: https://www.geeksforgeeks.org/k-means-clustering-introduction/

# A kmeans function
def Classify(means, item):
    # Classify item to the mean with minimum distance
    minimum = sys.maxsize;
    index = -1;

    for i in range(len(means)):

        # Find distance from item to mean
        dis = EuclideanDistance(item, means[i]);

        if (dis < minimum):
            minimum = dis;
            index = i;

    return index;

# References: https://www.geeksforgeeks.org/k-means-clustering-introduction/

# A kmeans function
def CalculateMeans(k, items, maxIterations=100000):
    # Find the minima and maxima for columns
    cMin, cMax = FindColMinMax(items);

    # Initialize means at random points
    means = InitializeMeans(items, k, cMin, cMax);

    # print('means \n')
    # print(means)

    # Initialize clusters, the array to hold
    # the number of items in a class
    clusterSizes = [0 for i in range(len(means))];

    # An array to hold the cluster an item is in
    belongsTo = [0 for i in range(len(items))];

    # Calculate means
    for e in range(maxIterations):

        # If no change of cluster occurs, halt
        noChange = True;
        for i in range(len(items)):

            item = items[i];

            # Classify item into a cluster and update the
            # corresponding means.
            index = Classify(means, item);

            clusterSizes[index] += 1;
            cSize = clusterSizes[index];
            means[index] = UpdateMean(cSize, means[index], item);

            # Item changed cluster
            if (index != belongsTo[i]):
                noChange = False;

            belongsTo[i] = index;

        # Nothing changed, return
        if (noChange):
            break;

    return means;

# References: https://www.geeksforgeeks.org/k-means-clustering-introduction/

# A kmeans function
def FindClusters(means, items):
    clusters = [[] for i in range(len(means))];  # Init clusters

    for item in items:
        # Classify item into a cluster
        index = Classify(means, item);

        # Add item to cluster
        clusters[index].append(item);

    return clusters;


def main():


    ###### command line section ####################

    filename = sys.argv[1]
    Cutoff = int(sys.argv[2])
    print('filename is %s' % (filename))
    NOF = 1
    print('NOF is %s' % (NOF))    
    #################################################
    
    # Extract lipid coordinates
    
    target = "DPPC.gro" # the name of the lipid is DPPC
    
    cmd = 'echo 8| gmx trjconv -f ' + str(filename) + ' -o ' +  target  + ' -s ' + str(filename) 
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    print (proc.communicate()[0])
    #####################################################################
    
    
    fp = open(target)
    line_count = 0
    for i, line in enumerate(fp):
        if i == 1:
            a = int(line.strip())
    fp.close()
    print('No. of elements in a single frame is %s' % (a))
    #####################################################################
    

    ######################################################################
    
    # K-means needs to repeated x times
    
    # Why ? These lipid based systems should ideally form one aggregate (one larger cluster of lipids). This way, K-means will divide 
    # the large aggregate into 2 equal halves, and hence we would get the most accurate estimation of the characteristic length. 
    # However, in reality, there is one large cluster that has ~95 % of the lipids, and the remaining are distributed randomly.
    # Due to a small concentration of randomly distributed lipids (noise for the algorithm), the k-means may not give incorrect
    # estimations of the characteristic length
    # However, if k-means is repeated sufficient number of times, we get the correct estimation of the characteristic length
    # We have repeated K-means 5 times for these systems.  
    # Number of repetitions is set up the "experiments" variable
    
    experiments=5
    store_dia = [0 for r in range(experiments)]
    for count in range(experiments):


        np.random.seed(19680801)
        print('Doing K Means of 2 in experiment %s' % (count+1))
        items = ReadData(filename,a,NOF)
        # print(items)

        means = CalculateMeans(2, items)

        clusters = FindClusters(means, items)

        centroid_x = [0 for r in range(2)]
        centroid_y = [0 for r in range(2)]
        centroid_z = [0 for r in range(2)]

        for i in range(len(means)):
            centroid_x[i] = float(means[i][0])
            centroid_y[i] = float(means[i][1])
            centroid_z[i] = float(means[i][2])

        ##################################################
        np.random.seed(19680801)
        print('Doing K Means of 1 to find COM in experiment %s '% (count+1))
        items = ReadData(filename,a,NOF)
        # print(items)

        means = CalculateMeans(1, items)

        clusters = FindClusters(means, items)

        com_x = [0 for r in range(2)]
        com_y = [0 for r in range(2)]
        com_z = [0 for r in range(2)]

        for i in range(len(means)):
            com_x = float(means[i][0])
            com_y = float(means[i][1])
            com_z = float(means[i][2])

        diameter_1 = math.sqrt(
            math.pow(centroid_x[0] - com_x, 2) + math.pow(centroid_y[0] - com_y, 2) + math.pow(centroid_z[0] - com_z,
                                                                                               2))
        diameter_2 = math.sqrt(
            math.pow(centroid_x[1] - com_x, 2) + math.pow(centroid_y[1] - com_y, 2) + math.pow(centroid_z[1] - com_z,
                                                                                               2))
        store_dia[count] = (diameter_1 + diameter_2)
        print("Result : % s \n" % (store_dia[count]))

    
    check = 0
    for i in range(experiments):
    
        if (store_dia[i] < Cutoff):
            check=1
              
    f = open("out.txt", "w")
    f.write(str(check))
    f.close()


if __name__ == "__main__":
    main()
