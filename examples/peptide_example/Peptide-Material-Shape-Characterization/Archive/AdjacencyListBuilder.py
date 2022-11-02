'''
@author: Kassandra Schmidt 

Interaction Count
v.2.01.2019
'''

import math

def main(data,dimensions):

    #Open a Dialog box to allow the user to choose the file they wish to use
    #root = Tk()
    #root.filepath = filedialog.askopenfilename()

    #If the file is not a .gro file, do not allow the program to run
    #if not root.filepath.endswith('.gro'):
    #    print('This is not a .gro file. Please try again with a .gro file.')
    #    return

    #Load the data in the file into a variable and remove the new line characters from the data
    #gro_file = open(root.filepath, 'r')
    #data_list = gro_file.readlines()
    #no_new_line_list = [line.rstrip('\n') for line in data_list]

    #nnll_length = len(no_new_line_list)

    #Store the simulation data not characterizing the particles in the system
    #box_dimensions = no_new_line_list[nnll_length - 1]
    #system_name = no_new_line_list[0]
    #number_of_particles = int(no_new_line_list[1])

    #Check to make sure that the system is big enough to actually run the code on
    #if number_of_particles < 2:
    #    print('There are fewer than two particles in the system. Therefore, there can be no neighbors.')
    #    return

    #Remove the simulation data not characterizing the particles in the system
    #no_new_line_list.remove(no_new_line_list[nnll_length - 1])
    #no_new_line_list.remove(no_new_line_list[0])
    #no_new_line_list.remove(no_new_line_list[0])

    #Transform the box dimensions from strings into floats
    #dimensions = getDimensions(box_dimensions)

    #Extract the particle data from the .gro file
    #residue_name_and_number_list, bead_placement_list, atom_id_list, x_position_list, y_position_list, z_position_list = getParticleData(
    #    no_new_line_list)

    #Check to make sure the box is big enough to have more than one bin
    x_distance2 = dimensions[0] * dimensions[0]
    y_distance2 = dimensions[1] * dimensions[1]
    z_distance2 = dimensions[2] * dimensions[2]

    #Define the lists
    x_position_list = data[:]['x']
    y_position_list = data[:]['y']
    z_position_list = data[:]['z']

    box_diagonal_distance = math.sqrt(x_distance2 + y_distance2 + z_distance2)

    #The cut-off distance
    r_max = 2*1.2  # nm
    neighbors = []

    number_of_peptides = len(x_position_list)

    #If the box is not big enough, every particle is a neighbor to every other particle
    if box_diagonal_distance < r_max:
        neighbors = atom_id_list
    #If the box is big enough, determine which particles are interacting with which other particles
    else:
        neighbors = binning(
            number_of_peptides, x_position_list, y_position_list, z_position_list, dimensions, r_max)

    #User can then create a variety of print statements to analyze their data

    #for particle in neighbors:
       #print(particle)
    return neighbors


def getDimensions(dimension_list_string):
    dimension_list = []

    split_box_dimensions = dimension_list_string.split()

    for x in split_box_dimensions:
       to_float = float(x)
       dimension_list.append(to_float)

    return dimension_list


def getParticleData(data_list_string):
    #Initialize empty lists for data storage
    residue_name_and_number_list = []
    bead_placement_list = []
    atom_id_list = []
    x_position_list = []
    y_position_list = []
    z_position_list = []

    data_list_string_length = len(data_list_string)

    #If the particle is a peptide, add the detail and position information to the appropriate lists
    for x in range(data_list_string_length):
        split_data_list = data_list_string[x].split()
        if 'PW' not in split_data_list[0] and 'ION' not in split_data_list[0]:
            residue_name_and_number_list.append(split_data_list[0])
            bead_placement_list.append(split_data_list[1])
            atom_id_list.append(split_data_list[2])
            x_position_list.append(float(split_data_list[3]))
            y_position_list.append(float(split_data_list[4]))
            z_position_list.append(float(split_data_list[5]))

    return residue_name_and_number_list, bead_placement_list, atom_id_list, x_position_list, y_position_list, z_position_list


def binning(number_of_beads, x_position_list, y_position_list, z_position_list, box_dimensions, r_max):
    #Initialize lists
    bin_list = []
    neighbors = []

    position_list = [x_position_list, y_position_list, z_position_list]

    #Calculate the number of bins. Must be a whole number. Fractions of bins makes no sense
    number_of_x_bins = int(box_dimensions[0] / r_max)
    number_of_y_bins = int(box_dimensions[1] / r_max)
    number_of_z_bins = int(box_dimensions[2] / r_max)

    number_of_bins_list = [number_of_x_bins, number_of_y_bins, number_of_z_bins]

    #Since the number of bins might not have been a whole number, adjust the dimensions of the bins accordingly
    x_bin_width = box_dimensions[0] / number_of_x_bins
    y_bin_width = box_dimensions[1] / number_of_y_bins
    z_bin_width = box_dimensions[2] / number_of_z_bins

    bin_width_list = [x_bin_width, y_bin_width, z_bin_width]

    #Initialize list for keeping track of what beads are in what bins
    for x in range(0, number_of_x_bins):
        bin_list.append([])
        for y in range(0, number_of_y_bins):
            bin_list[x].append([])
            for z in range(0, number_of_z_bins):
                bin_list[x][y].append([])

    #Put ith particle in a bin based upon its coordinates. Add a list to the neighbor list for each bead
    for i in range(number_of_beads):
        
        bins_for_particle_placement = determineBin(i, position_list, number_of_bins_list, bin_width_list)

        bin_list[bins_for_particle_placement[0]][bins_for_particle_placement[1]][bins_for_particle_placement[2]].append(i)

        neighbors.append([])
    
    for j in range(number_of_beads):
        #Determine what bin the particle we are looking at is in.         
        bins_for_particle_location = determineBin(j, position_list, number_of_bins_list, bin_width_list)

        #Check that bin plus the 26 bins around it for neighbors within the cut-off distance
        for adjusted_x, adjusted_y, adjusted_z in determineIndices(bins_for_particle_location[0], bins_for_particle_location[1], bins_for_particle_location[2]):
                   determineNeighbors(j, number_of_x_bins, number_of_y_bins, number_of_z_bins, adjusted_x, adjusted_y,
                                     adjusted_z, x_position_list, y_position_list, z_position_list, bin_list, box_dimensions, neighbors, r_max)
    
    #Return the list of all of the neighbors
    return neighbors

def determineBin(index, position_list, number_of_bins_list, bin_width_list):
    particle_x_bin_number = int(position_list[0][index] / bin_width_list[0])

    if particle_x_bin_number >= number_of_bins_list[0]:
        particle_x_bin_number = number_of_bins_list[0] - 1

    particle_y_bin_number = int(position_list[1][index] / bin_width_list[1])

    if particle_y_bin_number >= number_of_bins_list[1]:
        particle_y_bin_number = number_of_bins_list[1] - 1 

    particle_z_bin_number = int(position_list[2][index] / bin_width_list[2])

    if particle_z_bin_number >= number_of_bins_list[2]:
        particle_z_bin_number = number_of_bins_list[2] - 1
    
    particle_bin_number_list = [particle_x_bin_number, particle_y_bin_number, particle_z_bin_number]

    return particle_bin_number_list

def determineNeighbors(j, number_of_x_bins, number_of_y_bins, number_of_z_bins, particle_x_bin, particle_y_bin, particle_z_bin, x_position_list, y_position_list, z_position_list, bin_list, box_dimensions, neighbors, r_max):

    #Adjust which bin we are looking in based upon whether periodic boundary conditions need to be employed
    if particle_x_bin >= number_of_x_bins:
        particle_x_bin = 0
    elif particle_x_bin < 0:
        particle_x_bin = number_of_x_bins - 1

    if particle_y_bin >= number_of_y_bins:
        particle_y_bin = 0
    elif particle_y_bin < 0:
        particle_y_bin = number_of_y_bins - 1

    if particle_z_bin >= number_of_z_bins:
        particle_z_bin = 0
    elif particle_z_bin < 0:
        particle_z_bin = number_of_z_bins - 1

    #Pick a particle in the bin we are currently in
    for k in bin_list[particle_x_bin][particle_y_bin][particle_z_bin]:
        #Make sure that the particle we pick is not the particle we are trying to find neighbors for and that the new particle is not already in the neighbor list
        if j != k and k not in neighbors[j]:
            #Get the positions of both of the particles
            j_x_position = x_position_list[j]
            j_y_position = y_position_list[j]
            j_z_position = z_position_list[j]

            k_x_position = x_position_list[k]
            k_y_position = y_position_list[k]
            k_z_position = z_position_list[k]

            #Calculate the distance between the particles via the distance formula
            distance_between_j_k = distanceCalculation(
                j_x_position, j_y_position, j_z_position, k_x_position, k_y_position, k_z_position, box_dimensions)
            
            #If the distance between the two particles is less than the cut-off distance, add the particles to each other's neighbor lists
            if distance_between_j_k < r_max:
                neighbors[j].append(k)
                neighbors[k].append(j)


def distanceCalculation(j_x_position, j_y_position, j_z_position, k_x_position, k_y_position, k_z_position, box_dimensions):

    #Taking the absolute value of the distance prevents needing to checking for values that are less than or equal to -0.5*box_dimenions[q]
    x_value = abs(j_x_position - k_x_position)
    y_value = abs(j_y_position - k_y_position)
    z_value = abs(j_z_position - k_z_position)

    #Check for Periodic Boundary Conditions and Adjust Values Accordingly
    if x_value > (box_dimensions[0] * 0.5):
        x_value = x_value - box_dimensions[0]
    if y_value > (box_dimensions[1] * 0.5):
        y_value = y_value - box_dimensions[1]
    if z_value > (box_dimensions[2] * 0.5):
        z_value = z_value - box_dimensions[2]

    #Standard Distance Formula
    x2 = x_value * x_value
    y2 = y_value * y_value
    z2 = z_value * z_value

    distance = math.sqrt(x2 + y2 + z2)

    return distance


def determineIndices(x, y, z):
    #Calculate all combinations of boxes around the box the particle we are looking at is in
    for i in range(-1, 2):
        for j in range(-1, 2):
            for k in range(-1, 2):
                yield x + i, y + j, k + z
