'''
@author: Kassandra Schmidt
@author: Ethan Zang

Interaction Count
v.2.14.2019
'''

# KLS 2062019 - Ethan, going to have to change imports back to 2.7. VS Code is supposed to support 2.7 but doesn't work unless I change the inputs
# KLS 2062019 - changed particle type from 'water' to 'solvent' for genericness
# KLS 2062019 - added box dimensions to the end of both .dat files
# KLS 2062019 - Code works up to 9999 particles. 10,000 causes an alignment shift which messes up the generic algorithm (Discovered by Ethan on DATE).
#   Need to work on implementing a fix for this for next week; will probably have a check for 9999 particles and then a format shift when sorting the data
#   into their respective arrays in createParticleList()
# KLS 2062019 - need to add ION type particle to solvent list; this means changing peptide identification method from not water or ion to something else. Need
#   to make a decision on what this will be. Could just use the else but prefer the current set up.

# KLS 2072019 - Added ION to to if clause and put peptide identification as the outlier (removed elif statement)
# KLS 2072019 - Changed residue_name_and_number into two seperate variables: residue_name_and residue_number
# KLS 2072019 - Switched from splitting the data based upon white space (which was causing the error on 2062019 for > 9,999 particles) to just slicing the data based upon the
#   .gro file column sizes.
# KLS 2072019 - changed the formatting in the 'row' variable in printToText() to reflect how the data was being divided via slicing.


#from Tkinter import *
#from tkFileDialog import *
import math, sys, os


def main(file):
    #Open a Dialog box to allow the user to choose the file they wish to use
#    root = Tk()
    #root.filepath = filedialog.askopenfilename()
#    root.filepath = askopenfilename()
#    file = sys.argv[1]
    # If the file is not a .gro file, do not allow the program to run
#    if not root.filepath.endswith('.gro'):
#        print('This is not a .gro file. Please try again with a .gro file.')
#        return

    try:
        os.remove('solventFile.dat')
    except OSError:
        pass

    try:
        os.remove('peptideFile.dat')
    except OSError:
        pass

    # Load the data in the file into a variable and remove the new line characters from the data
    gro_file = open(file, 'r')
    data_list = gro_file.readlines()
    no_new_line_list = [line.rstrip('\n') for line in data_list]

    nnll_length = len(no_new_line_list)

    # Save the box dimensions, then remove the simulation data not characterizing the particles in the system
    box_dimensions = no_new_line_list[nnll_length - 1]
    # Remove box dimensions
    no_new_line_list.remove(no_new_line_list[nnll_length - 1])
    # Remove System Name
    no_new_line_list.remove(no_new_line_list[0])
    # Remove Total Number of Particles
    no_new_line_list.remove(no_new_line_list[0])

    # Extract the particle data from the .gro file
    createParticleList(no_new_line_list, box_dimensions)


def createParticleList(no_new_line_list, box_dimensions):
    # Initialiaze a list to store the data associated with each particle
    particle_list = []

    for line_of_data in no_new_line_list:
        # Use slicing (use the indices of the characters in the string) on the string to isolate the residue name.
        residue_name = line_of_data[5:10]
        # W - water; PW - Polarized Water; WF - Anti-Freeze Particle (Bigger Radius Water); ION - Counter ions for equilibrating system charge
        if 'W' in residue_name or 'PW' in residue_name or 'WF' in residue_name or 'ION' in residue_name:
            particle_type = 'solvent'
        else:
            particle_type = 'peptide'
        # Use slicing to get various columns of data. The size of the columns was calculated after studying .gro files made by different members of the group.
        particle = Particle(line_of_data[0:5], residue_name, line_of_data[10:15], line_of_data[15:20],
                            line_of_data[20:28], line_of_data[28:36], line_of_data[36:44], particle_type)

        particle_list.append(particle)

    printToText(particle_list, box_dimensions)


def printToText(list_of_particles, box_dimensions):
    # Loop through the lists and format the data in the list based upon where in the row it falls in a .gro file
    for particle in list_of_particles:
        # Format the data based upon the .gro file standard formatting
        row = '{:>5}{:5}{:>5}{:>5}{:8.3f}{:8.3f}{:8.3f}'.format(particle.residue_number, particle.residue_name, particle.bead_placement, particle.atom_id, float(
            particle.x_position), float(particle.y_position), float(particle.z_position))

        # Writes the formatted peptide information text to a .dat file
        if particle.type == 'peptide':
            with open('peptideFile.dat', 'a') as peptide_file:
                peptide_file.write(row+"\n")

        # Writes the formatted solvent information text to a .dat file
        elif particle.type == 'solvent':
            with open('solventFile.dat', 'a') as solvent_file:
                solvent_file.write(row+"\n")

    # Add in box dimensions to the end of the peptide file and the solvent file
    #with open('peptideTextFile.dat', 'a') as peptide_file:
    #    peptide_file.write(box_dimensions)

    with open('solventFile.dat', 'a') as solvent_file:
        solvent_file.write(box_dimensions)


class Particle(object):

    def __init__(self, residue_number, residue_name, bead_placement, atom_id, x_position, y_position, z_position, type_of_bead):

        self.residue_number = residue_number
        self.residue_name = residue_name
        self.bead_placement = bead_placement
        self.atom_id = atom_id
        self.x_position = x_position
        self.y_position = y_position
        self.z_position = z_position
        self.type = type_of_bead