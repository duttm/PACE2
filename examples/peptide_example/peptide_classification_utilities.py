import tensorflow as tf
import numpy as np
import glob
import random
import GroConverter
import matplotlib.pyplot as plt
import imageio.v2 as imageio
import os

from mpl_toolkits.mplot3d import Axes3D
from sklearn.model_selection import train_test_split
from tensorflow.keras import datasets, layers, models
from typing import Dict


def generate_images():
    root_dir = os.getcwd()
    os.path.join(root_dir, 'training_data')

    if not os.path.exists('images'):
        os.mkdir('images')
    for root, dirs, files in os.walk(root_dir):
        # Loop through all the subdirectories
        for cur_directory in dirs:
            temp = 'training_data/' + cur_directory + '/*.gro'

            # Get label from directory name
            label = cur_directory

            # Loop through current directory
            for file_name in glob.glob(temp):
                print(file_name)
                create_image(file_name, False)


def create_image(gro_file: str, predicting: bool) -> str:
    # Convert gro file to .dat peptide + solvent files and read data from peptide file
    GroConverter.main(gro_file)
    file = 'peptideFile.dat'
    peptide_file = open(file, 'r')
    peptide_data = peptide_file.readlines()
    data = _readPeptideData(peptide_data)

    # Reassign identifierIDs
    for i in range(len(data)):
        data[i]['identifierID'] = i + 1

    # Grab box dimensions
    box_dimensions = _getBoxDimensions()

    # Temporarily rename data for shifting
    shifting_data = data

    # Perform all of shifting 2 times
    for i in range(2):

        # Grab random peptide particle
        random_peptide_index = random.randint(0, len(shifting_data) - 1)

        x_shift_factor = shifting_data[random_peptide_index]['x'] - box_dimensions[0] / 2
        y_shift_factor = shifting_data[random_peptide_index]['y'] - box_dimensions[1] / 2
        z_shift_factor = shifting_data[random_peptide_index]['z'] - box_dimensions[2] / 2

        # Shift first based on random peptide particle
        data_firstShift = _shift(shifting_data, x_shift_factor, y_shift_factor, z_shift_factor, box_dimensions)

        # Temp reassign data to shift 40 more times
        shifting_data = data_firstShift

        # Find center of mass
        x_mean = np.mean(shifting_data[:]['x'])
        y_mean = np.mean(shifting_data[:]['y'])
        z_mean = np.mean(shifting_data[:]['z'])

        # Shift
        shifting_data = _shift(shifting_data, x_mean - box_dimensions[0] / 2, y_mean - box_dimensions[1] / 2,
                              z_mean - box_dimensions[2] / 2, box_dimensions)

        # Iteratively shift 39 more times based on center of mass
        for j in range(40):
            # Find center of mass
            x_mean = np.mean(shifting_data[:]['x'])
            y_mean = np.mean(shifting_data[:]['y'])
            z_mean = np.mean(shifting_data[:]['z'])

            # Shift
            shifting_data = _shift(shifting_data, x_mean + box_dimensions[0] / 2, y_mean + box_dimensions[1] / 2,
                                  z_mean + box_dimensions[2] / 2, box_dimensions)

    # Plot final post-shifting data
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    ax.scatter3D(shifting_data[:]['x'], shifting_data[:]['y'], shifting_data[:]['z'], cmap='hsv')

    if predicting:
        image_path = os.path.join('images', 'predicting', 'test_image')
        if not os.path.exists('images/predicting'):
            os.makedirs('images/predicting')
        fig.savefig(image_path)
    else:
        image_name_components = gro_file.split('/')
        file_name = image_name_components[2].split('.')[0]
        image_path = os.path.join('images', image_name_components[1] + '_' + file_name)
        fig.savefig(image_path)
    return image_path + '.png'


# Read peptide data (only backbone), returned as structured numpy array
def _readPeptideData(masterData):
    # Find pattern of backbones
    pattern_of_backbones = 0
    for i in range(1, len(masterData)):
        cur_line = masterData[i].split()
        pattern_of_backbones += 1

        if cur_line[1] == "BB":
            break

    # Populate structured numpy array using masterData
    dt = np.dtype(
        [('residueID', np.unicode_, 16), ('moleculeID', np.unicode_, 16), ('identifierID', int), ('x', np.float64),
         ('y', np.float64), ('z', np.float64)])
    data = np.empty(len(masterData) // pattern_of_backbones, dtype=dt)

    # Only record data that is backbone
    index = 0
    for i in range(len(masterData)):
        if (i % pattern_of_backbones == 0):
            cur_line = masterData[i].split()
            data[index]['residueID'] = cur_line[0]
            data[index]['moleculeID'] = cur_line[1]
            data[index]['identifierID'] = cur_line[2]
            data[index]['x'] = round(float(cur_line[3]), 3)
            data[index]['y'] = round(float(cur_line[4]), 3)
            data[index]['z'] = round(float(cur_line[5]), 3)
            index += 1

    return data


# Read box dimensions
def _getBoxDimensions():
    # Read solvent txt data as second argument
    file2 = "solventFile.dat"
    solvent_file = open(file2, 'r')
    solvent_data = solvent_file.readlines()

    # Extract box dimensions from end of file
    temp = solvent_data[len(solvent_data) - 1].split()

    box_dimensions = [round(float(temp[0]), 3), round(float(temp[1]), 3), round(float(temp[2]), 3)]

    return box_dimensions


# Read solvent data, building array with both peptide and solvent data and returning box dimensions
def _readSolventData(peptideData):
    # Read solvent txt data as second argument
    file2 = "solvent_file.dat"
    solvent_file = open(file2, 'r')
    solventData = solvent_file.readlines()

    solventData = np.delete(solventData, len(solventData) - 1)

    dt = np.dtype(
        [('residueID', np.unicode_, 16), ('moleculeID', np.unicode_, 16), ('identifierID', int), ('x', np.float64),
         ('y', np.float64), ('z', np.float64)])
    complete_data = np.empty(len(peptideData) + len(solventData), dtype=dt)

    for i in range(len(peptideData)):
        complete_data[i] = peptideData[i]

    for i in range(len(solventData)):
        cur_line = solventData[i].split()
        complete_data[i + len(peptideData)]['residueID'] = 'solvent'
        complete_data[i + len(peptideData)]['moleculeID'] = 'S'
        complete_data[i + len(peptideData)]['identifierID'] = i + len(peptideData)

        if len(cur_line) == 6:
            complete_data[i + len(peptideData)]['x'] = cur_line[3]
            complete_data[i + len(peptideData)]['y'] = cur_line[4]
            complete_data[i + len(peptideData)]['z'] = cur_line[5]
        else:
            complete_data[i + len(peptideData)]['x'] = cur_line[2]
            complete_data[i + len(peptideData)]['y'] = cur_line[3]
            complete_data[i + len(peptideData)]['z'] = cur_line[4]

    return complete_data


# Screen out outliers
def _dropOutliers(data):
    # Calculate Q1 and Q3 in each dimension
    x_q1 = np.percentile(data[:]['x'], 25)
    x_q3 = np.percentile(data[:]['x'], 75)
    y_q1 = np.percentile(data[:]['y'], 25)
    y_q3 = np.percentile(data[:]['y'], 75)
    z_q1 = np.percentile(data[:]['z'], 25)
    z_q3 = np.percentile(data[:]['z'], 75)

    # Calculate IQR in each dimension
    x_iqr = x_q3 - x_q1
    y_iqr = y_q3 - y_q1
    z_iqr = z_q3 - z_q1

    # Calculate allowed tolerance in each dimension
    x_max = x_q3 + 1.5 * x_iqr
    x_min = x_q1 - 1.5 * x_iqr
    y_max = y_q3 + 1.5 * y_iqr
    y_min = y_q1 - 1.5 * y_iqr
    z_max = z_q3 + 1.5 * z_iqr
    z_min = z_q1 - 1.5 * z_iqr

    # Populate list of indices that follow out of allowed tolerance in any dimension, and delete those indices
    outlier_indices = []
    for i in range(len(data)):
        if ((data[i]['x'] < x_min) | (data[i]['x'] > x_max) | (data[i]['y'] < y_min) | (data[i]['y'] > y_max) | (
                data[i]['z'] < z_min) | (data[i]['z'] > z_max)):
            outlier_indices.append(i)

    print('Indices where there are outliers:')
    for i in range(len(outlier_indices)):
        print(outlier_indices[i])

    new_data = np.delete(data, outlier_indices, 0)
    return new_data


# Shifting
def _shift(data, x_shift_factor, y_shift_factor, z_shift_factor, box_dimensions):
    for i in range(len(data)):
        data[i]['x'] = round(data[i]['x'] - x_shift_factor, 3)
        data[i]['y'] = round(data[i]['y'] - y_shift_factor, 3)
        data[i]['z'] = round(data[i]['z'] - z_shift_factor, 3)

        # Account for periodic boundary conditions
        if data[i]['x'] > box_dimensions[0]:
            data[i]['x'] = round(data[i]['x'] - box_dimensions[0], 3)
        elif data[i]['x'] < 0:
            data[i]['x'] = data[i]['x'] + box_dimensions[0]

        if data[i]['y'] > box_dimensions[1]:
            data[i]['y'] = round(data[i]['y'] - box_dimensions[1], 3)
        elif data[i]['y'] < 0:
            data[i]['y'] = round(data[i]['y'] + box_dimensions[1], 3)

        if data[i]['z'] > box_dimensions[2]:
            data[i]['z'] = round(data[i]['z'] - box_dimensions[2], 3)
        elif data[i]['z'] < 0:
            data[i]['z'] = round(data[i]['z'] + box_dimensions[2], 3)
    return data


def generate_new_model(label_map: Dict) -> None:
    images = []
    labels = []

    # Temporarily hard code file paths
    root_dir = os.getcwd()
    glob_path = os.path.join(root_dir, 'images', '*.png')
    files = glob.glob(glob_path)
    for myFile in files:
        print('Loading ', myFile)
        image = imageio.imread(myFile, pilmode="RGB")
        images.append(image)
        label = myFile.split('/images/')[-1].split('_')[0][:-1]
        print(label)
        labels.append(label_map[label])

    print('Images loaded.')
    train_images, test_images, train_labels, test_labels = train_test_split(images, labels, test_size=0.25,
                                                                            random_state=1)

    train_images = np.array(train_images)
    test_images = np.array(test_images)
    train_labels = np.array(train_labels)
    test_labels = np.array(test_labels)

    model = models.Sequential()
    model.add(layers.Conv2D(32, (3, 3), activation='relu', input_shape=(480, 640, 3)))
    model.add(layers.MaxPooling2D((2, 2)))
    model.add(layers.Conv2D(64, (3, 3), activation='relu'))
    model.add(layers.MaxPooling2D((2, 2)))
    model.add(layers.Conv2D(64, (3, 3), activation='relu'))

    model.add(layers.Flatten())
    model.add(layers.Dense(64, activation='relu'))
    model.add(layers.Dense(10))

    # Compile keras model
    model.compile(optimizer='adam',
                  loss=tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True),
                  metrics=['accuracy'])

    # Fit model
    model.fit(train_images, train_labels, epochs=5,
                        validation_data=(test_images, test_labels))

    # Evaluate model
    test_loss, test_acc = model.evaluate(test_images, test_labels, verbose=2)

    # Save model
    os.mkdir('model')
    model.save('model/trained_cnn_model')
