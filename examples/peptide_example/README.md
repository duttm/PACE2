# PACE2 Peptide Example

This example consists of the peptide assembly system consisting of two Candidates. In the analysis phase, a Convolutional Neural Network (CNN) image classification model is used. Code for the CNN model is included in this repository. Details of the original classifier code including authorship are below. See the Setup section there for additional Python requirements.

This example should run locally out of the box in a Python3.8+ environment. Config file templates are included here for execution on local host and remote HPC resources. Scripts associated with the CNN model have hardcoded filepaths in some places, so those should also be updated before running the example.

With a suitable Python environment, and after verifying filepaths are updated in all scripts, run the peptide example:
```
pace2.py simconfig_local.json resconfig_local.json
```

# Characterization and Analysis of Peptide Self-Assembilies

Peptide materials have various applications in biotechnology including drug delivery and electronics. Self-assembly of peptides has been shown to form ordered structures such as nanovesicles, nanotubes, and bilayers. With 20 amino acids, phase space for new molecules and thus, ordered structures, is vast. It is important to identify these molecules which yield materials with specific sequence-structure properties. This calls for a computational/ automated approach to peptide self-assembly structure characterization.

This analysis package uses a convolutional neural network-based machine learning approach to classifying peptide systems.

## Running the analysis

### Setup

This analysis package runs on Python 3. To setup your environment, please run the following commands:

```
python3 -m venv --system-site-packages ./venv
source ./venv/bin/activate
pip3 install numpy
pip3 install matplotlib
pip3 install imageio
pip3 install tensorflow
pip install -U scikit-learn
```

After cloning, all you need to do before running the analysis is to add your .gro test file to the working directory. 

### Analysis

Run peptide_classification_ml.py with your .gro test file as your first argument like so:

```
python3 peptide_classification_ml.py test_file.gro
```

This will continously generate two .dat files which are used to create the .png images used in the ML model. They will be overwritten in the process of generating the images.

```
peptideFile.dat
solventFile.dat
```

An `images` directory will also be created, where all of the training images will be stored. Using these images, the CNN model is trained and stored in the `model` directory. Note: this SAME analysis command as stated above can be used regardless of whether the model has already been trained or not. The script will automatically load the trained model if it already has been trained, and will create a new trained model if it has not. 

Using the specified gro file that the model is trying to predict for, the .gro file is converted to a .png image in the `images/predicting` directory. The analysis will  generate an `output.dat` file in the working directory containing a string corresponding to the three system types:
```
vesicle
nanotube
lamellae
```


## Authors

* **Ethan Zang** 

## Acknowledgments

* Srinivas Mushnoori, Kassandra Schmidt, Emran Reshid, Dr. Meenakshi Dutt
