#!/usr/bin/env python

import imageio.v2 as imageio
import numpy as np
import os
import sys
import tensorflow as tf

from peptide_classification_utilities import create_image, generate_images, generate_new_model

MODEL_PATHNAME='/home/mason/pace2/PACE2/examples/peptide_example/model/trained_cnn_model'

def main():
    label_map = {'vesicle': 0, 'nanotube': 1, 'lamellae': 2, 0: 'vesicle', 1: 'nanotube', 2: 'lamellae'}

    # Check if model has already been saved. If not, train a model.
    if not os.path.exists(MODEL_PATHNAME):
        generate_images()
        print('Training and generating new model.')
        generate_new_model(label_map)
    model = tf.keras.models.load_model(MODEL_PATHNAME)

    # Predict using the model.
    image_path = create_image(sys.argv[1], True)
    image = imageio.imread(image_path, pilmode="RGB")
    prediction = np.argmax(model.predict(np.array([image])))
    print(prediction)

    # Write the prediction to a file for data collection.
    with open('output.dat', 'w') as f:
        f.write(label_map[prediction])
    
    # Write the extension signal file.
    with open('out.txt','w') as f:
        if prediction in [0,1,2]:
            f.write(str(1))
        else: f.write(str(0))

if __name__ == '__main__':
    main()
