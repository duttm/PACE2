import imageio
import numpy as np
import os
import sys
import tensorflow as tf

from peptide_classification_utilities import create_image, generate_images, generate_new_model


def main():
    label_map = {'vesicle': 0, 'nanotube': 1, 'lamellae': 2, 0: 'vesicle', 1: 'nanotube', 2: 'lamellae'}

    # Check if model has already been saved
    if not os.path.exists('model'):
        generate_images()
        print('Training and generating new model.')
        generate_new_model(label_map)
    model = tf.keras.models.load_model('model/trained_cnn_model')

    # Predict using the model
    image_path = create_image(sys.argv[1], True)
    image = imageio.imread(image_path, as_gray=False, pilmode="RGB")
    prediction = np.argmax(model.predict(np.array([image])))
    print(prediction)
    with open('output.dat', 'w') as f:
        f.write(label_map[prediction])


if __name__ == '__main__':
    main()
