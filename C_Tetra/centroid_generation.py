import PIL
from PIL import Image
import numpy as np
import cPickle as p
import itertools as iter
import random
import sys

np.random.seed(0)
random.seed(0)
	
# number of test images to generate
num_images = 1000

# camera FOV in degrees
# should be less than 180
# should also be less than maximum catalog angle
fov_y = 10

# image size in pixels
num_pixels_x = 1024
num_pixels_y = 1024

# maximum number of stars allowed per image
max_num_stars_per_image = 25

# number of stars in star vectors file
num_star_vectors = 5904

# name of star vectors file
star_vectors_file = "stars"

# calculate diagonal FOV in degrees from fov_y
fov_diag = np.arctan(np.tan(fov_y * np.pi / 180) * np.sqrt(num_pixels_x ** 2 + num_pixels_y ** 2) / num_pixels_y) * 180 / np.pi

# star vectors data format
star_vectors_data_type = [("i", np.float64),("j", np.float64),("k", np.float64),("mag", np.float64),("id", np.uint32),("pad", np.uint32)]

# import catalog vectors
catalog_vectors = np.memmap(star_vectors_file, dtype=star_vectors_data_type, mode="r", shape=(num_star_vectors,))

# find vectors of matched stars
matched_vectors = [np.array((catalog_vectors[star][0], catalog_vectors[star][1], catalog_vectors[star][2])) for star in range(num_star_vectors)]

# create star hash table for fast nearby star lookup
star_hash = {}
star_hash_max_dot = np.cos((np.pi / 360) * fov_diag + 2*np.arcsin(np.sqrt(3) / 40))
for x in range(20):
	for y in range(20):
		for z in range(20):
			icv = ((float(x) / 10.0) - .95, (float(y) / 10.0) - .95, (float(z) / 10.0) - .95)
			icv = icv / np.sqrt(np.dot(icv,icv))
			star_hash[(x,y,z)] = [vector for vector in matched_vectors if np.dot(icv,vector) > star_hash_max_dot]

# create reverse vector catalog
star_catalog = {}
for star in range(num_star_vectors):
	star_catalog[str(catalog_vectors[star][0]) + "," + str(catalog_vectors[star][1]) + "," + str(catalog_vectors[star][2])] = star

# create memmap of image data for fast processing with C
image_data_file = 'image_data' + '.p'
image_data = np.memmap(image_data_file, dtype=np.uint16, mode="w+", shape=(max_num_stars_per_image * num_images,))
image_data_index = 0

# create memmap of rotation matrices
pointing_data_file = 'pointing_data' + '.p'
pointing_data = np.memmap(pointing_data_file, dtype=np.float32, mode="w+", shape=(num_images, 3, 3))
pointing_data_index = 0

# create memmap of centroid data for fast processing with C
centroid_data_file = 'centroid_data' + '.p'
centroid_data = np.memmap(centroid_data_file, dtype=np.float32, mode="w+", shape=(max_num_stars_per_image * num_images, 2))
centroid_data_index = 0

# precompute some useful values
center_x = float(num_pixels_x) / 2
center_y = float(num_pixels_y) / 2
scale_factor = np.tan(fov_y * np.pi / 360) / center_y
max_dot = np.cos((np.pi / 360) * fov_diag)

for image_number in range(num_images):

	# # yaw, declination, roll of camera
	# yaw = (random.random() - .5) * 360
	# declination = (random.random() - .5) * 180
	# roll = (random.random() - .5) * 360

	# # convert declination to pitch
	# pitch = -declination

	# # convert to radians
	# yaw = yaw * np.pi / 180
	# pitch = pitch * np.pi / 180
	# roll = roll * np.pi / 180

	# # obtain rotation matrix of camera frame from yaw, pitch, roll
	# rotation_matrix = np.matrix([[np.cos(yaw)*np.cos(pitch),np.cos(yaw)*np.sin(pitch)*np.sin(roll)-np.sin(yaw)*np.cos(roll),np.cos(yaw)*np.sin(pitch)*np.cos(roll)+np.sin(yaw)*np.sin(roll)],[np.sin(yaw)*np.cos(pitch),np.sin(yaw)*np.sin(pitch)*np.sin(roll)+np.cos(yaw)*np.cos(roll),np.sin(yaw)*np.sin(pitch)*np.cos(roll)-np.cos(yaw)*np.sin(roll)],[-np.sin(pitch),np.cos(pitch)*np.sin(roll),np.cos(pitch)*np.cos(roll)]]).T

	ax = np.array([random.gauss(0,10), random.gauss(0,10), random.gauss(0,10)]);
	ax = ax / np.linalg.norm(ax)
	x = ax[0]
	y = ax[1]
	z = ax[2]
	
	ro = np.random.uniform(0,2*np.pi)
	c = np.cos(ro)
	s = np.sin(ro)
	
	rotation_matrix = np.matrix([[c+x*x*(1-c), x*y*(1-c)-z*s, x*z*(1-c)+y*s],
								 [y*x*(1-c)+z*s, c+y*y*(1-c), y*z*(1-c)-x*s],
								 [z*x*(1-c)-y*s, z*y*(1-c)+x*s, c+z*z*(1-c)]])
	
	# find center of image vector
	image_center_vector = np.dot(rotation_matrix.T, np.matrix((1,0,0)).T).T
	icv = image_center_vector

	# filter stars outside field of view
	#filtered_vectors = [vector for vector in matched_vectors if np.dot(image_center_vector,vector) > max_dot]
	filtered_vectors = [vector for vector in star_hash[(int((icv.item(0)+1)*10),int((icv.item(1)+1)*10),int((icv.item(2)+1)*10))] if np.dot(image_center_vector,vector) > max_dot]
	
	# randomize their ordering
	randomly_ordered_vectors = [filtered_vectors[i] for i in random.sample(range(len(filtered_vectors)),len(filtered_vectors))]
	
	# rotate all catalog vectors to camera frame
	rotated_matched_vectors = [np.dot(rotation_matrix, vector) for vector in randomly_ordered_vectors]
  
	# set write location
	centroid_data_index = max_num_stars_per_image * image_number
	
	# fill image with stars (add centroids to centroid data file)
	for star_id in range(len(rotated_matched_vectors[:max_num_stars_per_image])):
		star = rotated_matched_vectors[star_id]
		i = star.item(0)
		j = star.item(1)
		k = star.item(2)
		x = -(j / i) / scale_factor
		y = (k / i) / scale_factor
		if (abs(x) > num_pixels_x / 2 or abs(y) > num_pixels_y / 2):
			randomly_ordered_vectors[star_id] = (9,9,9)
			continue
		centroid_data[centroid_data_index][0] = x
		centroid_data[centroid_data_index][1] = y
		centroid_data_index += 1
	
	# set write location
	image_data_index = max_num_stars_per_image * image_number
	
	# add star ids to star id file in same order as centroids
	for star in randomly_ordered_vectors[:max_num_stars_per_image]:#filtered_vectors[:max_num_stars_per_image]:
		if (star[0] == 9):
			continue
		image_data[image_data_index] = star_catalog[(str(star[0]) + "," + str(star[1]) + "," + str(star[2]))]
		image_data_index += 1
		
	pointing_data[pointing_data_index] = rotation_matrix
	pointing_data_index += 1

	if image_number % 10 == 0:
		sys.stdout.write("\r" + str(image_number))
