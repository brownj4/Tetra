"""
Copyright (c) 2016 Julian Brown

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

import numpy as np
import itertools
import shelve
import os
from PIL import Image
import scipy.ndimage
import scipy.optimize
import scipy.stats
import glob

# directory containing input images
image_directory = './pics'

# boolean for whether or not to display an annotated version
# of the image with identified stars circled in green and
# unmatched catalog stars circled in red
show_solution = True

# maximum fields of view of catalog patterns in degrees
# determines approximately what image fields of view
# can be identified by the algorithm
max_fovs = [10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]

# radius around image stars to search for matching catalog stars
# as a fraction of image field of view in x dimension
match_radius = .01

# maximum error in edge ratios
max_error = .01

# image downsampling factor for median filtering
downsample_factor = 2

# median filter window radius in down-sampled image pixels
filter_radius = 2
filter_width = filter_radius * 2 + 1

# percentage of pattern catalog that stores values
catalog_fill_factor = .5

# number of divisions along each dimension of the pattern catalog
num_catalog_bins = 25

# maximum number of pattern catalog stars within
# the maximum fov centered on any given pattern catalog star
max_stars_per_fov = 10

# minimum star brightness magnitude
magnitude_minimum = 5.0

# minimum allowable angle between star vectors in radians
# used to remove double stars
min_angle = .002

# number of stars to use in each pattern
pattern_size = 5

# minimum number of pixels in a group of bright pixels
# needed to classify the group as a star
min_pixels_in_group = 3

# centroiding window radius around a star's center pixel
# does not count the center pixel
window_radius = 2

# maximum number of bright stars to check against pattern catalog
max_pattern_checking_stars = 8

# maximum probability of mismatch for verifying an attitude determination
max_mismatch_probability = .00001

# percentage of fine sky map that stores values
fine_sky_map_fill_factor = .5

# number of divisions to break a single radius of
# the celestial sphere into for rapid star lookup
num_fine_sky_map_bins = 100

# percentage of course sky map that stores values
course_sky_map_fill_factor = .5

# number of divisions to break a single radius of
# the celestial sphere into for rapid star lookup
num_course_sky_map_bins = 4

# constant used for randomizing hash functions
avalanche_constant = 2654435761
  
# converts a hash_code into an index in the hash table
def hash_code_to_index(hash_code, bins_per_dimension, hash_table_size):
  # convert hashcode to python integers
  hash_code = [int(value) for value in hash_code]
  # represent hash code as a single integer
  integer_hash_code = sum(hash_code[i] * bins_per_dimension ** i for i in range(len(hash_code)))
  # randomize the hash code by multiplying by the avalanching constant
  # take the result modulo the table size to give a random index
  index = (integer_hash_code * avalanche_constant) % hash_table_size
  return index
  
# find all stars within a radius centered on the given vector using the compressed course sky map
def get_nearby_stars_compressed_course(vector, radius):
  # create list of nearby stars
  nearby_star_ids = []
  # given error of at most radius in each dimension, compute the space of hash codes to lookup in the sky map
  hash_code_space = [range(max(low,0), min(high+1,2*num_course_sky_map_bins)) for (low, high) in zip(((vector + 1 - radius) * num_course_sky_map_bins).astype(np.int),
                                                                                                     ((vector + 1 + radius) * num_course_sky_map_bins).astype(np.int))]
  # iterate over hash code space, looking up partitions of the sky map that are within range of the given vector
  for hash_code in itertools.product(*hash_code_space):
    hash_index = hash_code_to_index(hash_code, 2*num_course_sky_map_bins, compressed_course_sky_map_hash_table_size)
    # iterate over the star lists in the given partition, adding them to
    # the nearby stars list if they're in the correct bin and within range of the vector
    for index in ((2 * (hash_index + offset ** 2)) % compressed_course_sky_map_hash_table_size for offset in itertools.count()):
      # if the current slot is empty, the bin does not exist
      if not compressed_course_sky_map[index]:
        break
      # otherwise, check if the indices correspond to the correct bin
      indices = compressed_course_sky_map[index:index+2]
      # extract the sublist of star ids given by the indices
      star_id_list = compressed_course_sky_map[slice(*indices)]
      # check if the hash code for the first star matches the bin's hash code
      first_star_vector = star_table[star_id_list[0]]
      first_star_hash_code = tuple(((first_star_vector+1)*num_course_sky_map_bins).astype(np.int))
      if first_star_hash_code == hash_code:
        # iterate over the stars in the sublist, adding them to
        # the nearby stars list if they're within range of the vector
        for star_id in star_id_list:
          if np.dot(vector, star_table[star_id]) > np.cos(radius):
            nearby_star_ids.append(star_id)
  return nearby_star_ids

# open the pattern catalog and fine sky map and test whether they are fully
# generated with the following parameters and if not, regenerate them
parameters = (max_fovs, 
              num_catalog_bins, 
              max_stars_per_fov, 
              magnitude_minimum, 
              min_angle, 
              pattern_size, 
              fine_sky_map_fill_factor, 
              num_fine_sky_map_bins,
              course_sky_map_fill_factor, 
              num_course_sky_map_bins)
# try opening the database files
try:
  pattern_catalog = np.load('pattern_catalog.npy')
  fine_sky_map = np.load('fine_sky_map.npy')
  compressed_course_sky_map = np.load('compressed_course_sky_map.npy')
  compressed_course_sky_map_hash_table_size = compressed_course_sky_map[-1]
  star_table = np.load('star_table.npy')
  stored_parameters = open('params.txt', 'r').read()
  # if it got this far, the reads didn't fail
  read_failed = 0
except:
  # loading from the files failed, so they either didn't existing or weren't complete
  read_failed = 1
# there are two cases in which the catalog needs to be regenerated:
# reading the stored files failed or the stored parameters
# are different than those specified above
if read_failed or str(parameters) != stored_parameters:
  # number of stars in BSC5 catalog
  STARN = 9110

  # BSC5 data storage format
  bsc5_data_type = [("XNO", np.float32),
                    ("SRA0", np.float64),
                    ("SDEC0", np.float64),
                    ("IS", np.int16),
                    ("MAG", np.int16),
                    ("XRPM", np.float32),
                    ("XDPM", np.float32)
                   ]
  
  # open BSC5 catalog file for reading
  bsc5_file = open('BSC5', 'rb')
  # skip first 28 header bytes
  bsc5_file.seek(28)
  # read BSC5 catalog into array
  bsc5 = np.fromfile(bsc5_file, dtype=bsc5_data_type, count=STARN)

  # year to calculate catalog for
  # should be relatively close to 1950
  year = 2016

  # retrieve star positions, magnitudes and ids from BSC5 catalog
  stars = []
  for star_num in range(STARN):
    # only use stars brighter (i.e. lower magnitude)
    # than the minimum allowable magnitude
    mag = bsc5[star_num][4] / 100.0
    if mag <= magnitude_minimum:
      # retrieve RA in 1950
      ra = bsc5[star_num][1]
      # correct RA to modern day
      ra += bsc5[star_num][5] * (year - 1950)
      # retrieve DEC in 1950
      dec = bsc5[star_num][2]
      # correct DEC to modern day
      dec += bsc5[star_num][6] * (year - 1950)
      # skip blank star entries
      if ra == 0.0 and dec == 0.0:
        continue
      # convert RA, DEC to (x,y,z)
      vector = np.array([np.cos(ra)*np.cos(dec), np.sin(ra)*np.cos(dec), np.sin(dec)])
      # retrieve star ID number in BSC5
      star_id = int(bsc5[star_num][0])
      # add vector, magnitude pair to list of stars
      stars.append((vector, mag, star_id))
        
  # fast method for removing double stars using sort
  # sort list by x component of star vectors
  stars.sort(key=lambda star: star[0][0])
  # create boolean indicator list for which star vectors are double stars
  doubles = [0] * len(stars)
  for star_num1 in range(len(stars)):
    for star_num2 in xrange(star_num1 + 1, len(stars)):
      # skip checking star pairs with x component differences which
      # are already larger than the minimum allowable angle
      if stars[star_num2][0][0] - stars[star_num1][0][0] >= min_angle:
        break
      # if the star pair forms a double star, add both to the indicator list
      if np.dot(stars[star_num1][0], stars[star_num2][0]) > np.cos(min_angle):
        doubles[star_num1] = 1
        doubles[star_num2] = 1
        break
  # use the double star indicator list to create a new star vectors list without double stars
  stars_no_doubles = [stars[i] for i in range(len(stars)) if not doubles[i]]

  # output how many stars will be stored in the star table and the fine and course sky maps
  print("number of stars in star table and sky maps: " + str(len(stars)))
  
  # add all non-double stars brighter than magnitude_minimum to the star table and sky maps
  star_table = np.zeros((STARN+1, 3), dtype=np.float32)
  # create fine sky map hash table, which maps vectors to star ids
  fine_sky_map = np.zeros(len(stars) / fine_sky_map_fill_factor, dtype=np.uint16)
  # create course sky map hash table, which maps vectors to star ids
  course_sky_map = {}
  for (vector, mag, star_id) in stars_no_doubles:
    # add star vector to star table in position corresponding to its id
    star_table[star_id] = vector
    # find which partition the star occupies in the fine sky map hash table
    hash_code = ((vector+1)*num_fine_sky_map_bins).astype(np.int)
    hash_index = hash_code_to_index(hash_code, 2*num_fine_sky_map_bins, fine_sky_map.size)
    # use quadratic probing to find an open space in the fine sky map to insert the star in
    for index in ((hash_index + offset ** 2) % fine_sky_map.size for offset in itertools.count()):
      # if the current slot is empty, add the star
      # otherwise, move on to the next slot
      if not fine_sky_map[index]:
        fine_sky_map[index] = star_id
        break
    # find which partition the star occupies in the course sky map hash table
    hash_code = tuple(((vector+1)*num_course_sky_map_bins).astype(np.int))
    # if the partition is empty, create a new list to hold the star
    # if the partition already contains stars, add the star to the list
    course_sky_map[hash_code] = course_sky_map.pop(hash_code, []) + [star_id]

  # create compressed version of course sky map by indirectly mapping vectors to star ids
  # the map consists of a hash table, a superlist of stars, and a number representing the size of the hash table
  # the hash table consists of pairs of indices which slice the superlist into the output star id lists
  compressed_course_sky_map_hash_table_size = 2 * len(course_sky_map.keys()) / course_sky_map_fill_factor
  compressed_course_sky_map = np.zeros(compressed_course_sky_map_hash_table_size + len(stars_no_doubles) + 1, dtype=np.uint16)
  compressed_course_sky_map[-1] = compressed_course_sky_map_hash_table_size
  # add the items of the course sky map to the compressed course sky map one at a time
  first_open_slot_in_superlist = compressed_course_sky_map_hash_table_size
  for (hash_code, star_id_list) in course_sky_map.items():
    # compute the indices for the slice of the superlist the star id list will occupy
    slice_indices = (first_open_slot_in_superlist, first_open_slot_in_superlist + len(star_id_list))
    # add the star id list to the superlist
    compressed_course_sky_map[slice(*slice_indices)] = star_id_list
    # increment the counter for the first open slot in the superlist
    first_open_slot_in_superlist += len(star_id_list)
    hash_index = hash_code_to_index(hash_code, 2*num_course_sky_map_bins, compressed_course_sky_map_hash_table_size)
    # use quadratic probing to find an open space in the hash table to insert the star in
    for index in ((2 * (hash_index + offset ** 2)) % compressed_course_sky_map_hash_table_size for offset in itertools.count()):
      # if the current slot is empty, add the slice indices to the hash table
      # otherwise, move on to the next slot
      if not compressed_course_sky_map[index]:
        compressed_course_sky_map[index:index+2] = slice_indices
        break
  
  # sort list by star magnitude, from brightest to dimmest
  stars_no_doubles.sort(key=lambda star: star[1])

  # find all stars within a radius centered on the given vector using the pruned course sky map
  def get_nearby_stars_pruned_course(vector, radius):
    # create list of nearby stars
    nearby_star_ids = []
    # given error of at most radius in each dimension, compute the space of hash codes to lookup in the sky map
    hash_code_space = [range(max(low,0), min(high+1,2*num_course_sky_map_bins)) for (low, high) in zip(((vector + 1 - radius) * num_course_sky_map_bins).astype(np.int),
                                                                                                       ((vector + 1 + radius) * num_course_sky_map_bins).astype(np.int))]
    # iterate over hash code space, looking up partitions of the sky map that are within range of the given vector
    for hash_code in itertools.product(*hash_code_space):
      # iterate over the stars in the given partition, adding them to
      # the nearby stars list if they're within range of the vector
      for star_id in pruned_course_sky_map.get(hash_code, []):
        if np.dot(vector, star_table[star_id]) > np.cos(radius):
          nearby_star_ids.append(star_id)
    return nearby_star_ids
  
  # generate pattern catalog
  print("generating catalog, this may take an hour...")
  # create temporary list to store the patterns
  pattern_list = np.zeros((100000000, pattern_size), dtype=np.uint16)
  # create counter, which records how many patterns have been created
  num_patterns_found = 0

  # generate a piece of the catalog for each fov specified
  for max_fov in max_fovs:
    print("computing " + str(max_fov) + " degree fov patterns...")
  
    # change field of view from degrees to radians
    max_fov_rad = max_fov * np.pi / 180
    
    # fast method for pruning high density areas of the sky
    # create a hash table of the sky, which divides the
    # unit cube around the celestial sphere
    # up into (2*num_course_sky_map_bins)^3 partitions
    pruned_course_sky_map = {}
    # insert the stars into the hash table
    for (vector, mag, star_id) in stars_no_doubles:
      # skip stars that have too many closeby, brighter stars
      if len(get_nearby_stars_pruned_course(vector, max_fov_rad / 2)) >= max_stars_per_fov:
        continue
      # find which partition the star occupies in the hash table
      hash_code = tuple(((vector+1)*num_course_sky_map_bins).astype(np.int))
      # if the partition is empty, create a new list to hold the star
      # if the partition already contains stars, add the star to the list
      pruned_course_sky_map[hash_code] = pruned_course_sky_map.pop(hash_code, []) + [star_id]
    # create a list of stars without high density areas of the sky
    star_ids_pruned = [star_id for sublist in pruned_course_sky_map.values() for star_id in sublist]
      
    # initialize pattern, which will contain pattern_size star ids
    pattern = [None] * pattern_size
    for pattern[0] in star_ids_pruned:
      # find which partition the star occupies in the sky hash table
      hash_code = tuple(np.floor((star_table[pattern[0]]+1)*num_course_sky_map_bins).astype(np.int))
      # remove the star from the sky hash table
      pruned_course_sky_map[hash_code].remove(pattern[0])
      # iterate over all possible patterns containing the removed star
      for pattern[1:] in itertools.combinations(get_nearby_stars_pruned_course(star_table[pattern[0]], max_fov_rad), pattern_size-1):
        # retrieve the vectors of the stars in the pattern
        vectors = [star_table[star_id] for star_id in pattern]
        # verify that the pattern fits within the maximum field-of-view
        # by checking the distances between every pair of stars in the pattern
        if all(np.dot(*star_pair) > np.cos(max_fov_rad) for star_pair in itertools.combinations(vectors[1:], 2)):
          pattern_list[num_patterns_found] = pattern
          num_patterns_found += 1
  # truncate pattern list to only contain valid values
  pattern_list = pattern_list[:num_patterns_found]

  # insert star patterns into pattern catalog hash table
  print("inserting patterns into catalog...")
  pattern_catalog = np.zeros((num_patterns_found / catalog_fill_factor, pattern_size), dtype=np.uint16)
  for pattern in pattern_list:
    # retrieve the vectors of the stars in the pattern
    vectors = np.array([star_table[star_id] for star_id in pattern])
    # calculate and sort the edges of the star pattern, which are the distances between its stars
    edges = np.sort([np.linalg.norm(np.subtract(*star_pair)) for star_pair in itertools.combinations(vectors, 2)])
    # extract the largest edge
    largest_edge = edges[-1]
    # divide the edges by the largest edge to create dimensionless ratios
    edge_ratios = edges[:-1] / largest_edge
    # convert edge ratio float to hash code by binning
    hash_code = tuple((edge_ratios * num_catalog_bins).astype(np.int))
    hash_index = hash_code_to_index(hash_code, num_catalog_bins, pattern_catalog.shape[0])
    # use quadratic probing to find an open space in the pattern catalog to insert the pattern in
    for index in ((hash_index + offset ** 2) % pattern_catalog.shape[0] for offset in itertools.count()):
      # if the current slot is empty, add the pattern
      if not pattern_catalog[index][0]:
        pattern_catalog[index] = pattern
        break
      # if the current slot contains a previously inserted
      # copy of the same pattern, don't add the pattern
      elif sorted(pattern_catalog[index]) == sorted(pattern):
        break
      # otherwise, continue the search by moving on to the next slot
      else:
        continue

  # save star table, sky maps, pattern catalog, and parameters to disk
  np.save('star_table.npy', star_table)
  np.save('fine_sky_map.npy', fine_sky_map)
  np.save('compressed_course_sky_map.npy', compressed_course_sky_map)
  np.save('pattern_catalog.npy', pattern_catalog)
  parameters = open('params.txt', 'w').write(str(parameters))
  
# run the tetra star tracking algorithm on the given image
def tetra(image_file_name):
  # read image from file and convert to black and white
  image = np.array(Image.open(image_file_name).convert('L'))
  # extract height (y) and width (x) of image
  height, width = image.shape

  # flatten image by subtracting median filtered image
  # clip image so size is a multiple of downsample_factor
  # note that this may shift the center of the image
  height = height - height % downsample_factor
  width = width - width % downsample_factor
  image = image[:height, :width]
  # downsample image for median filtering
  downsampled_image = image.reshape((height//downsample_factor,downsample_factor,width//downsample_factor,downsample_factor)).mean(axis=3).mean(axis=1)
  # apply median filter to downsampled image
  median_filtered_image = scipy.ndimage.filters.median_filter(downsampled_image, size=filter_width, output=image.dtype)
  # upsample median filtered image back to original image size
  upsampled_median_filtered_image = median_filtered_image.repeat(downsample_factor, axis=0).repeat(downsample_factor, axis=1)
  # subtract the minimum of the image pixel and the local median to prevent values less than 0
  normalized_image = image - np.minimum.reduce([upsampled_median_filtered_image, image])

  # find all groups of pixels brighter than 5 sigma
  bright_pixels = zip(*np.where(normalized_image > 5 * np.std(normalized_image)))
  # group adjacent bright pixels together
  # create a dictionary mapping pixels to their group
  pixel_to_group = {}
  # iterate over the pixels from upper left to lower right
  for pixel in bright_pixels:
    # check whether the pixels above or to the left are part of
    # an existing group, which the current pixel will be added to
    left_pixel = (pixel[0]  , pixel[1]-1)
    up_pixel   = (pixel[0]-1, pixel[1]  )
    in_left_group = left_pixel in pixel_to_group
    in_up_group = up_pixel in pixel_to_group
    # if both are part of existing, disjoint groups, add the current pixel and combine the groups
    if in_left_group and in_up_group and id(pixel_to_group[left_pixel]) != id(pixel_to_group[up_pixel]):
      # add the current pixel to the upper pixel's group
      pixel_to_group[up_pixel].append(pixel)
      # append the upper pixel group onto the left pixel group
      pixel_to_group[left_pixel].extend(pixel_to_group[up_pixel])
      # replace all of the upper pixel group's dictionary entries
      # with references to the left pixel group
      for up_group_pixel in pixel_to_group[up_pixel]:
        pixel_to_group[up_group_pixel] = pixel_to_group[left_pixel]
    # if exactly one of the left pixel or upper pixels is part of an existing group,
    # add the current pixel to that group and add the current pixel to the dictionary
    elif in_left_group:
      pixel_to_group[left_pixel].append(pixel)
      pixel_to_group[pixel] = pixel_to_group[left_pixel]
    elif in_up_group:
      pixel_to_group[up_pixel].append(pixel)
      pixel_to_group[pixel] = pixel_to_group[up_pixel]
    # if neither of the left pixel or upper pixel are in an existing group,
    # add the current pixel to its own group and store it in the dictionary
    else:
      pixel_to_group[pixel] = [pixel]
  # iterate over the dictionary to extract all of the unique groups
  seen = set()
  groups = [seen.add(id(group)) or group for group in pixel_to_group.values() if id(group) not in seen]

  # find the brightest pixel for each group containing at least
  # the minimum number of pixels required to be classified as a star
  star_center_pixels = [max(group, key=lambda pixel: normalized_image[pixel]) for group in groups if len(group) > min_pixels_in_group]
  # find the centroid, or center of mass, of each star
  window_size = window_radius * 2 + 1
  # pixel values are weighted by their distances from the left (x) and top (y) of the window
  x_weights = np.fromfunction(lambda y,x:x+.5,(window_size, window_size))
  y_weights = np.fromfunction(lambda y,x:y+.5,(window_size, window_size))
  star_centroids = []
  for (y,x) in star_center_pixels:
    # throw out star if it's too close to the edge of the image
    if y < window_radius or y >= height - window_radius or \
       x < window_radius or x >= width  - window_radius:
      continue
    # extract the window around the star center from the image
    star_window = normalized_image[y-window_radius:y+window_radius+1, x-window_radius:x+window_radius+1]
    # find the total mass, or brightness, of the window
    mass = np.sum(star_window)
    # calculate the center of mass of the window in the x and y dimensions separately
    x_center = np.sum(star_window * x_weights) / mass - window_radius
    y_center = np.sum(star_window * y_weights) / mass - window_radius
    # correct the star center position using the calculated center of mass to create a centroid
    star_centroids.append((y + y_center, x + x_center))
  # sort star centroids from brightest to dimmest by comparing the total masses of their window pixels
  star_centroids.sort(key=lambda (y,x):-np.sum(normalized_image[y-window_radius:y+window_radius+1, x-window_radius:x+window_radius+1]))

  # compute list of (i,j,k) vectors given list of (y,x) star centroids and
  # an estimate of the image's field-of-view in the x dimension
  # by applying the pinhole camera equations
  def compute_vectors(star_centroids, fov):
    center_x = width / 2.
    center_y = height / 2.
    fov_rad = fov * np.pi / 180
    scale_factor = np.tan(fov_rad / 2) / center_x
    star_vectors = []
    for (star_y, star_x) in star_centroids:
      j_over_i = (center_x - star_x) * scale_factor
      k_over_i = (center_y - star_y) * scale_factor
      i = 1. / np.sqrt(1 + j_over_i**2 + k_over_i**2)
      j = j_over_i * i
      k = k_over_i * i
      star_vectors.append(np.array([i,j,k]))
    return star_vectors

  # generate star patterns in order of brightness
  def centroid_pattern_generator(star_centroids, pattern_size):
    # break if there aren't enough centroids to make even one pattern
    if len(star_centroids) < pattern_size:
      return
    star_centroids = np.array(star_centroids)
    # create a list of the pattern's centroid indices
    # add the lower and upper index bounds as the first
    # and last elements, respectively
    pattern_indices = [-1] + range(pattern_size) + [len(star_centroids)]
    # output the very brightest centroids before doing anything else
    yield star_centroids[pattern_indices[1:-1]]
    # iterate until the very dimmest centroids have been output
    # which occurs when the first pattern index has reached its maximum value
    while pattern_indices[1] < len(star_centroids) - pattern_size:
      # increment the pattern indices in order
      for index_to_change in range(1, pattern_size + 1):
        pattern_indices[index_to_change] += 1
        # if the current set of pattern indices is valid, use them
        if pattern_indices[index_to_change] < pattern_indices[index_to_change + 1]:
          break
        # otherwise, incrementing caused a conflict with the next pattern index
        # resolve the conflict by resetting the current pattern index and moving on
        else:
          pattern_indices[index_to_change] = pattern_indices[index_to_change - 1] + 1
      # output the centroids corresponding to the current set of pattern indices
      yield star_centroids[pattern_indices[1:-1]]
          
  # iterate over every combination of size pattern_size of the brightest max_pattern_checking_stars stars in the image
  for pattern_star_centroids in centroid_pattern_generator(star_centroids[:max_pattern_checking_stars], pattern_size):
    # iterate over possible fields-of-view
    for fov_estimate in max_fovs:
      # compute star vectors using an estimate for the field-of-view in the x dimension
      pattern_star_vectors = compute_vectors(pattern_star_centroids, fov_estimate)
      # calculate and sort the edges of the star pattern, which are the Euclidean distances between its stars' vectors
      pattern_edges = np.sort([np.linalg.norm(np.subtract(*star_pair)) for star_pair in itertools.combinations(pattern_star_vectors, 2)])
      # extract the largest edge
      pattern_largest_edge = pattern_edges[-1]
      # divide the pattern's edges by the largest edge to create dimensionless ratios for lookup in the catalog
      pattern_edge_ratios = pattern_edges[:-1] / pattern_largest_edge
      # given error of at most max_error in the edge_ratios, compute the space of hash codes to lookup in the catalog
      hash_code_space = [range(max(low,0), min(high+1,num_catalog_bins)) for (low, high) in zip(((pattern_edge_ratios - max_error) * num_catalog_bins).astype(np.int),
                                                                                                ((pattern_edge_ratios + max_error) * num_catalog_bins).astype(np.int))]
      # iterate over hash code space, only looking up non-duplicate codes that are in sorted order
      for hash_code in set([tuple(sorted(code)) for code in itertools.product(*hash_code_space)]):
        hash_code = tuple(hash_code)
        hash_index = hash_code_to_index(hash_code, num_catalog_bins, pattern_catalog.shape[0])
        # use quadratic probing to find all slots that patterns with the given hash code could appear in
        for index in ((hash_index + offset ** 2) % pattern_catalog.shape[0] for offset in itertools.count()):
          # if the current slot is empty, we've already
          # seen all patterns that match the given hash code
          if not pattern_catalog[index][0]:
            break
          # retrieve the star ids of possible match from pattern catalog
          catalog_pattern = pattern_catalog[index]
          # retrieve the vectors of the stars in the catalog pattern
          catalog_vectors = np.array([star_table[star_id] for star_id in catalog_pattern])
          # find the centroid, or average position, of the star pattern
          centroid = np.mean(catalog_vectors, axis=0)
          # calculate each star's radius, or Euclidean distance from the centroid
          radii = [np.linalg.norm(vector - centroid) for vector in catalog_vectors]
          # use the radii to uniquely order the catalog vectors
          catalog_sorted_vectors = catalog_vectors[np.argsort(radii)]
          # calculate and sort the edges of the star pattern, which are the distances between its stars
          catalog_edges = np.sort([np.linalg.norm(np.subtract(*star_pair)) for star_pair in itertools.combinations(catalog_vectors, 2)])
          # extract the largest edge
          catalog_largest_edge = catalog_edges[-1]
          # divide the edges by the largest edge to create dimensionless ratios
          catalog_edge_ratios = catalog_edges[:-1] / catalog_largest_edge
          # verify star patterns match to within the given maximum allowable error
          # note that this also filters out star patterns from colliding bins
          if any([abs(val) > max_error for val in (catalog_edge_ratios - pattern_edge_ratios)]):
            continue
          # compute the actual field-of-view using least squares optimization
          # compute the catalog pattern's edges for error estimation
          catalog_edges = np.append(catalog_edge_ratios * catalog_largest_edge, catalog_largest_edge)
          # helper function that calculates a list of errors in pattern edge lengths
          # with the catalog edge lengths for a given fov
          def fov_to_error(fov):
            # recalculate the pattern's star vectors given the new fov
            pattern_star_vectors = compute_vectors(pattern_star_centroids, fov)
            # recalculate the pattern's edge lengths
            pattern_edges = np.sort([np.linalg.norm(np.subtract(*star_pair)) for star_pair in itertools.combinations(pattern_star_vectors, 2)])
            # return a list of errors, one for each edge
            return catalog_edges - pattern_edges
          # find the fov that minimizes the squared error, starting with the given estimate
          fov = scipy.optimize.leastsq(fov_to_error, fov_estimate)[0][0]
          # convert newly computed fov to radians
          fov_rad = fov * np.pi / 180
          # find half diagonal fov of image in radians
          fov_half_diagonal_rad = fov_rad * np.sqrt(width ** 2 + height ** 2) / (2 * width)
          # recalculate star vectors using the new field-of-view
          pattern_star_vectors = compute_vectors(pattern_star_centroids, fov)
          # find the centroid, or average position, of the star pattern
          pattern_centroid = np.mean(pattern_star_vectors, axis=0)
          # calculate each star's radius, or Euclidean distance from the centroid
          pattern_radii = [np.linalg.norm(star_vector - pattern_centroid) for star_vector in pattern_star_vectors]
          # use the radii to uniquely order the pattern's star vectors so they can be matched with the catalog vectors
          pattern_sorted_vectors = np.array(pattern_star_vectors)[np.argsort(pattern_radii)]
          
          # calculate the least-squares rotation matrix from the catalog frame to the image frame
          def find_rotation_matrix(image_vectors, catalog_vectors):
            # find the covariance matrix H between the image vectors and catalog vectors
            H = np.sum([np.dot(image_vectors[i].reshape((3,1)), catalog_vectors[i].reshape((1,3))) for i in range(len(image_vectors))], axis=0)
            # use singular value decomposition to find the rotation matrix
            U, S, V = np.linalg.svd(H)
            rotation_matrix = np.dot(U, V)
            # correct reflection matrix if determinant is -1 instead of 1
            # by flipping the sign of the third column of the rotation matrix
            rotation_matrix[:,2] *= np.linalg.det(rotation_matrix)
            return rotation_matrix
          
          # use the pattern match to find an estimate for the image's rotation matrix
          rotation_matrix = find_rotation_matrix(pattern_sorted_vectors, catalog_sorted_vectors)
          # calculate all star vectors using the new field-of-view
          all_star_vectors = compute_vectors(star_centroids, fov)
          
          def find_matches(all_star_vectors, rotation_matrix):
            # rotate each of the star vectors into the catalog frame by
            # using the inverse (transpose) of the tentative rotation matrix
            rotated_star_vectors = [np.dot(rotation_matrix.T, star_vector) for star_vector in all_star_vectors]
            # retrieve matching catalog vectors for each image vector
            catalog_vectors = []
            for rotated_star_vector in rotated_star_vectors:
              hash_code_space = [range(max(low,0), min(high+1,2*num_fine_sky_map_bins)) for (low, high) in zip(((rotated_star_vector + 1 - match_radius) * num_fine_sky_map_bins).astype(np.int),
                                                                                                               ((rotated_star_vector + 1 + match_radius) * num_fine_sky_map_bins).astype(np.int))]
              # iterate over hash code space, only looking up non-duplicate codes that are in sorted order
              matching_stars = []
              for hash_code in [code for code in itertools.product(*hash_code_space)]:
                hash_index = hash_code_to_index(hash_code, 2*num_fine_sky_map_bins, fine_sky_map.size)
                # use quadratic probing to find an open space in the fine sky map to insert the star in
                for index in ((hash_index + offset ** 2) % fine_sky_map.size for offset in itertools.count()):
                  # if the current slot is empty, all of the matching stars have been found
                  # otherwise, move on to the next slot
                  if not fine_sky_map[index]:
                    break
                  # only accept stars within the match radius
                  elif np.dot(star_table[fine_sky_map[index]], rotated_star_vector) > np.cos(match_radius * fov_rad):
                    matching_stars.append(star_table[fine_sky_map[index]])
              catalog_vectors.append(matching_stars)
            # stars must uniquely match a catalog star brighter than magnitude_minimum
            matches = [(image_vector, catalog_star[0]) for (image_vector, catalog_star) in zip(all_star_vectors, catalog_vectors) if len(catalog_star) == 1]
            return matches
          
          matches = find_matches(all_star_vectors, rotation_matrix)
          # calculate loose upper bound on probability of mismatch assuming random star distribution
          # find number of catalog stars appear in a circumscribed circle around the image
          image_center_vector = np.dot(rotation_matrix.T, np.array((1,0,0)))
          num_nearby_catalog_stars = len(get_nearby_stars_compressed_course(image_center_vector, fov_half_diagonal_rad))
          # calculate probability of a single random image centroid matching to a catalog star
          single_star_match_probability = num_nearby_catalog_stars * match_radius ** 2 * width / height
          # apply binomial theorem to calculate probability upper bound on this many mismatches
          mismatch_probability_upper_bound = 1 - scipy.stats.binom.cdf(len(matches)-1, num_nearby_catalog_stars, single_star_match_probability)
          # if a high probability match has been found, recompute the attitude using all matching stars
          if mismatch_probability_upper_bound < max_mismatch_probability:
            # recalculate the rotation matrix using the newly identified stars
            rotation_matrix = find_rotation_matrix(*zip(*matches))
            # recalculate matched stars given new rotation matrix
            matches = find_matches(all_star_vectors, rotation_matrix)
            # extract right ascension, declination, and roll from rotation matrix and convert to degrees
            ra = (np.arctan2(rotation_matrix[0][1], rotation_matrix[0][0]) % (2 * np.pi)) * 180 / np.pi
            dec = np.arctan2(rotation_matrix[0][2], np.sqrt(rotation_matrix[1][2]**2 + rotation_matrix[2][2]**2)) * 180 / np.pi
            roll = (np.arctan2(rotation_matrix[1][2], rotation_matrix[2][2]) % (2 * np.pi)) * 180 / np.pi
            # print out attitude and field-of-view to 4 decimal places
            print("RA:   %.4f" % ra)
            print("DEC:  %.4f" % dec)
            print("ROLL: %.4f" % roll)
            print("FOV:  %.4f" % fov)
            # display input image with green circles around matched catalog stars
            # and red circles around unmatched catalog stars
            if show_solution:
              # draws circles around where the given vectors appear in an image
              def draw_circles(image, vectors, color, circle_fidelity, circle_radius):
                # calculate the pixel position of the center of the image
                image_center_x = width / 2.
                image_center_y = height / 2.
                # calculate conversion ratio between pixels and distance in the unit celestial sphere
                scale_factor = image_center_x / np.tan(fov_rad / 2)
                # iterate over the vectors, adding a circle for each one that appears in the image frame
                for (i, j, k) in vectors:
                  # find the center pixel for the vector's circle
                  circle_center_x = np.floor(image_center_x - (j / i) * scale_factor)
                  circle_center_y = np.floor(image_center_y - (k / i) * scale_factor)
                  # draw a circle of the given color with the given fidelity
                  for angle in np.array(range(circle_fidelity)) * 2 * np.pi / circle_fidelity:
                    # find the x and y coordinates for the pixel that will be drawn
                    pixel_x = int(circle_center_x + circle_radius * np.sin(angle))
                    pixel_y = int(circle_center_y + circle_radius * np.cos(angle))
                    # verify the pixel is within the image bounds
                    if pixel_x < 0 or pixel_x >= width or pixel_y < 0 or pixel_y >= height:
                      continue
                    # draw the pixel
                    image.putpixel((pixel_x, pixel_y), color)
              # plot the image with green circles around matched stars
              # and red circles around stars that weren't matched
              rgb_image = Image.fromarray(image).convert('RGB')
              # the circle is drawn using the corners of an n-gon, where the circle fidelity is n
              circle_fidelity = 100
              # star centroids that appear within the circle radius would match with the circle's corresponding catalog vector
              circle_radius = match_radius * width + 1
              # find which catalog stars could appear in the image
              image_center_vector = np.dot(rotation_matrix.T, np.array((1,0,0)))
              nearby_catalog_stars = get_nearby_stars_compressed_course(image_center_vector, fov_half_diagonal_rad)
              # rotate the vectors of all of the nearby catalog stars into the image frame
              rotated_nearby_catalog_vectors = [np.dot(rotation_matrix, star_table[star_id]) for star_id in nearby_catalog_stars]
              # color all of the circles red by default
              color_all = (255, 0, 0)
              draw_circles(rgb_image, rotated_nearby_catalog_vectors, color_all, circle_fidelity, circle_radius)
              # rotate the matched catalog stars into the image frame
              matched_rotated_catalog_vectors = [np.dot(rotation_matrix, catalog_vector) for (image_vector, catalog_vector) in matches]
              # recolor matched circles green
              color_matched = (0, 255, 0)
              draw_circles(rgb_image, matched_rotated_catalog_vectors, color_matched, circle_fidelity, circle_radius)
              rgb_image.show()
            return

  # print failure message
  print("failed to determine attitude")

for image_file_name in glob.glob(image_directory + '/*'):
  print image_file_name
  tetra(image_file_name)
