#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

/* Number of stars in the stars file. */
#define num_stars 5904
/* Number of patterns locations cached when the catalog is accessed, specifying the */
/* maximum distance matching patterns can be from the original offset in catalog. */
/* Too small a value can cause false negatives.  Too large a value can cause an */
/* unnecessary increase in runtime.  Should be a triangular number plus one */
/* (i.e. 1,2,4,7,11,16,22,29,37,46,56,67,79,92...) as the catalog uses */
/* quadratic probing, so matches appear at triangular probe depths. */
/* Should be less than or equal to max_probe_depth. */
#define pattern_cache_size 16
/* Maximum distance matching pattern can be from original offset in catalog. */
/* Also the number of extra catalog positions past the end of the catalog, */
/* due to the catalog hash table being non-cyclic. */
#define max_probe_depth 4278
/* Number of stars in star pattern. */
/* Should be greater than 2.  Recommended number is 4. */
#define num_stars_in_pattern 4
/* Size of tetra catalog in Tetras. */
/* Hard coded for speed.  Must match generated catalog. */
#define catalog_size_in_patterns 770708495
/* Maximum Field of View for catalog in radians.  Must exactly */
/* match the max_fov used to generate the catalog . */
#define max_fov .247
/* Ratio of bin size to error range, where bins are the discretization of */
/* Pattern data values to allow them to be hashed and where the error range covers */
/* a range of twice the data value's maximum error.  When a data value's error range */
/* overlaps two bins, it's replicated into both.  By linearity of expectation, */
/* the expected number of replicas of a given Pattern is: */
/* (1+1/bin_size_ratio)^num_dims, where num_dims is the number of data values */
/* stored in Patterns.  The expected ratio of Patterns with matching bins to */
/* Patterns with matching values is: (1 + bin_size_ratio)^num_dims. */
/* The bin_size_ratio represents a tradeoff between catalog size and */
/* runtime, as more replicas means a larger catalog and more mismatching results */
/* means more time spent checking for a match.  In most cases, Patterns with */
/* matching values are rare, so a larger bin_size_ratio is better.  Must be */
/* greater than 1.0 without changes to Pattern struct, recommended value of 2.0. */
#define bin_size_ratio 3.0
/* Maximum star coordinate centroiding error as fraction of maximum FOV. */
/* .001 is .1% of the max FOV or 1.414 pixels in a 1000x1000 image. */
// 1 / (1024*sqrt(2)) < .00069054
#define max_centroid_error .00069054
/* Maximum error in imager FOV estimate as fraction of true imager FOV. */
/* Must be greater than zero to prevent division by zero errors. */
/* max_fov*(1+max_fov_error) must be less than pi, but should be much less. */
/* Should be less than or equal to the max_fov_error used to generate the catalog. */
/* .01 for a 10 degree FOV imager covers estimates from 9.9 to 10.1 degrees. */
#define max_fov_error 0.01
/* Maximum number of stars to process. */
#define max_stars 12
/* Maximum number of stars per image. */
#define max_stars_per_image 25
/* Number of test images. */
#define num_images 100000
/* Mathematical constant pi's approximate value. */
#define PI 3.1415926

/* The following values are not user defined constants and should not be changed. */
/* Maximum scaling of image caused by FOV error. */
float max_scale_factor = fmax(tan(max_fov*(1+max_fov_error)/2.0)/tan(max_fov/2.0),
                              1-tan(max_fov*(1-max_fov_error)/2.0)/tan(max_fov/2.0));
/* Largest edge error is proportional to the largest edge ratio by the image scale. */
/* For example, doubling the largest edge length doubles the scaling error. */
#define le_error_slope (max_scale_factor-1)
/* Largest edge error has a constant factor determined by the centroiding error. */
/* This value is the worst case of maximum FOV and centroid error.  It corresponds to */
/* The FOV taking on its minimum value and both centroids having maximal error inwards. */
#define le_error_offset (2*max_centroid_error/(2-max_scale_factor))
/* Largest possible largest_edge_length given the maximum FOV. */
#define max_le_length (2*sin(max_fov*(1+max_fov_error)/2.0))

/* Feature struct */
/* Data format of what's stored in a given feature of a catalog Pattern. */
/* A coordinate system is constructed by centering the largest edge along the x axis, */
/* where edges are straight lines between stars in the pattern. */
/* The x and y coordinates of every star are divided by the length of the largest edge. */
/* This places two stars at the (x,y) coordinates (-.5,0) and (.5,0) and the remaining */
/* stars in one of two places due to a 180 degree rotational ambiguity. */
/* Each Feature encodes each of the remaining stars' bins, coordinates, and id. */
typedef struct Feature
{
  /* Implicitly encoded doubles are used to store the coordinates. */
  /* This cuts down on disc usage by instead using integers with an */
  /* implicit division factor, as all coordinates are between -1 and 1. */
  int x : 15;
  /* Bins are identifying values for the Pattern's catalog position.  */
  /* They discretize Pattern data values, allowing them to be hashed. */
  /* A bin offset is the offset of this Pattern's bin from the */
  /* lowest bin the Pattern can occupy.  With a bin_size_ratio */
  /* greater than 1.0, bin offsets can be stored in a single bit. */
  unsigned int x_bin_offset : 1;
  int y : 15;
  unsigned int y_bin_offset : 1;
  /* ID of the Feature's star.  Unsigned 15 bit integers support 32768 unique stars. */
  unsigned int star_id : 15;
  /* Extra, unused bit. */
  unsigned int pad : 1;
} Feature;

/* Pattern struct */
/* Data format of what's stored in a given catalog position for a given star pattern. */
/* The rotation that results in the largest x bin (y bin as tie-breaker) */
/* is chosen to deal with the pattern's 180 degree rotational ambiguity. */
/* The largest edge length is also stored to allow for sanity-checking of FOV/image scale. */
typedef struct Pattern
{
  /* Features are stored in order of increasing x bin, with increasing y bin */
  /* as a tie breaker.  Two Features cannot contain the same bin pair, as this */
  /* would cause ambiguities while matching image stars to catalog stars. */
  Feature features[num_stars_in_pattern-2];
  /* Length of the largest edge in the star pattern.  Used to recompute the FOV as well as */
  /* sanity check matches.  Also stored as an implicitly encoded double between 0 and 1, */
  /* as edge lengths must be between 0 and max_le_length, the maximum largest edge length. */
  uint16_t largest_edge;
  /* As the largest edge length is always greater than zero, it is also used as */
  /* a boolean to indicate whether the hash table slot contains a star pattern. */
  /* Note that this is therefore not explicitly set during catalog creation. */
  #define has_pattern largest_edge
  unsigned int le_bin_offset : 1;
  unsigned int fixed_star_id_1 : 15;
  unsigned int fixed_star_id_2 : 15;
  /* Boolean value indicating whether this catalog position contains the last */
  /* matching pattern in the catalog.  This avoids having to probe to the next slot. */
  unsigned int is_last : 1;
} Pattern;

/* Star struct */
/* Data format of what's stored in the stars array for a given star. */
typedef struct Star
{
  /* Unit vector x,y,z coordinates pointing to the star from the center of the celestial sphere. */
  double vec[3];
  /* Star magnitude */
  double mag;
  /* Star id, also known as the Hipparchos number */
  unsigned int star_id;
} Star;

/* Euclidean difference evaluator for a pair of 3D vectors. */
/* Output is the vector from the end of i to the end of j. */
static void diff(double v1[3],
          double v2[3],
          double output[3]){
  output[0] = v1[0]-v2[0];
  output[1] = v1[1]-v2[1];
  output[2] = v1[2]-v2[2];
}

/* Vector magnitude evaluator for a 3D vector. */
static double mag(double v[3]){
  return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

/* Normalize a 3D vector by dividing by its magnitude. */
static void normalize(double v[3]){
  double v_magnitude = mag(v);
  v[0] /= v_magnitude;
  v[1] /= v_magnitude;
  v[2] /= v_magnitude;
}

/* Euclidean distance evaluator for a pair of 3D vectors. */
static double dist(double v1[3],
           double v2[3]){
  double diff_vector[3];
  diff(v1,v2,diff_vector);
  return mag(diff_vector);
}

/* Dot product evaluator for a pair of 3D vectors. */
static double dot_prod(double v1[3],
               double v2[3]){
  return v1[0]*v2[0]+ 
         v1[1]*v2[1]+ 
         v1[2]*v2[2];
}

/* Cross product magnitude evaluator for a pair of 3D vectors. */
static void cross_prod(double v1[3],
                double v2[3],
                double output[3]){
  output[0] = v1[1]*v2[2]-v1[2]*v2[1];
  output[1] = v1[2]*v2[0]-v1[0]*v2[2];
  output[2] = v1[0]*v2[1]-v1[1]*v2[0];
}

/* Calculate the exponent base for logarithmic binning. */
/* Also bounds check error values. */
static double get_base(double error_slope,
                      double error_offset){
  /* If the fixed error is non-positive, return error and exit. */
  if(error_offset <= 0){
    printf("\nNon-positive error value detected: increase error values.\n");
    exit(EXIT_FAILURE);
  }
  /* Calculate and return base of logarithmic binning function. */
  double base = (1+error_slope)/fmax(1-error_slope, 0);
  return base;
}

/* Logarithmically bin doubles that have error linear in their value.  Binning function */
/* is designed to ensure each value's error range overlaps at most two bins. */
static int log_bin(double input,
                   double error_slope,
                   double error_offset){
  /* If either of the error values is infinite, all values share the 0 bin. */
  if(!isfinite(error_slope) || !isfinite(error_offset)){
    return 0;
  }
  /* Calculate base of logarithmic binning function. */
  double base = get_base(error_slope, error_offset);
  int bin;
  /* Slopes much smaller than the offset result in linear binning. */
  if(base <= 1+error_offset*bin_size_ratio/10.){
    bin = input/(2*(error_slope+error_offset)*bin_size_ratio);
  }
  /* Otherwise, calculate logarithmic bin. */
  /* Slopes greater than or equal to one result in all values sharing the 0 bin. */
  else{
    bin = (log(input*error_slope/error_offset+1)/log(base))/bin_size_ratio;
  }
  return bin;
}

/* Retrieve minimum possible pre-binned value from a logarithmically binned value. */
static double log_unbin(int bin,
                       double error_slope,
                       double error_offset){
  /* If either of the error values is infinite, 0 is the minimum value. */
  if(!isfinite(error_slope) || !isfinite(error_offset)){
    return 0;
  }
  /* Calculate base of logarithmic binning function. */
  double base = get_base(error_slope, error_offset);
  double min_input;
  /* Slopes much smaller than the offset result in linear binning. */
  if(base <= 1+error_offset*bin_size_ratio/10.){
    min_input = bin*2*(error_slope+error_offset)*bin_size_ratio;
  }
  /* Otherwise, calculate minimum possible logarithmic pre-binned value. */
  else{
    min_input = (pow(base, bin*bin_size_ratio)-1)*error_offset/error_slope;
  }
  return min_input;
}

/* Evenly bin largest edge length values using the fact that the maximum error is linear */
/* with respect to length.  Logarithmic functions satisfy this bin size constraint. */
static int bin_largest_edge(unsigned int largest_edge,
                            int error_ratio){
  /* Convert largest edge length back to double between 0 and 1 from implicitly divided double. */
  /* Resulting value is a ratio of max_le_length, the maximum largest edge length. */
  double le_ratio = largest_edge/((1<<16)-1.0);
  /* Adjust error in le_ratio based on error_ratio between -1 and 1. */
  /* An error_ratio of -1 will give the lower bin, +1 will give the upper bin. */
  le_ratio += error_ratio*(le_ratio*le_error_slope+le_error_offset);
  /* Find and return largest edge bin using logarithmic binning function. */
  return log_bin(le_ratio, le_error_slope, le_error_offset);
}

/* Transform largest edge bin into a largest edge ratio. */
/* Returns the minimum largest edge ratio within the bin. */
static double unbin_largest_edge(unsigned int bin){
  /* Invert logarithmic binning function to retrieve minimum largest edge ratio. */
  double min_le_ratio = log_unbin(bin, le_error_slope, le_error_offset);
  return min_le_ratio;
}

/* Evenly bin coordinate using the largest edge bin to find the image scale. */
static int bin_y(int y,
                 unsigned int le_bin,
                 int error_ratio){
  /* Calculate the minimum size the largest edge ratio could have, given its bin. */
  double min_le_ratio = unbin_largest_edge(le_bin);
  /* Coordinate error is affected by both the largest edge error and the FOV error. */
  double error_constant = le_error_offset/(2-max_scale_factor);
  /* Coordinate error has a component proportional to the coordinate value. */
  /* This value is the worst case of maximum FOV and centroid error.  It corresponds */
  /* to the FOV taking on its minimum value, both the largest edge centroids having */
  /* maximal error in the same direction, and the coordinate centroid having maximal */
  /* error in the opposite direction of the centroiding errors of the largest edge. */
  double error_slope = error_constant/fmax(min_le_ratio-error_constant, 0);
  /* Coordinate error also has a constant factor determined by the centroiding error. */
  double error_offset = error_slope;
  /* Convert coordinate back to double between -1 and 1 from implicitly divided double. */
  /* Resulting value is a ratio of image_scale. */
  double y_ratio = y/((1<<14)-1.0);
  /* Adjust error in coordinate based on error_ratio between -1 and 1. */
  /* An error_ratio of -1 will give the lower bin, +1 will give the upper bin. */
  /* The error is inversely scaled by min_le_ratio, as the coordinate is a ratio */
  /* of the largest edge length, while the error is a ratio of max_le_length. */
  y_ratio += error_ratio*copysign(fabs(y_ratio)*error_slope+error_offset, y_ratio);
  /* Maintain symmetry of coordinate in bins by mirroring the bins across the y axis. */
  int bin = log_bin(fabs(y_ratio), error_slope, error_offset);
  if(y_ratio < 0){
    bin = ~bin;
  }
  return bin;
}

/* Transform y coordinate bin into a y coordinate ratio absolute value. */
/* Returns the maximum y coordinate ratio absolute value within the bin. */
static double unbin_y(int bin,
                      unsigned int le_bin){
  /* Calculate the minimum size the largest edge ratio could have, given its bin. */
  double min_le_ratio = unbin_largest_edge(le_bin);
  /* Coordinate error is affected by both the largest edge error and the FOV error. */
  double error_constant = le_error_offset/(2-max_scale_factor);
  /* Coordinate error has a constant factor determined by the centroiding error. */
  /* This value is the worst case of maximum FOV and centroid error.  It corresponds */
  /* to the FOV taking on its minimum value, both the largest edge centroids having */
  /* maximal error in the same direction, and the coordinate centroid having maximal */
  /* error in the opposite direction of the centroiding errors of the largest edge. */
  double error_slope = error_constant/fmax(min_le_ratio-error_constant, 0);
  /* Coordinate error also has a constant factor determined by the centroiding error. */
  double error_offset = error_slope;
  /* Invert logarithmic binning function to retrieve maximum y coordinate ratio. */
  /* Returns the absolute value of the y coordinate ratio, not the actual value. */
  double max_y_ratio = log_unbin(bin>=0?bin+1:(~bin)+1, error_slope, error_offset);
  return max_y_ratio;
}

/* Evenly bin x coordinate using the y and largest edge bins. */
static int bin_x(int x,
                 unsigned int le_bin,
                 int y_bin,
                 int error_ratio){
  /* Calculate the minimum size the largest edge ratio could have, given its bin. */
  double min_le_ratio = unbin_largest_edge(le_bin);
  /* Calculate the maximum size the y coordinate ratio could have, given its bin. */
  double max_y_ratio = unbin_y(y_bin, le_bin);
  /* Coordinate error is affected by both the largest edge error and the FOV error. */
  double error_constant = le_error_offset/(2-max_scale_factor);
  /* Coordinate error has a component proportional to the coordinate value. */
  /* This value is the worst case of maximum FOV and centroid error.  It corresponds */
  /* to the FOV taking on its minimum value, both the largest edge centroids having */
  /* maximal error in the same direction, and the coordinate centroid having maximal */
  /* error in the opposite direction of the centroiding errors of the largest edge. */
  double error_slope = error_constant/fmax(min_le_ratio-error_constant, 0);
  /* Coordinate error also has a constant factor determined by the centroiding error. */
  double error_offset = error_slope*(1+2*sqrt((1.0/4)+max_y_ratio*max_y_ratio))/2;
  /* Convert coordinate back to double between -1 and 1 from implicitly divided double. */
  /* Resulting value is a ratio of image_scale. */
  double x_ratio = x/((1<<14)-1.0);
  /* Adjust error in coordinate based on error_ratio between -1 and 1. */
  /* An error_ratio of -1 will give the lower bin, +1 will give the upper bin. */
  /* The error is inversely scaled by min_le_ratio, as the coordinate is a ratio */
  /* of the largest edge length, while the error is a ratio of max_le_length. */
  x_ratio += error_ratio*copysign(fabs(x_ratio)*error_slope+error_offset, x_ratio);
  /* Maintain symmetry of coordinate in bins by mirroring the bins across the y axis. */
  int bin = log_bin(fabs(x_ratio), error_slope, error_offset);
  if(x_ratio < 0){
    bin = ~bin;
  }
  return bin;
}

/* Constant chosen for optimal avalanching.  It is the closest prime to 2^64 divided */
/* by the golden ratio.  This is a modified Knuth multiplicative hash, as it uses the */
/* closest prime instead.  A prime is chosen to prevent an unfortunate catalog size */
/* (i.e. a factor of 2^64 / phi) from effectively canceling out this step.  Note again */
/* that the result is automatically taken modulo 2^64 as it is stored as a uint64_t. */
static uint64_t hash_int(uint64_t old_hash, uint64_t key){
  key = key*11400714819323198549ULL;
  /* XOR with higher order bits to properly mix in the old hash. */
  return old_hash^(old_hash >> 13)^key;
}

/* Hash function that takes a Pattern as input and produces a deterministically */
/* random catalog position as output based on its bins.  Note that this */
/* hash function was chosen for simplicity/speed.  Other functions may have */
/* somewhat better (e.g. more even) distributions of hashed values. */
static uint64_t hash_pattern(Pattern pattern_instance){
  /* Initialize hash value to the largest edge bin. */
  unsigned int le_bin = bin_largest_edge(pattern_instance.largest_edge, 0);
  uint64_t hash = hash_int(0, le_bin);
  /* Each feature is hashed one at a time through modular multiplication. */
  for(int i=0;i<num_stars_in_pattern-2;i++){
    /* The bins are used to update the hash by multiplying it.  The hash is stored */
    /* with an unsigned 64 bit value, so the result is automatically taken modulo 2^64. */
    int y_bin = bin_y(pattern_instance.features[i].y, le_bin, 0);
    hash = hash_int(hash, y_bin+(1<<31));
    hash = hash_int(hash, bin_x(pattern_instance.features[i].x, le_bin, y_bin, 0)+(1<<31));
  }
  /* Hash value is taken modulo the catalog size to convert it to a location in the catalog. */
  /* Due to hash collisions, the pattern may not be stored in this exact location, but nearby. */
  return hash%catalog_size_in_patterns;
}

/* Verifies bin pairs are the same for the new and */
/* stored Patterns, in case a collision occurred. */
static int hash_same(Pattern new_pattern,
                     Pattern cat_pattern){
  /* Check whether the Patterns share the same largest edge and fixed star magnitude bins. */
  unsigned int new_le_bin = bin_largest_edge(new_pattern.largest_edge, 0);
  unsigned int cat_le_bin = bin_largest_edge(cat_pattern.largest_edge,
                                             2*cat_pattern.le_bin_offset-1);
  if(new_le_bin != cat_le_bin){
    /* Largest edge bins don't match, so it must be a collision. */
    return 0;
  }
  /* Iterate over both Patterns' Features, checking whether their bins match. */
  for(int i=0;i<num_stars_in_pattern-2;i++){
    Feature new_feature = new_pattern.features[i];
    Feature cat_feature = cat_pattern.features[i];
    int new_y_bin = bin_y(new_feature.y, new_le_bin, 0);
    int cat_y_bin = bin_y(cat_feature.y, cat_le_bin, 2*cat_feature.y_bin_offset-1);
    int new_x_bin = bin_x(new_feature.x, new_le_bin, new_y_bin, 0);
    int cat_x_bin = bin_x(cat_feature.x, cat_le_bin, cat_y_bin, 2*cat_feature.x_bin_offset-1);
    if((new_y_bin != cat_y_bin) || (new_x_bin != cat_x_bin)){
      /* Coordinate bins don't match, so it must be a collision. */
      return 0;
    }
  }
  /* All bins matched, so it must not be a match and not a collision. */
  return 1;
}

/* Determines whether two Patterns match by checking if their corresponding */
/* x,y coordinate pairs, magnitudes and largest edge lengths match. */
static int is_match(Pattern new_pattern,
             Pattern cat_pattern){
  /* Verify the image Patterns' largest edge length is within range of the catalog Pattern. */
  double new_le_ratio = new_pattern.largest_edge/((1<<16)-1.0);
  double cat_le_ratio = cat_pattern.largest_edge/((1<<16)-1.0);
  double max_le_error = cat_le_ratio*le_error_slope+le_error_offset;
  if(fabs(new_le_ratio-cat_le_ratio) > max_le_error){
    return 0;
  }
  /* Coordinate error is affected by both the largest edge error and the FOV error. */
  double coord_error_constant = le_error_offset/(2-max_scale_factor);
  /* Coordinate error has a constant factor determined by the centroiding error. */
  double coord_error_slope = coord_error_constant/fmax(new_le_ratio-coord_error_constant, 0);
  /* Coordinate error also has a component proportional to the coordinate value. */
  double coord_error_offset_y = coord_error_slope;
  /* Iterate over both Patterns' Features, checking that their x,y coordinate pairs */
  /* and their magnitudes are within range. */
  for(int i=0;i<num_stars_in_pattern-2;i++){
    /* Retrieve both Features' y coordinates and verify they're within range. */
    double new_y = new_pattern.features[i].y/((1<<14)-1.0);
    double cat_y = cat_pattern.features[i].y/((1<<14)-1.0);
    double max_y_error = fabs(cat_y)*coord_error_slope+coord_error_offset_y;
    if(fabs(new_y-cat_y) > max_y_error){
      /* If a coordinate is out of range, the Patterns do not match. */
      return 0;
    }
  }
  unsigned int cat_le_bin = bin_largest_edge(cat_pattern.largest_edge,
                                             2*cat_pattern.le_bin_offset-1);
  for(int i=0;i<num_stars_in_pattern-2;i++){
    int cat_y_bin = bin_y(cat_pattern.features[i].y, cat_le_bin, 2*cat_pattern.features[i].y_bin_offset-1);
    double max_y_ratio = unbin_y(cat_y_bin, cat_le_bin);
    double coord_error_offset_x = coord_error_slope*(1+2*sqrt((1.0/4)+max_y_ratio*max_y_ratio))/2;
    /* Retrieve both Features' x coordinates and verify they're within range. */
    double new_x = new_pattern.features[i].x/((1<<14)-1.0);
    double cat_x = cat_pattern.features[i].x/((1<<14)-1.0);
    double max_x_error = fabs(cat_x)*coord_error_slope+coord_error_offset_x;
    if(fabs(new_x-cat_x) > max_x_error){
      /* If a coordinate is out of range, the Patterns do not match. */
      return 0;
    }
  }
  /* Every pair of values matched, so the Patterns match. */
  return 1;
}

/* Quadratically probes the Pattern cache with bounds checking. */
static int increment_offset(FILE *pattern_catalog,
                            Pattern catalog_pattern_cache[pattern_cache_size],
                            uint64_t *offset,
                            int *cache_offset,
                            int *probe_step){
  /* If probe goes out of probe bounds, return failure. */
  if(((*probe_step)*(*probe_step+1))/2 > max_probe_depth){
    return 0;
  }
  /* Update cache offset to next probe offset. */
  *cache_offset += *probe_step;
  /* If cache_offset goes out of cache bounds, move cache to next probe offset. */
  if (*cache_offset >= pattern_cache_size) {
    /* Update offset of cache within catalog to next probe offset. */
    *offset += *cache_offset;
    /* Reset offset within cache to the start of the cache. */
    *cache_offset = 0;
    /* Cache section of catalog beginning at the next probe offset. */
    /* Note that this may read beyond the end of the file, but */
    /* will not use any values from beyond the end of the file. */
    /* If reading causes a crash, it can be fixed by taking a max */
    /* of the number of read Patterns with the sum of the catalog */
    /* size and the maximum probe depth minus the offset. */
    _fseeki64(pattern_catalog, (*offset)*sizeof(Pattern), SEEK_SET);
    fread(catalog_pattern_cache, sizeof(Pattern), pattern_cache_size, pattern_catalog);
  }
  *probe_step += 1;
  /* If probe stayed within bounds, return success. */
  return 1;
}

/* Retrieves catalog pattern matching image Pattern.  Returns 1 if a unique */
/* match is found.  Returns 0 if multiple matches or no matches are found. */
static int get_matching_pattern(Pattern image_pattern,
                                Pattern *catalog_pattern,
                                FILE *pattern_catalog){
	/* Cache of Pattern instances from catalog. */
  static Pattern catalog_pattern_cache[pattern_cache_size];
  /* Initialize offset of the Pattern's probing sequence in the pattern cache. */
  int cache_offset = 0;
  /* Spacing between Patterns in the same probing sequence.  Grows linearly, */
  /* resulting in quadratic probing. (i.e. 0,1,3,6,10,15...) */
  int probe_step = 1;
  /* Boolean representing whether or not a catalog match has been found yet. */
  /* A match is a catalog Pattern within max_coord_error of the given image Pattern. */
  int found_match = 0;
  /* Initialize offset of the beginning of the Pattern cache in the pattern catalog. */
  uint64_t offset = hash_pattern(image_pattern);
  /* Initialize cache of catalog Patterns. */
  _fseeki64(pattern_catalog, offset*sizeof(Pattern), SEEK_SET);
  fread(catalog_pattern_cache, sizeof(Pattern), pattern_cache_size, pattern_catalog);
  /* Iterate over catalog locations in the image Pattern's probing sequence until */
  /* a catalog location without a Pattern is found, which means no matches exist. */
  /* If the probing sequence contains Patterns with the same sub-bins, iterate up */
  /* to the last one, returning success if a single matching Pattern is found, */
  /* and returning failure if no matching Patterns are found before the last one. */
  /* If two or more matching catalog Patterns are found, exit early and */
  /* return failure, as a unique identification cannot be made. */
  while(catalog_pattern_cache[cache_offset].has_pattern){
    /* Only examine catalog Patterns with the same sub-bins as the image Pattern, */
    /* as all matches will also have the same sub-bins as the image Pattern. */
    if(hash_same(image_pattern, catalog_pattern_cache[cache_offset])){
      /* Check whether the image and catalog Patterns are a match by checking if */
      /* their corresponding Features' coordinates are all within max_coord_error. */
      if(is_match(image_pattern, catalog_pattern_cache[cache_offset])){
        /* If a match has already been found previously, this must be the second */
        /* matching catalog Pattern found.  Return failure, as a unique */
        /* identification cannot be made without checking other stars. */
        /* Note that it may be worth checking all possible matches if */
        /* verifying a match is less costly than another catalog access. */
        if(found_match){
          /* Multiple matching Patterns were found.  Return failure. */
          return 0;
        }
        /* This must be the first matching catalog Pattern found.  Store it as the */
        /* output Pattern so it will be output once its uniqueness has been verified. */
        *catalog_pattern = catalog_pattern_cache[cache_offset];
        /* Set the flag indicating a matching Pattern has already been found. */
        found_match = 1;
      }
      /* If this is the last Pattern with the same sub-bins, there is no need to */
      /* keep searching until an empty catalog location is found, as any Patterns */
      /* beyond this point cannot be matches.  Exit early to save time. */
      if(catalog_pattern_cache[cache_offset].is_last){
        break;
      }
    }
    /* Advance to the next catalog location given by quadratic probing. */
    /* If cache_offset indexes beyond pattern_cache_size, return failure. */
    if(!increment_offset(pattern_catalog,
                         catalog_pattern_cache,
                         &offset,
                         &cache_offset,
                         &probe_step)){
      return 0;
    }
  }
  /* Exactly one matching Pattern was found.  Return success. */
  if(found_match){
    return 1;
  }
  /* No matching Patterns were found.  Return failure. */
  return 0;
}

/* Identifies stars given normalized camera vectors as input.  Uses (i,j,k) vectors */
/* where the image center is at (1,0,0).  One normalized vector consists of */
/* (forward, left, up).  Each vector points at a star on the celestial sphere in */
/* the camera reference frame.  Returns 1 on identification success, placing */
/* matching pairs of camera vector indices and catalog star ids into the matches */
/* array in the same order they appear in their Patterns' Features. */
static int identify_stars(double image_stars[max_stars][3],
                    int image_star_ids[num_stars_in_pattern],
                    FILE *pattern_catalog,
                    int matches[num_stars_in_pattern][2]){
  /* Dummy variable used for iteration. */
  int i,j;
  /* Image Pattern created from the input camera vectors. */
	Pattern new_pattern;
	/* Catalog Pattern which uniquely matches the image Pattern. */
	Pattern catalog_pattern;
  /* Iterate over all pairs of stars in the pattern to */
  /* find the length of the largest edge. */
  double largest_edge_length = 0.0;
  for(i=0;i<num_stars_in_pattern;i++){
    for(j=i+1;j<num_stars_in_pattern;j++){
      /* Track the largest edge so far by comparing once with each edge. */
      double new_edge_length = dist(image_stars[image_star_ids[i]],
                                   image_stars[image_star_ids[j]]);
      /* Check if the new edge is larger than the current largest edge. */
      if(new_edge_length > largest_edge_length){
        /* Set the new edge as the current largest edge length. */
        largest_edge_length = new_edge_length;
        /* Store the star ids composing the largest edge in the image Pattern. */
        /* Their order may be swapped later to resolve the 180 degree ambiguity. */
        new_pattern.fixed_star_id_1 = image_star_ids[i];
        new_pattern.fixed_star_id_2 = image_star_ids[j];
      }
    }
  }
  /* Set Pattern's largest_edge_length and largest_edge_subbin.  Both are */
  /* implicitly encoded as unsigned integers with a range from 0 up to */
  /* the sine of the maximum catalog FOV to keep the catalog compact. */
  new_pattern.largest_edge = (largest_edge_length/max_le_length)*((1<<16)-1);
  /* Calculate vector along x axis of Pattern's coordinate system. */
  /* The vector points from the first fixed Pattern star to the second. */
  double x_axis_vector[3];
  diff(image_stars[new_pattern.fixed_star_id_2],
       image_stars[new_pattern.fixed_star_id_1],
       x_axis_vector);
  /* Calculate vector along y axis of Pattern's coordinate system. */
  double y_axis_vector[3];
  cross_prod(image_stars[new_pattern.fixed_star_id_2],
             image_stars[new_pattern.fixed_star_id_1],
             y_axis_vector);
  /* Normalize axis vectors to unit length by dividing by their magnitudes. */
  normalize(x_axis_vector);
  normalize(y_axis_vector);
  /* Use the remaining stars to initialize the Pattern's Features. */
  int feature_index = 0;
  for(i=0;i<num_stars_in_pattern;i++){
    /* Skip the fixed star ids, as they don't have their own Features. */
    if(image_star_ids[i] != new_pattern.fixed_star_id_1 &&
       image_star_ids[i] != new_pattern.fixed_star_id_2){
      /* Set the Feature's star id to match its corresponding star. */
      new_pattern.features[feature_index].star_id = image_star_ids[i];
      /* Calculate the normalized x and y coordinates using vector projection. */
      double x = dot_prod(x_axis_vector, image_stars[image_star_ids[i]])/largest_edge_length;
      double y = dot_prod(y_axis_vector, image_stars[image_star_ids[i]])/largest_edge_length;
      /* Set Feature's coordinates by converting to implicitly divided integers. */
      new_pattern.features[feature_index].x = x*((1<<14)-1);
      new_pattern.features[feature_index].y = y*((1<<14)-1);
      /* Disallow 0, as rotational ambiguity correction would fail. */
      if(new_pattern.features[feature_index].x == 0){
        new_pattern.features[feature_index].x = 1;
      }
      if(new_pattern.features[feature_index].y == 0){
        new_pattern.features[feature_index].y = 1;
      }
      feature_index++;
    }
  }
  /* Variable encoding which 180 degree rotation will be inserted into the catalog. */
  /* A negative value means the current rotation will be inserted. */
  /* A positive value means the opposite rotation will be inserted. */
  /* A value of zero means both rotations will be inserted into the catalog. */
  int pattern_rotation;
  /* Compute largest edge bin for use in sorting Features based on x and y bins. */
  unsigned int le_bin = bin_largest_edge(new_pattern.largest_edge, 0);
  /* Helper function for sorting Features.  Sorts by x bin, then by y bin. */
  /* Returns a positive number if the first Feature has larger bin values, */
  /* returns a negative number if the second Feature has larger bin values, */
  /* and raises an error if both Features have the same bin values. */
  int compare_bins(const void *p, const void *q) {
    /* Compare the Features' x bins first, then their y bins. */
    int p_y_bin = bin_y(((Feature*)p)->y, le_bin, 0);
    int q_y_bin = bin_y(((Feature*)q)->y, le_bin, 0);
    int p_x_bin = bin_x(((Feature*)p)->x, le_bin, p_y_bin, 0);
    int q_x_bin = bin_x(((Feature*)q)->x, le_bin, q_y_bin, 0);
    /* If the x bins have different values, the y bins don't need to be computed. */
    if(p_x_bin != q_x_bin){
      return p_x_bin-q_x_bin;
    }
    return p_y_bin-q_y_bin;
  }
  /* Sort Pattern's Features based on coordinate bins to give a unique ordering. */
  qsort(new_pattern.features, num_stars_in_pattern-2, sizeof(Feature), compare_bins);
  /* Create a copy of the first Feature of the Pattern. */
  Feature first_feature = new_pattern.features[0];
  /* Rotate the copy by 180 degrees by taking complements of its coordinates. */
  first_feature.x = -first_feature.x;
  first_feature.y = -first_feature.y;
  /* Compare with the last Feature's bins to determine which has the largest */
  /* x bin (with y bin as a tie breaker).  This will determine the */
  /* 180 degree rotation of the Pattern's coordinate system.  The orientation */
  /* which gives the larger Feature a positive x bin value is chosen. */
  /* Put another way, the Feature furthest from the y-axis is placed on the right. */
  /* In the case that the first and last Features' bins are ambiguous after */
  /* rotating the first Feature by 180 degrees, both orientations are inserted. */
  pattern_rotation = compare_bins((void*)&first_feature,
                                  (void*)&(new_pattern.features[num_stars_in_pattern-3]));
  /* If the current rotation is incorrect, rotate the Pattern by 180 degrees by taking */
  /* the complement of its Features' bin offsets and coordinates, reversing the order */
  /* of its Features, and swapping its fixed stars before inserting it into the catalog. */
  if(pattern_rotation >= 0){
    for(i=0;i<num_stars_in_pattern-2;i++){
      /* Take the complement of each Feature's coordinates. */
      new_pattern.features[i].x = -new_pattern.features[i].x;
      new_pattern.features[i].y = -new_pattern.features[i].y;
    }
    /* Reverse the order of the Pattern's Features by swapping across the middle. */
    for(i=0;i<(num_stars_in_pattern-2)/2;i++){
      Feature feature_swap = new_pattern.features[i];
      new_pattern.features[i] = new_pattern.features[num_stars_in_pattern-3-i];
      new_pattern.features[num_stars_in_pattern-3-i] = feature_swap;
    }
    /* Swap the order of the Pattern's fixed star ids and magnitudes. */
    unsigned int fixed_star_id_swap = new_pattern.fixed_star_id_1;
    new_pattern.fixed_star_id_1 = new_pattern.fixed_star_id_2;
    new_pattern.fixed_star_id_2 = fixed_star_id_swap;
  }

  /* Check cached section of catalog for Patterns matching image Pattern. */
	if(!get_matching_pattern(new_pattern, &catalog_pattern, pattern_catalog)){
		return 0;
	}
  /* Create matching pairs of stars by corresponding fixed_star_ids and */
  /* Feature star ids between the image Pattern and catalog Pattern.  */
  matches[0][0] = new_pattern.fixed_star_id_1;
  matches[1][0] = new_pattern.fixed_star_id_2;
  matches[0][1] = catalog_pattern.fixed_star_id_1;
  matches[1][1] = catalog_pattern.fixed_star_id_2;
  for(i=0;i<num_stars_in_pattern-2;i++){
    matches[i+2][0] = new_pattern.features[i].star_id;
    matches[i+2][1] = catalog_pattern.features[i].star_id;
  }
	return 1;
}

/* Identifies an image given normalized camera vectors as input.  Uses (i,j,k) vectors */
/* where the image center is at (1,0,0).  One normalized vector consists of */
/* (forward, left, up).  Each vector points at a star on the celestial sphere in */
/* the camera reference frame.  Writes matching star ids to file. */
static int identify_image(double image_stars[max_stars][3],
                   FILE *pattern_catalog,
                   int num_image_stars,
                   int matches[num_stars_in_pattern][2],
                   int num_stars_selected){
  /* Array of image star ids for a given Pattern.  Being static avoids needing */
  /* to pass the array between recursive calls of the function. */
  static int image_star_ids[num_stars_in_pattern];
  /* If some image Pattern stars still need to be selected, select the next star and recurse. */
  if(num_stars_selected < num_stars_in_pattern){
    for(image_star_ids[num_stars_selected] = num_stars_in_pattern-num_stars_selected-1;
        image_star_ids[num_stars_selected] < num_image_stars;
        image_star_ids[num_stars_selected]++){
      /* The next star to be selected must have an id lower than the previous selection. */
      if(identify_image(image_stars,
                       pattern_catalog,
                       image_star_ids[num_stars_selected],
                       matches,
                       num_stars_selected+1)){
        return 1;
      }
    }
  }
  /* Once all of the image Pattern stars have been selected, identify the stars. */
  else{
    if(identify_stars(image_stars, image_star_ids, pattern_catalog, matches)){
      return 1;
    }
  }
  return 0;
}

int main(int argc,char *argv[] ) {
	clock_t begin, end;
	int ms_spent;
	/* Tetra catalog file pointer */
	FILE *pattern_catalog;
  int matches[num_stars_in_pattern][2];
	/* Size of sides of image in pixels. */
	int num_pixels_x = 1024;
	int num_pixels_y = 1024;
	/* Open pattern catalog file. */
	pattern_catalog = fopen("pattern_catalog","rb");
	if (!pattern_catalog){
		printf("Unable to open pattern catalog file!\n");
		return 1;
	};
  FILE *centroid_data;
	/* open centroid memmap */
	centroid_data = fopen("centroid_data.p","rb");
	if (!centroid_data){
		printf("Unable to open centroid file!");
		return 1;
	};
	/* load image data from file for answer checking */
	uint16_t *image_data = (uint16_t *)malloc(sizeof(uint16_t)*num_images*max_stars_per_image);
	FILE *image_data_file = fopen("image_data.p", "rb");
	if (!image_data_file){
		printf("Unable to open image file!\n");
		return 1;
	};
	fread(image_data, sizeof(uint16_t), num_images * max_stars_per_image, image_data_file);	
	/* Close Image Data File */
	fclose(image_data_file);
	/* fov of image along one axis */
	for (double fov = 10.1; fov < 10.1001; fov += .05) {
    double fov_factor = (tan(fov*PI/360)*2)/num_pixels_x;
    /* centroiding error in 1/100ths of pixels */
    for (int centroid_error = 100; centroid_error < 101; centroid_error += 10) {
      begin = clock();
      int right = 0;
      int failed = 0;
      int wrong = 0;
      int wrong_order = 0;
      int too_few_stars = 0;
      srand(0);
      int count = 0;
      for (int file_index = 0; file_index < num_images; file_index++) {
        int i;
        int num_image_stars;
        float centroids[max_stars][2];
        fseek(centroid_data, sizeof(float) * 2 * max_stars_per_image * file_index, SEEK_SET);
        fread(centroids, sizeof(float), max_stars*2, centroid_data);
        for (num_image_stars = 0; num_image_stars < max_stars; num_image_stars++) {
          if (centroids[num_image_stars][0] == 0 && centroids[num_image_stars][1] == 0) {
            break;
          }
        }
        if (num_image_stars < num_stars_in_pattern) {
          too_few_stars += 1;
          continue;
        }
        /* (i,j,k,mag) vectors with magnitude where the image center is at (1,0,0). */
        /* One normalized vector consists of (forward, left, up). */
        /* Each vector points at a star on the celestial sphere in */
        /* the camera reference frame. */
        double image_stars[max_stars][3] = {{0}};
        /* convert centroids to camera vectors */
        double x;
        double y;
        for (i = 0; i < num_image_stars; i++) {
          float rand1 = (float)((rand() % 10000)+1) / 10000;
          float rand2 = (float)((rand() % 10000)+1) / 10000;
          count += 1;
          float x_err = 2 * (rand1 - .5);
          float y_err = 2 * (rand2 - .5);
          while (sqrt(x_err*x_err + y_err*y_err) > 1.0) {
            rand1 = (float)((rand() % 10000)+1) / 10000;
            rand2 = (float)((rand() % 10000)+1) / 10000;
            count += 1;
            x_err = 2 * (rand1 - .5);
            y_err = 2 * (rand2 - .5);
          }
          x = (centroids[i][0] + (x_err * centroid_error / 100)) * fov_factor;
          y = (centroids[i][1] + (y_err * centroid_error / 100)) * fov_factor;
          image_stars[i][0] = 1 / sqrt(1 + x*x + y*y);
          image_stars[i][1] = -image_stars[i][0] * x;
          image_stars[i][2] = image_stars[i][0] * y;
        }
        if(identify_image(image_stars, pattern_catalog, num_image_stars, matches, 0)){
          /* Check matches. */
          int num_verified_stars = 0;
          for (int v = 0; v < max_stars_per_image; v++) {
            uint16_t image_data_star = image_data[file_index*max_stars_per_image + v];
            for(int match_id = 0; match_id < num_stars_in_pattern; match_id++){
              if(matches[match_id][1] == image_data_star){
                num_verified_stars += 1;
                break;
              }
            }
          }
          if (num_verified_stars < num_stars_in_pattern) {
            wrong += 1;
          }
          else {
            int mismatch = 0;
            for(int match_id = 0; match_id < num_stars_in_pattern; match_id++){
              if(matches[match_id][1] != image_data[file_index*max_stars_per_image + matches[match_id][0]]){
                mismatch = 1;
                break;
              }
            }
            if(mismatch){
              wrong_order += 1;
            }
            else{
              right += 1;
            }
          }							
        }
        else{
          failed += 1;
        }
      }
      end = clock();
      ms_spent = (int) 1000 * (end - begin) / CLOCKS_PER_SEC;
      /* Print to Screen */
      printf("fov: %.2f\n", roundf(fov * 100) / 100);
      printf("centroid_error: %.2f\n", ((float) centroid_error / 100));
      printf("number right: %d\n", right);
      printf("number failed: %d\n", failed);
      printf("num wrngordr: %d\n", wrong_order);
      printf("number wrong: %d\n", wrong);
      printf("num 2 few stars: %d\n", too_few_stars);
      printf("ms taken: %d\n", ms_spent);
      /* Write to Output File */
      FILE *output_file = fopen("output.txt", "a");
      fprintf(output_file, "fov: %.2f\n", roundf(fov * 100) / 100);
      fprintf(output_file, "centroid_error: %.2f\n", ((float) centroid_error / 100));
      fprintf(output_file, "number right: %d\n", right);
      fprintf(output_file, "number failed: %d\n", failed);
      fprintf(output_file, "num wrngordr: %d\n", wrong_order);
      fprintf(output_file, "number wrong: %d\n", wrong);
      fprintf(output_file, "ms taken: %d\n", ms_spent);
      fclose(output_file);
    }
  }
  /* Close Dictionary Memmap */
  fclose(pattern_catalog);	
  /* Close centroid Memmap */
  fclose(centroid_data);
};