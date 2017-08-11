#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

/* User defined constants. */
/* Number of stars in star pattern. */
/* Must be at least 3.  Recommended number is 4. */
#define num_stars_in_pattern 4
/* Minimum star brightness (in magnitude) for inclusion in catalog. */
/* Note that lower magnitude values are brighter. */
#define min_magnitude 6.2
/* Maximum Field of View for catalog in radians. */
/* Also the maximum angle between any two stars in a tetrahedron. */
/* Typically equal to the angle subtended by the imager's diagonal. */
/* Must be less than pi, but should be much less to prevent invalidation */
/* of small angle approximation and to minimize non-linearity of FOV error. */
/* .247 radians is about 14.1 degrees or the diagonal of a 10 degree FOV */
#define max_fov .247
/* Maximum star coordinate centroiding error as fraction of maximum FOV. */
/* .001 is .1% of the max FOV or 1.414 pixels in a 1000x1000 image. */
// 1 / (1024*sqrt(2)) < .00069054
#define max_centroid_error .00069054
/* Maximum error in imager FOV estimate as fraction of true imager FOV. */
/* max_fov*(1+max_fov_error) must be less than pi, but should be much less. */
/* .01 for a 10 degree FOV imager covers estimates from 9.9 to 10.1 degrees. */
#define max_fov_error 0.01
/* Minimum angle between stars in catalog in radians. */
/* Optionally used to remove double stars, set to 0 otherwise. */
/* .001 is about 5.7 pixels distance with 1000 pixels over a 10 degree FOV. */
#define min_star_separation .004
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
/* Ratio specifying how much of the catalog is occupied by Patterns rather than */
/* empty space.  Values below .5 create unnecessarily large catalog sizes, while */
/* values above .7 cause longer lookup times due to increased collision frequency. */
#define catalog_density .6
/* Size of catalog cache in number of patterns.  Generating the catalog in cached */
/* pieces reduces disk I/O and greatly speeds up catalog generation. */
#define max_cat_cache_size 100000000
/* Maximum distance matching pattern can be from original offset in catalog. */
/* Measured in number of patterns.  Determines overlap between cached catalog */
/* pieces during catalog generation and the number of extra catalog positions */
/* needed at the end of the catalog due to the catalog hash table being non-cyclic. */
#define max_probe_depth 50000
/* Number of entries in the Hipparchos catalog. */
#define STARN 9110
/* Mathematical constant pi's approximate value. */
#define PI 3.1415926
/* The current calendar year. */
#define current_year 2017

/* The following values are not user defined constants and should therefore not be changed. */
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

/* BSC_Entry struct */
/* Data format of what's stored in a single star entry of the BSC5. */
/* The BSC5 is Yale's Bright Star Catalog 5, containing 9110 stars. */
typedef struct BSC_Entry
{
  /* Catalog number of star. */
  float XNO;
  /* B1950 Right Ascension (radians) */
  double SRA0;
  /* B1950 Declination (radians) */
  double SDEC0;
  /* Spectral type (2 characters) */
  short IS;
  /* V Magnitude * 100 */
  short MAG;
  /* R.A. proper motion (radians per year) */
  float XRPM;
  /* Dec. proper motion (radians per year) */
  float XDPM;
} BSC_Entry;

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
static uint64_t hash_pattern(Pattern pattern_instance,
                             uint64_t catalog_size_in_patterns){
  /* Initialize hash value to the largest edge bin. */
  unsigned int le_bin = bin_largest_edge(pattern_instance.largest_edge,
                                         2*pattern_instance.le_bin_offset-1);
  uint64_t hash = hash_int(0, le_bin);
  /* Each feature is hashed one at a time through modular multiplication. */
  for(int i=0;i<num_stars_in_pattern-2;i++){
    /* The bins are used to update the hash by multiplying it.  The hash is stored */
    /* with an unsigned 64 bit value, so the result is automatically taken modulo 2^64. */
    int y_bin = bin_y(pattern_instance.features[i].y, le_bin,
                                    2*pattern_instance.features[i].y_bin_offset-1);
    hash = hash_int(hash, y_bin+(1<<31));
    hash = hash_int(hash, bin_x(pattern_instance.features[i].x, le_bin, y_bin,
                                2*pattern_instance.features[i].x_bin_offset-1)+(1<<31));
  }
  /* Hash value is taken modulo the catalog size to convert it to a location in the catalog. */
  /* Due to hash collisions, the pattern may not be stored in this exact location, but nearby. */
  return hash%catalog_size_in_patterns;
}

/* Verifies bin pairs are the same for the new and */
/* stored Patterns, in case a collision occurred. */
static int hash_same(Pattern new_pattern,
                     Pattern cat_pattern){
  /* Check whether the Patterns share the same largest edge bins. */
  unsigned int new_le_bin = bin_largest_edge(new_pattern.largest_edge,
                                             2*new_pattern.le_bin_offset-1);
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
    int new_y_bin = bin_y(new_feature.y, new_le_bin, 2*new_feature.y_bin_offset-1);
    int cat_y_bin = bin_y(cat_feature.y, cat_le_bin, 2*cat_feature.y_bin_offset-1);
    int new_x_bin = bin_x(new_feature.x, new_le_bin, new_y_bin, 2*new_feature.x_bin_offset-1);
    int cat_x_bin = bin_x(cat_feature.x, cat_le_bin, cat_y_bin, 2*cat_feature.x_bin_offset-1);
    if((new_y_bin != cat_y_bin) || (new_x_bin != cat_x_bin)){
      /* Coordinate bins don't match, so it must be a collision. */
      return 0;
    }
  }
  /* All bins matched, so it must not be a match and not a collision. */
  return 1;
}

/* Quadratically probes the Pattern cache with bounds checking. */
static int increment_offset(uint64_t *cache_offset,
                            int *probe_step){
  *cache_offset += *probe_step;
  static int deepest_probe = 0;
  /* Keep track of catalog's deepest probe depth. */
  if(((*probe_step)*(*probe_step+1))/2 > deepest_probe){
    deepest_probe = ((*probe_step)*(*probe_step+1))/2;
  }
  /* If cache_offset goes out of bounds, terminate program with failure code. */
  /* The problem can be fixed by increasing the max_probe_depth, lowering the */
  /* catalog_density, or otherwise lowering the collision frequency. */
  if(((*probe_step)*(*probe_step+1))/2 > max_probe_depth){
    printf("\nMaximum probe depth exceeded: increase max_probe_depth.\n");
    exit(EXIT_FAILURE);
  }
  *probe_step += 1;
  return deepest_probe;
}

/* Inserts a Pattern instance into the catalog cache based on its bins. */
/* Also sorts the Pattern's Features before inserting it into the catalog. */
/* If the pattern would be placed outside of cache range, returns without doing anything. */
static void insert_pattern(int just_count,
                           uint64_t *num_patterns,
                           Pattern pattern_catalog_cache[],
                           uint64_t *catalog_size_in_patterns,
                           int *cache_index,
                           Pattern new_pattern){
  /* If set to "just count," increment the counter of valid Patterns and return. */
  if(just_count){
    (*num_patterns)++;
    return;
  }
  /* Offset of the beginning of the Pattern's probing sequence in the pattern catalog. */
  uint64_t offset = hash_pattern(new_pattern, *catalog_size_in_patterns);
  /* If the Pattern's initial offset isn't in the cached section */
  /* of memory, return without doing anything. */
  if(offset/max_cat_cache_size != *cache_index){
    return;
  }
  /* Spacing between Patterns in the same probing sequence.  Grows linearly, */
  /* resulting in quadratic probing. (i.e. 0,1,3,6,10,15...) */
  int probe_step = 1;
  /* The new pattern will always be the last match in the probing sequence so far. */
  new_pattern.is_last = 1;
  /* Get offset in catalog cache. */
  uint64_t cache_offset = offset%max_cat_cache_size;
  /* If there is already a Pattern instance at this catalog position and either its */
  /* bin pairs don't match with the new Pattern or it isn't the last matching catalog */
  /* Pattern, continue on to the next catalog position given by quadratic probing. */
  while(pattern_catalog_cache[cache_offset].has_pattern &&
        (!hash_same(new_pattern, pattern_catalog_cache[cache_offset]) ||
         !pattern_catalog_cache[cache_offset].is_last)){
    /* Advance to the next catalog location given by quadratic probing. */
    increment_offset(&cache_offset, &probe_step);
  }
  /* If there is a Pattern instance already occupying this catalog position, */
  /* it must be the last matching Pattern, so set it to not be the last one */
  /* and place the new Pattern in the next unoccupied catalog position given by */
  /* quadratic probing as the new last matching Pattern.  If there is no Pattern */
  /* instance already occupying this catalog position, just add the new Pattern */
  /* as there must be no matching Patterns already in the catalog. */
  if(pattern_catalog_cache[cache_offset].has_pattern){
    pattern_catalog_cache[cache_offset].is_last = 0;
    while(pattern_catalog_cache[cache_offset].has_pattern){
      /* Advance to the next catalog location given by quadratic probing. */
      increment_offset(&cache_offset, &probe_step);
    }
  }
  pattern_catalog_cache[cache_offset] = new_pattern;
}

/* Duplicates the Pattern across the catalog given ambiguities in its x,y coordinate bins. */
/* Called by disambiguate_rotation, and calls insert_pattern. */
static void disambiguate_feature_order(int just_count,
                                       uint64_t *num_patterns,
                                       Pattern pattern_catalog_cache[],
                                       uint64_t *catalog_size_in_patterns,
                                       int *cache_index,
                                       Pattern new_pattern,
                                       int feature_index){
  /* If all of the ambiguously ordered Features have been disambiguated, insert the Pattern. */
  if(feature_index >= num_stars_in_pattern-3){
    insert_pattern(just_count,
                   num_patterns,
                   pattern_catalog_cache,
                   catalog_size_in_patterns,
                   cache_index,
                   new_pattern);
  }
  /* Otherwise, search for and disambiguate ambiguously ordered Features. */
  else{
    /* Helper function for permuting Features.  Produces all permutations */
    /* of Features with index >= from_index and <= to_index. */
    void permute_features(int from_index, int to_index){
      /* If from_index is equal to to_index, recurse down starting from the next index. */
      if(from_index == to_index){
        disambiguate_feature_order(just_count,
                                   num_patterns,
                                   pattern_catalog_cache,
                                   catalog_size_in_patterns,
                                   cache_index,
                                   new_pattern,
                                   to_index + 1);
      }
      /* Otherwise, there must be ambiguous Features, which need to be permuted. */
      else{
        /* Helper function for swapping two Features. */
        void swap_features(int index_1, int index_2){
          Feature feature_swap = new_pattern.features[index_1];
          new_pattern.features[index_1] = new_pattern.features[index_2];
          new_pattern.features[index_2] = feature_swap;
        }
        /* Permute ambiguous Features using recursion and swapping. */
        for(int j=from_index;j<=to_index;j++){
          swap_features(from_index, j);
          permute_features(from_index+1, to_index);
          swap_features(from_index, j);
        }
      }
    }
    /* Compute largest edge bin for use in calculating Features' x, y coordinate bins. */
    unsigned int le_bin = bin_largest_edge(new_pattern.largest_edge,
                                           2*new_pattern.le_bin_offset-1);
    /* Find the index of the last Feature with bins that match the feature_index Feature. */
    int bin_y_1 = bin_y(new_pattern.features[feature_index].y, le_bin,
                            2*(new_pattern.features[feature_index].y_bin_offset)-1);
    int bin_x_1 = bin_x(new_pattern.features[feature_index].x, le_bin, bin_y_1,
                            2*(new_pattern.features[feature_index].x_bin_offset)-1);
    int i;
    for(i=feature_index+1;i<num_stars_in_pattern-2;i++){
      int bin_y_2 = bin_y(new_pattern.features[i].y, le_bin,
                              2*(new_pattern.features[i].y_bin_offset)-1);
      int bin_x_2 = bin_x(new_pattern.features[i].x, le_bin, bin_y_2,
                              2*(new_pattern.features[i].x_bin_offset)-1);
      if(bin_x_1 != bin_x_2 ||
         bin_y_1 != bin_y_2){
        break;
      }
    }
    /* Insert all permutations of the ambiguous Features into the catalog and recurse. */
    permute_features(feature_index, i-1);
  }
}

/* Duplicates the Pattern across the catalog given ambiguities in its 180 degree */
/* rotation.  Called by disambiguate_bins, and calls disambiguate_feature_order. */
static void disambiguate_rotation(int just_count,
                                  uint64_t *num_patterns,
                                  Pattern pattern_catalog_cache[],
                                  uint64_t *catalog_size_in_patterns,
                                  int *cache_index,
                                  Pattern new_pattern){
  /* Dummy variable used for iteration. */
  int i;
  /* Variable encoding which 180 degree rotation will be inserted into the catalog. */
  /* A negative value means the current rotation will be inserted. */
  /* A positive value means the opposite rotation will be inserted. */
  /* A value of zero means both rotations will be inserted into the catalog. */
  int pattern_rotation;
  /* Compute largest edge bin for use in sorting Features based on x and y bins. */
  unsigned int le_bin = bin_largest_edge(new_pattern.largest_edge,
                                         2*new_pattern.le_bin_offset-1);
  /* Helper function for sorting Features.  Sorts by x bin, then by y bin. */
  /* Returns a positive number if the first Feature has larger bin values, */
  /* returns a negative number if the second Feature has larger bin values, */
  /* and raises an error if both Features have the same bin values. */
  int compare_bins(const void *p, const void *q) {
    /* Compare the Features' x bins first, then their y bins. */
    int p_y_bin = bin_y(((Feature*)p)->y, le_bin, 2*(((Feature*)p)->y_bin_offset)-1);
    int q_y_bin = bin_y(((Feature*)q)->y, le_bin, 2*(((Feature*)q)->y_bin_offset)-1);
    int p_x_bin = bin_x(((Feature*)p)->x, le_bin, p_y_bin, 2*(((Feature*)p)->x_bin_offset)-1);
    int q_x_bin = bin_x(((Feature*)q)->x, le_bin, q_y_bin, 2*(((Feature*)q)->x_bin_offset)-1);
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
  /* If the current rotation is correct, insert the Pattern as-is. */
  if(pattern_rotation <= 0){
    disambiguate_feature_order(just_count,
                               num_patterns,
                               pattern_catalog_cache,
                               catalog_size_in_patterns,
                               cache_index,
                               new_pattern,
                               0);
  }
  /* If the current rotation is incorrect, rotate the Pattern by 180 degrees by taking */
  /* the complement of its Features' bin offsets and coordinates, reversing the order */
  /* of its Features, and swapping its fixed stars before inserting it into the catalog. */
  if(pattern_rotation >= 0){
    for(i=0;i<num_stars_in_pattern-2;i++){
      /* Take the complement of each Feature's bin offsets and coordinates. */
      new_pattern.features[i].x = -new_pattern.features[i].x;
      new_pattern.features[i].y = -new_pattern.features[i].y;
    }
    /* Reverse the order of the Pattern's Features by swapping across the middle. */
    for(i=0;i<(num_stars_in_pattern-2)/2;i++){
      Feature feature_swap = new_pattern.features[i];
      new_pattern.features[i] = new_pattern.features[num_stars_in_pattern-3-i];
      new_pattern.features[num_stars_in_pattern-3-i] = feature_swap;
    }
    /* Swap the order of the Pattern's fixed star ids. */
    unsigned int fixed_star_id_swap = new_pattern.fixed_star_id_1;
    new_pattern.fixed_star_id_1 = new_pattern.fixed_star_id_2;
    new_pattern.fixed_star_id_2 = fixed_star_id_swap;
    /* Insert the 180 degree rotated Pattern into the catalog. */
    disambiguate_feature_order(just_count,
                               num_patterns,
                               pattern_catalog_cache,
                               catalog_size_in_patterns,
                               cache_index,
                               new_pattern,
                               0);
  }
}

/* Calculate the difference between the maximum and minimum bins. */
/* Also verify minimum and maximum bins are adjacent.  Invalid bins should never */
/* occur during normal operation.  They could be caused by invalid settings, */
/* doubleing point arithmetic errors, or another unforeseen bug. */
static int bin_diff(int min_bin,
                    int max_bin){
  int diff = abs(max_bin - min_bin);
  /* If the maximum bin isn't within one of the minimum bin, */
  /* the bins are invalid, so return an error and exit. */
  if(diff > 1){
    printf("\nCoordinate overlaps 3 or more bins: raise bin_size_ratio.\n");
    exit(EXIT_FAILURE);
  }
  return diff;
}

/* Recursive function that iterates over the Features in the given Pattern and */
/* duplicates the Pattern across the catalog given ambiguities in bin placement. */
/* Called by disambiguate_largest_edge, and calls disambiguate_rotation. */
static void disambiguate_bins(int just_count,
                              uint64_t *num_patterns,
                              Pattern pattern_catalog_cache[],
                              uint64_t *catalog_size_in_patterns,
                              int *cache_index,
                              Pattern new_pattern,
                              int feature_index){
  /* Before the Pattern can have its Features' bins disambiguated, it must first */
  /* have its fixed stars' largest edge bin disambiguated. */
  if(feature_index < 0){
    /* Iterate over all possible largest edge bin ambiguities. */
    int min_le_bin = bin_largest_edge(new_pattern.largest_edge, -1);
    int max_le_bin = bin_largest_edge(new_pattern.largest_edge, 1);
    for(int le_bin_offset=0;
        le_bin_offset<=bin_diff(min_le_bin, max_le_bin);
        le_bin_offset++){
      /* Set the Pattern's largest edge bin offset value. */
      new_pattern.le_bin_offset = le_bin_offset;
      /* Recursively calls itself to disambiguate the next bin. */
      disambiguate_bins(just_count,
                        num_patterns,
                        pattern_catalog_cache,
                        catalog_size_in_patterns,
                        cache_index,
                        new_pattern,
                        feature_index+1);
    }
  }
  /* If the Pattern hasn't had all of its coordinate bins disambiguated, */
  /* recurse down to select a bin for the next Feature in the Pattern. */
  else if(feature_index < num_stars_in_pattern-2){
    /* Compute largest edge bin for use in creating x and y coordinate bins. */
    unsigned int le_bin = bin_largest_edge(new_pattern.largest_edge,
                                           2*new_pattern.le_bin_offset-1);
    /* Iterate over all possible bin ambiguities in the x and y coordinate bins. */
    int min_y_bin = bin_y(new_pattern.features[feature_index].y, le_bin, -1);
    int max_y_bin = bin_y(new_pattern.features[feature_index].y, le_bin, 1);
    for(int y_bin_offset=0;
        y_bin_offset<=bin_diff(min_y_bin, max_y_bin);
        y_bin_offset++){
      /* Set the Feature's y bin offset. */
      new_pattern.features[feature_index].y_bin_offset = y_bin_offset;
      int y_bin = bin_y(new_pattern.features[feature_index].y, le_bin, 2*y_bin_offset-1);
      int min_x_bin = bin_x(new_pattern.features[feature_index].x, le_bin, y_bin, -1);
      int max_x_bin = bin_x(new_pattern.features[feature_index].x, le_bin, y_bin, 1);
      for(int x_bin_offset=0;
          x_bin_offset<=bin_diff(min_x_bin, max_x_bin);
          x_bin_offset++){
        /* Set the Feature's x bin offset. */
        new_pattern.features[feature_index].x_bin_offset = x_bin_offset;
        /* Recursively calls itself to disambiguate the next Feature. */
        disambiguate_bins(just_count,
                          num_patterns,
                          pattern_catalog_cache,
                          catalog_size_in_patterns,
                          cache_index,
                          new_pattern,
                          feature_index+1);
      }
    }
  }
  /* Once the Pattern has had its largest edge bin and coordinate bins disambiguated, */
  /* pass the Pattern on to disambiguate_rotation for insertion into the catalog. */
  else{
    /* Disambiguate the Pattern's 180 degree rotation, sort the Pattern's features, */
    /* and insert the fully disambiguated Pattern into the catalog. */
    disambiguate_rotation(just_count,
                          num_patterns,
                          pattern_catalog_cache,
                          catalog_size_in_patterns,
                          cache_index,
                          new_pattern);
  }
}

/* Duplicates the Pattern across the catalog given ambiguities in the largest edge. */
/* Called by iterate_patterns, and calls disambiguate_bins. */
static void disambiguate_largest_edge(Star stars[],
                                      int star_indices[num_stars_in_pattern],
                                      int just_count,
                                      uint64_t *num_patterns,
                                      Pattern pattern_catalog_cache[],
                                      uint64_t *catalog_size_in_patterns,
                                      int *cache_index){
  /* Dummy variables used for iteration. */
  int i,j,k;
  /* Pattern instance to be inserted into catalog. */
  Pattern new_pattern;
  /* Iterate over all pairs of stars in the pattern to */
  /* find the length of the largest edge. */
  double largest_edge_length = 0.0;
  for(i=0;i<num_stars_in_pattern;i++){
    for(j=i+1;j<num_stars_in_pattern;j++){
      /* Track the largest edge so far by comparing once with each edge. */
      largest_edge_length = fmax(largest_edge_length,
                                 dist(stars[star_indices[i]].vec,
                                      stars[star_indices[j]].vec));
    }
  }
  /* Disambiguate largest edge length by considering all choices of star */
  /* pairs within the pattern that are within error range of the largest edge. */
  for(i=0;i<num_stars_in_pattern;i++){
    for(j=i+1;j<num_stars_in_pattern;j++){
      /* Test ambiguity of each edge length with the largest edge. */
      double edge_length = dist(stars[star_indices[i]].vec,stars[star_indices[j]].vec);
      /* The error in each of the two stars sharing the edge may contribute up to */
      /* le_error_offset to the length of the edge.  An edge is therefore considered */
      /* ambiguous if its length is within 2*le_error_offset of the largest edge. */
      if(edge_length >= largest_edge_length-2*le_error_offset){
        /* Set the pattern's fixed star ids to those of the ambiguous edge. */
        /* Their order may be swapped later to resolve the 180 degree ambiguity. */
        new_pattern.fixed_star_id_1 = stars[star_indices[i]].star_id;
        new_pattern.fixed_star_id_2 = stars[star_indices[j]].star_id;
        /* Set the Pattern's largest edge length with the ambiguous edge, as it */
        /* will be used to determine the coordinates of the remaining stars. */
        /* Implicitly encodes the range from 0 up to the sine of the maximum */
        /* camera FOV as a 16 bit unsigned integer to keep the catalog compact. */
        new_pattern.largest_edge = (edge_length/max_le_length)*((1<<16)-1);
        /* Calculate vector along x axis of Pattern's coordinate system. */
        /* The vector points from Pattern star i to Pattern star j. */
        double x_axis_vector[3];
        diff(stars[star_indices[j]].vec,stars[star_indices[i]].vec,x_axis_vector);
        /* Calculate vector along y axis of Pattern's coordinate system. */
        double y_axis_vector[3];
        cross_prod(stars[star_indices[j]].vec,stars[star_indices[i]].vec,y_axis_vector);
        /* Normalize axis vectors to unit length by dividing by their magnitudes. */
        normalize(x_axis_vector);
        normalize(y_axis_vector);
        /* Use the remaining stars to initialize the Pattern's Features. */
        int feature_index = 0;
        for(k=0;k<num_stars_in_pattern;k++){
          /* Skip the fixed star ids, as they don't have their own Features. */
          if(k != i && k != j){
            /* Set the Feature's star id to match its corresponding star. */
            new_pattern.features[feature_index].star_id = stars[star_indices[k]].star_id;
            /* Calculate the normalized x and y coordinates using vector projection. */
            double x = dot_prod(x_axis_vector, stars[star_indices[k]].vec)/edge_length;
            double y = dot_prod(y_axis_vector, stars[star_indices[k]].vec)/edge_length;
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
        /* Disambiguate bins and 180 degree rotation before insertion into catalog. */
        disambiguate_bins(just_count,
                          num_patterns,
                          pattern_catalog_cache,
                          catalog_size_in_patterns,
                          cache_index,
                          new_pattern,
                          -1);
      }
    }
  }
}

/* Recursive function that iterates over all Patterns with largest edges */
/* smaller than the maximum Field of View.  If set to "just count," each */
/* valid Pattern is only used to increment a counter of valid Patterns. */
/* Otherwise, the function uses the Patterns to populate the catalog. */
/* just_count is a boolean that tells the function only to count */
/* the number of patterns and not to add them to the catalog. */
static void iterate_patterns(Star stars[],
                             int num_stars,
                             int just_count,
                             uint64_t *num_patterns,
                             Pattern pattern_catalog_cache[],
                             uint64_t *catalog_size_in_patterns,
                             int *cache_index,
                             int star_indices_index){
  /* Array of star ids for a given Pattern.  Being static avoids needing */
  /* to pass the array between recursive calls of the function. */
  static int star_indices[num_stars_in_pattern];
  /* For the given Pattern star, iterate over its possible star ids, */
  /* starting with the star id one greater than the Pattern star before */
  /* it in the list or 0 if it's the first Pattern star. */
  for(star_indices[star_indices_index]=star_indices_index>0?star_indices[star_indices_index-1]+1:0;
      star_indices[star_indices_index]<num_stars;
      star_indices[star_indices_index]++){
    /* If the first Pattern star is being changed, update the progress percentage. */
    if(star_indices_index == 0){
      /* Progress is measured as the percentage of star ids the first Pattern star */
      /* has covered.  Since Patterns are produced in sorted order by star id, later */
      /* star ids have fewer Patterns to iterate over, progress speeds up as it goes. */
      printf("\r%.2f%%", 100.0*star_indices[0]/(num_stars-1));
    }
    /* Check whether the selected star id is within the maximum Field of View of */
    /* the previously selected stars in the Pattern.  If any of the previously */
    /* selected stars are "out of range" of the newly selected star id, the */
    /* Pattern will be too big to fit within the Field of View, so skip over */
    /* the newly selected star completely and go to the next one. */
    int within_range = 1;
    for(int i=0;i<star_indices_index;i++){
      if(dist(stars[star_indices[i]].vec,
              stars[star_indices[star_indices_index]].vec) > max_le_length){
        within_range = 0;
        break;
      }
    }
    if(within_range){
      /* If the Pattern is valid but hasn't had all of its star ids selected, */
      /* recurse down to select a star id for the next star in the Pattern. */
      if(star_indices_index < num_stars_in_pattern-1){
        iterate_patterns(stars,
                         num_stars,
                         just_count,
                         num_patterns,
                         pattern_catalog_cache,
                         catalog_size_in_patterns,
                         cache_index,
                         star_indices_index+1);
      }
      /* If the Pattern is valid and has all of its star ids selected, call the */
      /* disambiguation functions to insert the Pattern into the catalog, duplicating */
      /* the Pattern based on its possible ambiguities.  The possible ambiguities */
      /* are in the largest edge length, bin placement, and 180 degree rotation. */
      else{
        disambiguate_largest_edge(stars,
                                  star_indices,
                                  just_count,
                                  num_patterns,
                                  pattern_catalog_cache,
                                  catalog_size_in_patterns,
                                  cache_index);
      }
    }
  }
}

/* Helper function that calls the iterate_patterns function in order */
/* to count the number of Patterns that will appear in the catalog. */
static uint64_t count_patterns(Star stars[],
                               int num_stars){
  /* Initialize Pattern counter to zero. */
  uint64_t num_patterns = 0;
  /* The 1 tells iterate_patterns it should only count patterns */
  /* without adding them to the catalog. */
  /* The NULLs are stand-ins for pointers to the catalog size, */
  /* the unused catalog cache, and cache_index. */
  /* The 0 tells iterate_patterns it needs to select the first star. */
  iterate_patterns(stars,
                   num_stars,
                   1,
                   &num_patterns,
                   NULL,
                   NULL,
                   NULL,
                   0);
  return num_patterns;
}

/* Populate pattern catalog file with patterns in pieces.  The pieces are */
/* kept in memory in the catalog cache until they're completed.  Then they're */
/* written to the catalog file before the next piece is made.  The pieces */
/* are created in order from the beginning of the catalog to the end. */
static void populate_catalog(FILE *pattern_catalog_file,
                             uint64_t catalog_size_in_patterns,
                             Star stars[],
                             int num_stars){
  /* Allocate memory for catalog cache. */
  /* By generating the catalog in pieces, it can be constructed in memory, */
  /* which speeds up the process significantly even for large catalogs. */
  Pattern *pattern_catalog_cache = calloc(max_cat_cache_size+max_probe_depth, sizeof(Pattern));
  /* Initialize catalog cache to its maximum value.  It is only */
  /* changed when generating the last piece of the catalog. */
  uint64_t cat_cache_size = max_cat_cache_size;
  /* Generate the catalog one piece at a time.  Cache_index tracks */
  /* which piece of the catalog is currently being generated. */
  int cache_index;
  /* Calculate number of cached pieces that will be needed to build entire catalog. */
  int num_catalog_pieces = catalog_size_in_patterns/max_cat_cache_size+1;
  /* Iterate over the pieces.  The last piece is truncated at the end of the catalog. */
  for (cache_index=0;cache_index<num_catalog_pieces;cache_index++) {
    printf("Generating catalog piece %d of %d...\n", cache_index+1, num_catalog_pieces);
    /* When generating the last piece of the catalog, resize the cache. */
    if (cache_index == num_catalog_pieces - 1){
      cat_cache_size = catalog_size_in_patterns%max_cat_cache_size;
    }
    /* Copy the overlapping max_probe_depth catalog locations from the previous cache. */
    /* This maintains the catalog's continuity despite generating it in pieces. */
    memcpy(pattern_catalog_cache,
           pattern_catalog_cache+max_cat_cache_size,
           max_probe_depth*sizeof(Pattern));
    /* Initialize this piece of the catalog by setting the rest to all zeros. */
    memset(pattern_catalog_cache+max_probe_depth,
           0,
           max_cat_cache_size*sizeof(Pattern));
    /* Populate catalog cache by iterating over all patterns, inserting only those which */
    /* have initial offsets which place them within the cached section of the catalog. */
    /* The NULL is a stand-in for the unused pointer to the Pattern counter. */
    iterate_patterns(stars,
                     num_stars,
                     0,
                     NULL,
                     pattern_catalog_cache,
                     &catalog_size_in_patterns,
                     &cache_index,
                     0);
    /* Write cached piece of the catalog to disk. */
    printf("\nWriting piece to disk...\n");
    fseeko64(pattern_catalog_file,
              max_cat_cache_size*cache_index*sizeof(Pattern),
              SEEK_SET);
    fwrite(pattern_catalog_cache,
           sizeof(Pattern),
           cat_cache_size+max_probe_depth,
           pattern_catalog_file);
  }
}

int main(int argc, char *argv[]) {
  /* Hipparchos catalog, pattern catalog, and star vectors file pointers. */
  FILE *bsc_file,*pattern_catalog_file,*stars_file,*hip_numbers_file;
  /* BSC5 cache */
  BSC_Entry bsc_cache[STARN];
  /* Allocate array of normalized star vectors, magnitudes, and star ids */
  /* in x,y,z,mag,id format.  Double stars have been removed from this array. */
  Star *stars_temp = malloc(STARN*sizeof(Star));
  Star *stars = malloc(STARN*sizeof(Star));
  /* Number of vectors/stars in stars array. */
  int num_stars_temp = 0;
  int num_stars = 0;
  /* Right ascension, Declination of star after correction to the current year. */
  double ra, dec;
  /* Number of Patterns in output catalog.  Stars in these Patterns are non-double stars */
  /* satisfying both the minimum brightness contraint and the max_fov constraint. */
  uint64_t num_patterns = 0;
  /* Size of pattern catalog in Patterns. */
  uint64_t catalog_size_in_patterns;
  /* Size of pattern catalog in bytes. */
  uint64_t catalog_size_in_bytes;

  /* Open BSC5 file for reading. */
  bsc_file = fopen("BSC5","rb");
  if(!bsc_file){
    printf("Unable to open BSC5 file!");
    return 1;
  }
  /* Offset read position in file by 28 bytes to ignore header data. */
  fseeko64(bsc_file, 28, SEEK_SET);
  /* Read BSC5 into cache. */
  for(int i=0;i<STARN;i++){
    fread(&bsc_cache[i].XNO, 4, 1, bsc_file);
    fread(&bsc_cache[i].SRA0, 8, 1, bsc_file);
    fread(&bsc_cache[i].SDEC0, 8, 1, bsc_file);
    fread(&bsc_cache[i].IS, 2, 1, bsc_file);
    fread(&bsc_cache[i].MAG, 2, 1, bsc_file);
    fread(&bsc_cache[i].XRPM, 4, 1, bsc_file);
    fread(&bsc_cache[i].XDPM, 4, 1, bsc_file);
  }
  /* Close BSC5 file. */
  fclose(bsc_file);

  /* Read BSC5 catalog and create temporary array of normalized vectors pointing */
  /* at each star.  Also filter out stars dimmer than the minimum magnitude. */
  for(int i=0;i<STARN;i++){
    /* If the star is at least as bright as the minimum magnitude, add its vector to the array. */
    if(bsc_cache[i].MAG/100.0 <= min_magnitude){
      /* Correct right ascension by adding proper motion times the number of years since 1950. */
      double ra = bsc_cache[i].SRA0+bsc_cache[i].XRPM*(current_year-1950);
      /* Correct declination by adding proper motion times the number of years since 1950 */
      double dec = bsc_cache[i].SDEC0+bsc_cache[i].XDPM*(current_year-1950);
      if(ra == 0.0 && dec == 0.0){
        continue;
      }
      double magnitude = bsc_cache[i].MAG/100.0;
      /* x component (forward) */
      stars_temp[num_stars_temp].vec[0] = cos(ra)*cos(dec);
      /* y component (left) */
      stars_temp[num_stars_temp].vec[1] = sin(ra)*cos(dec);
      /* z component (up) */
      stars_temp[num_stars_temp].vec[2] = sin(dec);
      num_stars_temp += 1;
    }
  }
  /* Filter double stars out of star vectors array. */
  for(int i=0;i<num_stars_temp;i++){
    int is_double_star = 0;
    for(int j=0;j<num_stars_temp;j++){
      if (j == i) {
        continue;
      }
      /* If the star vector is too close to any other star vector, */
      /* don't add it to the vector array. */
      if(dot_prod(stars_temp[i].vec,stars_temp[j].vec) > cos(min_star_separation)){
        is_double_star = 1;
        break;
      }
    }
    /* If it isn't too close to any other star vector, add it to the vector array. */
    if(!is_double_star) {
      stars[num_stars] = stars_temp[i];
      stars[num_stars].star_id = num_stars;
      num_stars += 1;
    }
  }

  /* Close Hipparchos catalog file. */
  fclose(bsc_file);
  /* Print the number of non-double stars at least as bright as the minimum magnitude. */
  printf("%d non-double stars at least magnitude %.2f found\n", num_stars, min_magnitude);

  /* Open stars file for writing. */
  stars_file = fopen("stars", "wb+");
  if(!stars_file){
    printf("Unable to open stars file!");
    return 1;
  }
  /* Generate stars file, which stores star vectors, magnitudes, and star ids. */
  printf("Inserting the %d stars into stars file...\n", num_stars);
  /* Write array of stars to the stars file. */
  fwrite(stars, sizeof(Star), num_stars, stars_file);
  /* Close stars file. */
  fclose(stars_file);

  /* Count number of Patterns that will be inserted into catalog. */
  /* The number of Patterns is needed prior to catalog generation */
  /* to set the catalog size based on the user specified density. */
  printf("Counting the number of Patterns the catalog will contain...\n");
  num_patterns = count_patterns(stars, num_stars);
  printf("\n%llu patterns counted\n", num_patterns);

  catalog_size_in_patterns = num_patterns*(1.0/catalog_density);
  catalog_size_in_bytes = catalog_size_in_patterns * sizeof(Pattern);
  printf("Generating catalog of size %llu patterns and %llu bytes...\n",
         catalog_size_in_patterns,
         catalog_size_in_bytes);
  /* Open pattern catalog file for writing. */
  pattern_catalog_file = fopen("pattern_catalog", "wb+");
  if(!pattern_catalog_file){
    printf("Unable to open pattern catalog file!");
    return 1;
  }
  /* Populate pattern catalog file with patterns. */
  populate_catalog(pattern_catalog_file,
                   catalog_size_in_patterns,
                   stars,
                   num_stars);
  /* Retrieve deepest probe depth reached while generating catalog. */
  uint64_t cache_offset = 0;
  int probe_step = 0;
  printf("Deepest probe depth: %d\n", increment_offset(&cache_offset, &probe_step));
  printf("Catalog generation complete!");
  return 0;
};
