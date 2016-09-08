# Tetra

##Try Tetra Out
Test Tetra before you download it with the [Web Demo](https://s3.amazonaws.com/startracking/index.html).

##Run Tetra on Your Personal Computer

1. Create a directory (i.e. Tetra).
2. Place tetra.py in the directory.
3. Place Yale's Bright Star Catalog in the directory: http://tdc-www.harvard.edu/catalogs/BSC5
4. Create a subdirectory called 'pics'
5. Place images in 'pics' such as this one: http://stars.astro.illinois.edu/sow/cas.jpg
6. Run 'python tetra.py'

The first run may take a while as it needs to generate the catalog.  From then on, the majority of the runtime will be taken up by loading the catalog into memory and image processing.
