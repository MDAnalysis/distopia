distopia
--------

Faster distance calculations for the year 2020.

To build::
  mkdir build
  cd build
  cmake ..
  make
  
To generate data::
  python generate_coords.py 1000000
  
Will generate 1000000 random coordinates.

To benchmark methods::
  Xdist data.txt