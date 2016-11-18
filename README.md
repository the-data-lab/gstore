# gstore
Semi-external graph engine for trillion-edge graphs.

==========Help file==========
We will be updating this file as and when required. 

Building:
  make gstored
  make gstoreu
  
gstored is for directed graphs.
gstoreu is for undirected graphs.

Things to tune:
gstore.h
  NUM_THDS : Total number of threads to spawn using openMP. We have set its value as 56. Change it prior to compilation.
  
gstore.cpp  
  total_size: total memory size that you want to set for streaming and caching. Default value is set as 8GB. There is a run time input for this variable.
  memory: Size of segment. There is two segment used for IO and processing alternatly. Make sure that you set memory and total_size in such a way that total_size/memory gives you an integer.
