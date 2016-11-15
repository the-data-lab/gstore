
/*
 * Copyright 2016 The George Washington University
 * Written by Pradeep Kumar 
 * Directed by Prof. Howie Huang
 *
 * https://www.seas.gwu.edu/~howie/
 * Contact: iheartgraph@gmail.com
 *
 * 
 * Please cite the following paper:
 * 
 * Pradeep Kumar and H. Howie Huang. 2016. G-Store: High-Performance Graph Store for Trillion-Edge Processing. In Proceedings of the International Conference for High Performance Computing, Networking, Storage and Analysis (SC '16).
 
 *
 * This file is part of G-Store.
 *
 * G-Store is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * G-Store is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with G-Store.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <omp.h>
#include <fstream>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <assert.h>
#include <algorithm>
#include <errno.h>

#define handle_error(msg) \
                   do { perror(msg); exit(EXIT_FAILURE); } while (0)

void *map_anon_memory(uint64_t size, 
			     bool mlocked,
			     bool zero = false)
{
  size = size >0 ? size:4096;
  void *space = mmap(NULL, size > 0 ? size:4096, 
		     PROT_READ|PROT_WRITE,
		     MAP_ANONYMOUS|MAP_SHARED, -1, 0);
  if(space == MAP_FAILED) {
    handle_error("map_anon_memory:mmap:");
    exit(-1);
  }
  if(mlocked) {
    if(mlock(space, size) < 0) {
    handle_error("map_anon_memory:mlock:");
    }
  }
  if(zero) {
    memset(space, 0, size);
  }
  return space;
}

int main()
{
		//uint64_t size = 33706807296L;
		//uint64_t size = 32006807296; //1 GB left
		  uint64_t size = 32206807296; //1 GB left
		//uint64_t size =   25184698368;  
		map_anon_memory(size, true);
		while(1) pause();
}
