/*
  construct_sym_graph.c: assumes coincident source and target locations
  and generates the adaptive tree and associated lists used in the new-version fmm.
  Copyright (c) 2012 Bo Zhang, Jingfang Huang, Nikos P. Pitsianis, Xiaobai Sun

  This program is free software: you can redistribute it and/or modify 
  it under the terms of the GNU General Public Licenses as published by 
  the Free Software Foundation, either version 3 of the Licenses, or any 
  later version. 

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of 
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program. If not, see http://www.gnu.org/licenses/.
 */

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "fmm_utils.h"
#include "fmm_ds.h"
#include "fmm_global.h"

int create_all_box(const double *ploc, const int nparts, const int s, 
		   int *perm, int *nboxes, int *nlev, double *size, 
		   fmmbox **boxes, int **content)
{
  // create_all_box assumes that source and target ensembles are the same.
  // It partitions a box with more than s points equally along each direction.
  int iter, L;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  const int xshift[] = {-1,1,-1,1,-1,1,-1,1};
  const int yshift[] = {-1,-1,1,1,-1,-1,1,1};
  const int zshift[] = {-1,-1,-1,-1,1,1,1,1};

  // setup root box
  xmin = DBL_MAX; ymin = xmin; zmin = xmin;
  xmax = -xmin; ymax = xmax; zmax = xmax;

  for ( iter = 0; iter < nparts; iter++ ) {
    double px = ploc[3*iter];
    double py = ploc[3*iter+1];
    double pz = ploc[3*iter+2];
    xmin = ( xmin < px? xmin: px );
    xmax = ( xmax > px? xmax: px );
    ymin = ( ymin < py? ymin: py );
    ymax = ( ymax > py? ymax: py );
    zmin = ( zmin < pz? zmin: pz );
    zmax = ( zmax > pz? zmax: pz );
  }

  double sidex = xmax - xmin;
  double sidey = ymax - ymin;
  double sidez = zmax - zmin;

  *size = MAX(sidex, sidey);
  *size = MAX(*size, sidez);

  double bcenter[3];
  bcenter[0] = ( xmax + xmin )*0.5;
  bcenter[1] = ( ymax + ymin )*0.5;
  bcenter[2] = ( zmax + zmin )*0.5;

  // initialize perm with identity mapping
  for ( iter = 0; iter <= nparts; iter++ ) 
    perm[iter] = iter;

  // allocate working space to perform partition
  const int maxlev = 200;
  fmmbox **work = (fmmbox **)calloc(maxlev, sizeof(fmmbox *));
  int *nnodes = (int *)calloc(maxlev, sizeof(int));
  if ( work == 0 || nnodes == 0) {
    printf("Error in %s, line %d: unable to allocate memory\n", __FILE__, __LINE__);
    return -1;
  }

  // setup root level 
  *nboxes = 1; 
  work[0] = (fmmbox *)calloc(1, sizeof(fmmbox));
  if ( work[0] == 0 ) {
    printf("Error in %s, line %d: unable to allocate memory\n", __FILE__, __LINE__);
    return -1;
  }
  work[0][0] = (fmmbox){{bcenter[0],bcenter[1],bcenter[2]},0,1,0,0, 
			{0,0,0,0,0,0,0,0},nparts,1};
  nnodes[0] = 1;

  int *swap = (int *)calloc(nparts, sizeof(int));
  if ( swap == 0 ) {
    printf("Error in %s, line %d: unable to allocate memory\n", __FILE__, __LINE__);
    return -1;
  }

  for ( L = 1; L < maxlev; L++ ) {
    int prevnboxes, box2divide, boxdivided;
    fmmbox *prev, *curr;

    prev = work[L-1];
    prevnboxes = nnodes[L-1];
    box2divide = 0;
    boxdivided = 0;

    // scan boxes at previous level, mark ones with more than s points for partition
    for ( iter = 0; iter < prevnboxes; iter++ ) {
      if ( prev[iter].npts > s ) {
	prev[iter].nchild = 1;
	box2divide += 1;
      }
    }
    
    if ( box2divide ) {
      int *addrs = (int *)calloc(8*box2divide, sizeof(int));
      int *counts = (int *)calloc(8*box2divide, sizeof(int));
      if ( addrs == 0 || counts == 0 ) {
	printf("Error in %s, line %d: unable to allocate memory\n", __FILE__, __LINE__);
	return -1;
      }
    
      for ( iter = 0; iter < prevnboxes; iter++ ) {
	if ( prev[iter].nchild == 1 ) {
	  int ptr = 8*boxdivided;
	  divide_box(prev[iter].center, ploc, perm+prev[iter].addr,
		     prev[iter].npts, swap, &addrs[ptr], &counts[ptr]);
	  boxdivided += 1;
	}
      }
  
      int boxcreated = 0;
      for ( iter = 0; iter < 8*box2divide; iter++) 
	boxcreated += ( counts[iter] > 0 );

      work[L] = (fmmbox *)calloc(boxcreated, sizeof(fmmbox));
      if ( work[L] == 0 ) {
	printf("Error in %s, line %d: unable to allocate memory\n", __FILE__, __LINE__);
	return -1;
      }
      
      curr = work[L];
      nnodes[L] = boxcreated;
      boxcreated = 0;
      boxdivided = 0;

      double h = *size/pow(2,L+1);
      for ( iter = 0; iter < prevnboxes; iter++ ) {
	if ( prev[iter].nchild == 1 ) { 
	  prev[iter].nchild = 0; 
	  int j; 
	  int offset = 8*boxdivided; 
	  for ( j = 0; j < 8; j++ ) {
	    if ( counts[offset+j] ) {
	      *nboxes += 1; 
	      curr[boxcreated].level = L; 
	      curr[boxcreated].boxid = *nboxes;
	      curr[boxcreated].parent = prev[iter].boxid; 
	      curr[boxcreated].npts = counts[offset+j];
	      curr[boxcreated].addr = prev[iter].addr + addrs[offset+j];
	      curr[boxcreated].center[0] = prev[iter].center[0] + xshift[j]*h;
	      curr[boxcreated].center[1] = prev[iter].center[1] + yshift[j]*h;
	      curr[boxcreated].center[2] = prev[iter].center[2] + zshift[j]*h;
	      prev[iter].nchild += 1; 
	      prev[iter].child[j] = curr[boxcreated].boxid;
	      boxcreated += 1;
	    }
	  }	  
	  boxdivided += 1;
	}
      }      
      free(addrs); 
      free(counts);
    } else {
      *nlev = L-1;
      break;
    }
  }

  if ( L == 200) {
    printf("Error in %s, line %d: too many levels attemped\n", __FILE__, __LINE__);
    return -1;
  }
  
  // write box information into one dimensional array
  fmmbox *_boxes = (fmmbox *)calloc(1+*nboxes, sizeof(fmmbox));
  int *_content = (int *)calloc((1+*nlev)*2, sizeof(int));
  if ( _boxes == 0 || _content == 0) {
    printf("Error in %s, line %d: unable to allocate memory\n", __FILE__, __LINE__);
    return -1;
  }

  for ( L = 0; L <= *nlev; L++ ) {
    int nnodes_L = nnodes[L];
    for ( iter = 0; iter < nnodes_L; iter++ ) {
      int j = work[L][iter].boxid;
      _boxes[j] = work[L][iter];
    }
    _content[2*L] = work[L][0].boxid;
    _content[2*L+1] = nnodes_L;
  }

  *boxes = _boxes;
  *content = _content;  


  for ( L = 0; L <= *nlev; L++ ) 
    free (work[L]);
  free(work);
  free(nnodes);

  free(swap);

  return 0;
}

void divide_box(double *bcenter, const double *ploc, int *start,
		const int npoints, int *swap, int *addrs, int *counts)
{
  int n1, n2, n3, n4, n5, n6, n7, n8;
  int n12, n34, n56, n78, n1234, n5678;
  int direction;
  double cutoff;

  n1 = 0; n2 = 0; n3 = 0; n4 = 0;
  n5 = 0; n6 = 0; n7 = 0; n8 = 0;
  n12 = 0; n34 = 0;
  n56 = 0; n78 = 0;
  n1234 = 0; 
  n5678 = 0;

  direction = 2;
  cutoff = bcenter[direction];
  assign_particle(ploc, start, npoints, direction,cutoff, swap, &n1234);
  n5678 = npoints - n1234;

  direction = 1;
  cutoff = bcenter[direction];
  if ( n1234 > 0 )
    assign_particle(ploc, start, n1234, direction, cutoff, swap, &n12);
  n34 = n1234 - n12;

  if ( n5678 > 0 )
    assign_particle(ploc, start+n1234, n5678, direction, cutoff, swap, &n56);		     
  n78 = n5678 - n56;

  direction = 0;
  cutoff = bcenter[direction];
  if ( n12 > 0 ) 
    assign_particle(ploc, start, n12, direction, cutoff, swap, &n1);
  n2 = n12 - n1;

  if ( n34 > 0 ) 
    assign_particle(ploc, start+n12, n34, direction, cutoff, swap, &n3);
  n4 = n34 - n3;

  if ( n56 > 0 ) 
    assign_particle(ploc, start+n1234, n56, direction, cutoff, swap, &n5);
  n6 = n56 - n5;

  if ( n78 > 0 ) 
    assign_particle(ploc, start+n1234+n56, n78, direction, cutoff, swap, &n7);
  n8 = n78 - n7;

  counts[0] = n1;
  counts[1] = n2;
  counts[2] = n3;
  counts[3] = n4;
  counts[4] = n5;
  counts[5] = n6;
  counts[6] = n7;
  counts[7] = n8;

  addrs[0] = 0;
  addrs[1] = addrs[0] + n1;
  addrs[2] = addrs[1] + n2;
  addrs[3] = addrs[2] + n3;
  addrs[4] = addrs[3] + n4;
  addrs[5] = addrs[4] + n5;
  addrs[6] = addrs[5] + n6;
  addrs[7] = addrs[6] + n7;
}

void assign_particle(const double *ploc, int *start, const int npoints,
		     const int direction, const double cutoff, int *swap, int *undercutoff)
{
  int stat1, stat2, partid, i;
  stat1 = 0;
  stat2 = 0;

  for ( i = 0; i < npoints; i++ ) {
    partid = start[i];
    if ( ploc[ 3*(partid-1) + direction ] <= cutoff ) {
      start[stat1] = partid;
      stat1 += 1;
    } else {
      swap[stat2] = partid;
      stat2 += 1;
    }
  }

  for ( i = 0; i < stat2; i++ ) 
    start[stat1 + i] = swap[i];

  *undercutoff = stat1;
}
 
void free_list(fmmlist *list, const int nboxes)
{
  int i;
  for ( i = 0; i <= nboxes; i++ ) {
    if ( list[i].list1 != 0 ) 
      free(list[i].list1);    
    if ( list[i].list3 != 0 ) 
      free(list[i].list3);   
    if ( list[i].list4 != 0 ) 
      free(list[i].list4);   
    if ( list[i].colleague != 0 ) 
      free(list[i].colleague);    
  }
  free(list);
}

int create_colleague(const fmmbox *boxes, const int *content, const int nlev, 
		 const double size, fmmlist *list)
{
  int i, j, k, m, boxnumber, listsize, parent, member, 
    box2sort, this, stat, stat5, indicator, current,  
    workzone[217], workzone5[28];

  // store colleague information of the root node
  list[1].colleague = (int *)calloc(2, sizeof(int));
  if ( list[1].colleague == 0 ) {
    printf("Error in %s, line %d: unable to allocate memory\n", 
	   __FILE__, __LINE__);
    return -1;
  }
  list[1].colleague[0] = 1;
  list[1].colleague[1] = 1;

  for ( i = 1; i <= nlev; i++ ) {
    current = content[2*i]; // id of the first box of level i
    boxnumber = content[2*i+1]; // number of boxes at level i

    // process each box at level i, generate its list 5
    for ( j = 1; j <= boxnumber; j++ ) {
      stat = 0; 
      stat5 = 0; 
      parent = boxes[current+j-1].parent;
      listsize = list[parent].colleague[0]; 

      /*
	check parent box's list 5, and any member who is a parent box 
	require further examiniation of its children
      */
      for ( k = 1; k <= listsize; k++ ) {
	member = list[parent].colleague[k];
	if ( boxes[member].nchild ) {
	  for ( m = 0; m < 8; m++ ) {
	    if ( boxes[member].child[m] > 0 ) {
	      stat += 1;
	      workzone[stat] = boxes[member].child[m];
	    }
	  }
	}
      }

      for ( m = 1; m <= stat; m++ ) {
	box2sort = workzone[m];       
	indicator = if_touch(size, i, boxes[current+j-1].center, i, boxes[box2sort].center);
	if ( indicator == 1 ) {
	  stat5 += 1;
	  workzone5[stat5] = box2sort;
	}
      }

      this = boxes[current+j-1].boxid;
      list[this].colleague = (int *)calloc(1+stat5, sizeof(int));
      if ( list[this].colleague == 0 ) {
	printf("Error in %s, line %d: unable to allocate memory\n", 
	       __FILE__, __LINE__);
	return -1;
      }

      for ( m = 1; m <= stat5; m++ ) 
	list[this].colleague[m] = workzone5[m];

      list[this].colleague[0] = stat5;
    }
  }

  return 0;
}

int if_touch(const double size, const int levela, const double *centera, 
	     const int levelb, const double *centerb)
{
  double sidea, halfsidea, sideb, halfsideb, eps;
  double xmina, ymina, zmina, xmaxa, ymaxa, zmaxa;
  double xminb, yminb, zminb, xmaxb, ymaxb, zmaxb;

  sidea = size/pow(2,levela);
  sideb = size/pow(2,levelb);

  halfsidea = sidea/2.0;
  halfsideb = sideb/2.0;

  xmina = centera[0] - halfsidea;
  xmaxa = centera[0] + halfsidea;

  ymina = centera[1] - halfsidea;
  ymaxa = centera[1] + halfsidea;

  zmina = centera[2] - halfsidea;
  zmaxa = centera[2] + halfsidea;

  xminb = centerb[0] - halfsideb;
  xmaxb = centerb[0] + halfsideb;

  yminb = centerb[1] - halfsideb;
  ymaxb = centerb[1] + halfsideb;

  zminb = centerb[2] - halfsideb;
  zmaxb = centerb[2] + halfsideb;


  eps = ( sidea < sideb ? sidea : sideb );

  eps /= 10000;

  if ( xmina > xmaxb + eps )
    return 0;

  if ( xminb > xmaxa + eps )
    return 0;

  if ( ymina > ymaxb + eps )
    return 0;

  if ( yminb > ymaxa + eps )
    return 0;

  if ( zmina > zmaxb + eps )
    return 0;

  if ( zminb > zmaxa + eps )
    return 0;

  return 1;
}

int create_list134(const fmmbox *boxes, const int *content, const int nlev, 
		   const int nboxes, const double size, fmmlist *list)
{
  int i, j, k, m, boxnumber, listsize, todivide, indicator, 
    current, this, box2proc;
  listnode *workzone1, *workzone3, *workzone4, head, *push, *pop, 
    *dual, *coarse, *child; 

  /* 
     we use a stack structure to sort the boxes into lists 1 and 3. 
     variable head is the top of the stack, whose data field is the number
     of entries in the stack
  */

  head.next = 0;
  head.data = 0; 

  // allocate working space
  workzone1 = (listnode *)calloc(1+nboxes, sizeof(listnode));
  workzone3 = (listnode *)calloc(1+nboxes, sizeof(listnode));
  workzone4 = (listnode *)calloc(1+nboxes, sizeof(listnode));

  if ( workzone1 == 0 || workzone3 == 0 || workzone4 == 0 ) {
    printf("Error in %s, line %d: unable to allocate memory\n", 
	   __FILE__, __LINE__);
    return -1;
  }


  for ( i = 0; i <= nlev; i++ ) {
    current = content[2*i]; // id of first box at level i
    boxnumber = content[2*i+1]; // number of boxes at level i

    for ( j = 1; j <= boxnumber; j++ ) {
      if ( boxes[current+j-1].nchild == 0 ) {
	this = boxes[current+j-1].boxid; 
	listsize = list[this].colleague[0];

	// process each box of the colleague. 
	for ( k = 1; k <= listsize; k++ ) {
	  todivide = list[this].colleague[k];
	  push = (listnode *)calloc(1, sizeof(listnode));
	  if ( push == 0 ) {
	    printf("Error in %s, line %d: unable to allocate memory\n",
		   __FILE__, __LINE__);
	    return -1;
	  }
	  push->data = todivide;
	  head.next = push;
	  head.data = 1;

	  while ( head.data>0 ) {
	    pop = head.next;
	    box2proc = pop->data;

	    if ( box2proc == this ) {
	      // each childless box contains itself in its list1
	      pop->next = workzone1[this].next;
	      workzone1[this].next = pop;
	      workzone1[this].data += 1;
	      break;
	    }

	    indicator = if_touch (size, i, boxes[current+j-1].center, 
				  boxes[box2proc].level, boxes[box2proc].center);

	    if ( indicator == 0 ) {
	      // the box is a list3 member
	      head.next = pop->next;
	      head.data -= 1;
	      pop->next = workzone3[this].next;
	      workzone3[this].next = pop;
	      workzone3[this].data += 1;

	      // update list4 due to its dual relationship with list3
	      dual = (listnode *)calloc(1, sizeof(listnode));
	      if ( dual == 0 ) {	
		printf("Error in %s, line %d: unable to allocate memory\n", 
		       __FILE__, __LINE__);
		return -1;
	      }
	      dual->data = this;
	      dual->next = workzone4[box2proc].next;
	      workzone4[box2proc].next = dual;
	      workzone4[box2proc].data += 1;

	    } else {
	      // two boxes touch
	      if ( boxes[box2proc].nchild == 0 ) {
		// then update the list1
		head.next = pop->next;
		head.data -= 1;

		if ( i < boxes[box2proc].level ) {
		  // the box resides at a finer level, update its list1 
		  coarse = (listnode *)calloc(1, sizeof(listnode));

		  if ( coarse == 0 ) {
		    printf("Error in %s, line %d: unable to allocate memory\n", 
			   __FILE__, __LINE__);
		    return -1;
		  }

		  coarse->next = workzone1[box2proc].next;
	          coarse->data = this;
		  workzone1[box2proc].next = coarse;
		  workzone1[box2proc].data += 1;
		}

		pop->next = workzone1[this].next;
		workzone1[this].next = pop;
		workzone1[this].data += 1;
	      } else {
		/*
		  we encounter a touching parent box. remove it from the stack 
		  and push its children in for examination
		*/
		head.next = pop->next;
		free(pop);
		head.data -= 1;

		for ( m = 0; m < 8; m++ ) {
		  if ( boxes[box2proc].child[m] > 0 ) {
		    child = (listnode *)calloc(1, sizeof(listnode));
		    if ( child == 0 ) {
		      printf("Error in %s, line %d: unable to allocate memory\n",
			     __FILE__, __LINE__);
		      return -1;
		    }
		    
		    child->data = boxes[box2proc].child[m];
		    child->next = head.next;
		    head.next = child;
		    head.data += 1;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }

  head.next = 0;
  head.data = 0;

  for ( i = 0; i <= nlev; i++ ) {
    current = content[2*i]; 
    boxnumber = content[2*i+1];

    for ( j = 1; j <= boxnumber; j++ ) {
      this = boxes[current+j-1].boxid;

      listsize = workzone1[this].data;
      if ( listsize > 0 ) {
	list[this].list1 = (int *)calloc(1+listsize, sizeof(int));
	if ( list[this].list1 == 0 ) {
	  printf("Error in %s, line %d: unable to allocate memory\n",
		 __FILE__, __LINE__);
	  return -1;
	}

	list[this].list1[0] = listsize;
	for ( k = listsize; k >= 1; k-- ) {
	  pop = workzone1[this].next;
	  workzone1[this].next = pop->next;
	  list[this].list1[k] = pop->data;
	  free(pop);
	}
      }


      listsize = workzone3[this].data;

      if ( listsize > 0 ) {
	list[this].list3 = (int *)calloc(1+listsize, sizeof(int));
	if ( list[this].list3 == 0 ) {
	  printf("Error in %s, line %d: unable to allocate memory\n",
		 __FILE__, __LINE__);
	  return -1;
	}

	list[this].list3[0] = listsize;
	for ( k = listsize; k >= 1; k-- ) {
	  pop = workzone3[this].next;
	  workzone3[this].next = pop->next;
	  list[this].list3[k] = pop->data;
	  free(pop);
	}
      }

      listsize = workzone4[this].data;
      if ( listsize > 0 ) {
	list[this].list4 = (int *)calloc(1+listsize, sizeof(int));	
	if ( list[this].list4 == 0 ) {
	  printf("Error in %s, line %d: unable to allocate memory\n",
		 __FILE__, __LINE__);
	  return -1;
	}

	list[this].list4[0] = listsize;
	for ( k = listsize; k >= 1; k-- ) {
	  pop = workzone4[this].next;
	  workzone4[this].next = pop->next;
	  list[this].list4[k] = pop->data;
	  free(pop);
	}
      }
    }
  }


  free(workzone1);
  free(workzone3);
  free(workzone4);

  return 0;

}

void create_list_up(const int mybox, int *uall, int *nuall, int *xuall, int *yuall, 
		    int *u1234, int *nu1234, int *x1234, int *y1234)
{
  int *colleague = &LIST[mybox].colleague[1];
  int ncolleague = LIST[mybox].colleague[0];
  int i, j, ix, iy, member, child;
  double side1, side2, temp; 

  *nuall = 0; 
  *nu1234 = 0; 

  for ( i = 0; i < 36; i++ ) {
    uall[i] = 0; 
    xuall[i] = 0; 
    yuall[i] = 0; 
  }

  for ( i = 0; i < 16; i++ ) {
    u1234[i] = 0; 
    x1234[i] = 0; 
    y1234[i] = 0; 
  }

  side1 = SIZE/pow(2,BOXES[mybox].level);
  side2 = side1*0.5;

  for ( i = 0; i < ncolleague; i++ ) {
    member = colleague[i];
    if ( BOXES[member].nchild > 0 && 
	 BOXES[member].center[2] - BOXES[mybox].center[2] > side2 ) {
      for ( j = 0; j < 8; j++ ) {
	child = BOXES[member].child[j];
	if ( child > 0 ) {
	  temp = BOXES[child].center[0] - BOXES[mybox].center[0] + side2/2;
	  temp = temp/side2;
	  ix = NINT(temp);
	  temp = BOXES[child].center[1] - BOXES[mybox].center[1] + side2/2;
	  temp = temp/side2;
	  iy = NINT(temp);

	  if ( BOXES[child].center[2] - BOXES[mybox].center[2] > side1 ) {
	    uall[*nuall] = child;
	    xuall[*nuall] = ix;
	    yuall[*nuall] = iy;
	    *nuall += 1;
	  } else {
	    if ( ix == -2 || ix == 3 || iy == -2 || iy ==3 ) {
	    } else {
	      u1234[*nu1234] = child;
	      x1234[*nu1234] = ix;
	      y1234[*nu1234] = iy;
	      *nu1234 += 1;
	    }
	  }
	}
      }
    }
  }
}

void create_list_down(const int mybox, int *dall, int *ndall, int *xdall, int *ydall, 
		      int *d5678, int *nd5678, int *x5678, int *y5678)
{
  int *colleague = &LIST[mybox].colleague[1];
  int ncolleague = LIST[mybox].colleague[0];
  int i, j, ix, iy, member, child;
  double side1, side2, temp;

  *ndall = 0;
  *nd5678 = 0;

  for ( i = 0; i < 36; i++ ) {
    dall[i] = 0;
    xdall[i] = 0;
    ydall[i] = 0;
  }

  for ( i = 0; i < 16; i++ ) {
    d5678[i] = 0;
    x5678[i] = 0;
    y5678[i] = 0;
  }

  side1 = SIZE/pow(2,BOXES[mybox].level);
  side2 = side1*0.5;

  for ( i = 0; i < ncolleague; i++ ) {
    member = colleague[i];
    if ( BOXES[member].nchild > 0 && 
	 BOXES[member].center[2] - BOXES[mybox].center[2] < -side2 ) {
      for ( j = 0; j < 8; j++ ) {
	child = BOXES[member].child[j];
	if ( child > 0 ) {
	  temp = BOXES[child].center[0] - BOXES[mybox].center[0] + side2/2;
	  temp = temp/side2;
	  ix = NINT(temp);
	  temp = BOXES[child].center[1] - BOXES[mybox].center[1] + side2/2;
	  temp = temp/side2;
	  iy = NINT(temp);

	  if ( BOXES[child].center[2] - BOXES[mybox].center[2] <= -side1 ) {
	    dall[*ndall] = child;
	    xdall[*ndall] = ix;
	    ydall[*ndall] = iy;
	    *ndall += 1;
	  } else {
	    if ( ix == -2 || ix == 3 || iy == -2 || iy == 3 ) {
	    } else {
	      d5678[*nd5678] = child;
	      x5678[*nd5678] = ix;
	      y5678[*nd5678] = iy;
	      *nd5678 += 1;
	    }
	  }
	}
      }
    }
  }
}


void create_list_north(const int mybox, int *nall, int *nnall, int *xnall, int *ynall,
		       int *n1256, int *nn1256, int *x1256, int *y1256,
		       int *n12, int *nn12, int *x12, int *y12,
		       int *n56, int *nn56, int *x56, int *y56)
{
  int *colleague = &LIST[mybox].colleague[1];
  int ncolleague = LIST[mybox].colleague[0];
  int i, j, ix, iy, member,child;
  double side1, side2, temp;

  *nnall = 0;
  *nn1256 = 0;
  *nn12 = 0;
  *nn56 = 0;

  for ( i = 0; i < 24; i++ ) {
    nall[i] = 0;
    xnall[i] = 0;
    ynall[i] = 0;
  }

  for ( i = 0; i < 8; i++ ) {
    n1256[i] = 0;
    x1256[i] = 0;
    y1256[i] = 0;
  }

  for ( i = 0; i < 4; i++ ) {
    n12[i] = 0;
    x12[i] = 0;
    y12[i] = 0;
  }

  for ( i = 0; i < 4; i++ ) {
    n56[i] = 0;
    x56[i] = 0;
    y56[i] = 0;
  }

  side1 = SIZE/pow(2, BOXES[mybox].level); 
  side2 = side1*0.5;

  for ( i = 0; i < ncolleague; i++ ) {
    member = colleague[i];
    if ( BOXES[member].nchild > 0 && 
	 BOXES[member].center[1] - BOXES[mybox].center[1] > side2 ) {
      for ( j = 0; j < 8; j++ ) {
	child = BOXES[member].child[j];
	if ( child > 0 ) {
	  temp = BOXES[child].center[2] - BOXES[mybox].center[2] + side2/2;
	  temp = temp/side2;
	  ix = NINT(temp);
	  temp = BOXES[child].center[0] - BOXES[mybox].center[0] + side2/2;
	  temp = temp/side2;
	  iy = NINT(temp);
 
	  if ( BOXES[child].center[1] - BOXES[mybox].center[1] > side1 ) {
	    if ( ix != -2 && ix != 3 ) {
	      nall[*nnall] = child;
	      xnall[*nnall] = ix;
	      ynall[*nnall] = iy;
	      *nnall += 1;
	    }
	  } else {
	    if ( (ix == 0 || ix == 1) && iy >= -1 && iy <= 2 ) {
	      n1256[*nn1256] = child;
	      x1256[*nn1256] = ix;
	      y1256[*nn1256] = iy;
	      *nn1256 += 1;
	    } else if ( ix == -1 && (iy >= -1 && iy <= 2) ) {
	      n12[*nn12] = child;
	      x12[*nn12] = ix;
	      y12[*nn12] = iy;
	      *nn12 += 1;
	    } else if ( ix == 2 && (iy >= -1 && iy <= 2) ) {
	      n56[*nn56] = child;
	      x56[*nn56] = ix;
	      y56[*nn56] = iy;
	      *nn56 += 1;
	    }
	  }
	}
      }
    }
  }
}

void create_list_south(const int mybox, int *sall, int *nsall, int *xsall, int *ysall,
		       int *s3478, int *ns3478, int *x3478, int *y3478,
		       int *s34, int *ns34, int *x34, int *y34,
		       int *s78, int *ns78, int *x78, int *y78)
{
  int *colleague = &LIST[mybox].colleague[1];
  int ncolleague = LIST[mybox].colleague[0];
  int i, j, ix, iy, member, child;
  double side1, side2, temp;

  *nsall = 0;
  *ns3478 = 0;
  *ns34 = 0;
  *ns78 = 0;

  for ( i = 0; i < 24; i++ ) {
    sall[i] = 0;
    xsall[i] = 0;
    ysall[i] = 0;
  }

  for ( i = 0; i < 8; i++ ) {
    s3478[i] = 0;
    x3478[i] = 0;
    y3478[i] = 0;
  }

  for ( i = 0; i < 4; i++ ) {
    s34[i] = 0;
    x34[i] = 0;
    y34[i] = 0;
  }

  for ( i = 0; i < 4; i++ ) {
    s78[i] = 0;
    x78[i] = 0;
    y78[i] = 0;
  }

  side1 = SIZE/pow(2, BOXES[mybox].level);
  side2 = side1 * 0.5;

  for ( i = 0; i < ncolleague; i++ ) {
    member = colleague[i];
    if ( BOXES[member].nchild > 0 &&
	 BOXES[member].center[1] - BOXES[mybox].center[1] < -side2 ) {      
      for ( j = 0; j < 8; j++ ) {
	child = BOXES[member].child[j];
	if ( child > 0 ) {
	  temp = BOXES[child].center[2] - BOXES[mybox].center[2] + side2/2;
	  temp = temp/side2;
	  ix = NINT(temp);
	  temp = BOXES[child].center[0] - BOXES[mybox].center[0] + side2/2;
	  temp = temp/side2;
	  iy = NINT(temp);

	  if ( BOXES[child].center[1] - BOXES[mybox].center[1] < -side1 ) {
	    if ( ix != -2 && ix != 3 ) {
	      sall[*nsall] = child;
	      xsall[*nsall] = ix;
	      ysall[*nsall] = iy;
	      *nsall += 1;
	    }	    
	  } else {	    
	    if ( (ix == 0 || ix == 1) && iy >= -1 && iy <= 2) {
	      s3478[*ns3478] = child;
	      x3478[*ns3478] = ix;
	      y3478[*ns3478] = iy;
	      *ns3478 += 1;
	    } else if ( ix == -1 && (iy >= -1 && iy <= 2) ) {
	      s34[*ns34] = child;
	      x34[*ns34] = ix;
	      y34[*ns34] = iy;
	      *ns34 += 1;	      
	    } else if ( ix == 2 && (iy >= -1 && iy <= 2) ) {
	      s78[*ns78] = child;
	      x78[*ns78] = ix;
	      y78[*ns78] = iy;
	      *ns78 += 1;	      
	    }
	  }
	}
      }
    }
  }
}

void create_list_east(const int mybox, int *eall, int *neall, int *xeall, int *yeall,
		      int *e1357, int *ne1357, int *x1357, int *y1357,
		      int *e13, int *ne13, int *x13, int *y13,
		      int *e57, int *ne57, int *x57, int *y57,
		      int *e1, int *ne1, int *x1, int *y1,
		      int *e3, int *ne3, int *x3, int *y3,
		      int *e5, int *ne5, int *x5, int *y5,
		      int *e7, int *ne7, int *x7, int *y7)
{
  int *colleague = &LIST[mybox].colleague[1];
  int ncolleague =LIST[mybox].colleague[0];
  int i, j, ix, iy, member, child;
  double side1, side2, temp;

  *neall = 0;
  *ne1357 = 0;
  *ne13 = 0; *ne57 = 0;
  *ne1 = 0; *ne3 = 0; *ne5 = 0; *ne7 = 0;

  for ( i = 0; i < 16; i++ ) {
    eall[i] = 0;
    xeall[i] = 0;
    yeall[i] = 0;
  }

  for ( i = 0; i < 4; i++ ) {
    e1357[i] = 0;
    x1357[i] = 0;
    y1357[i] = 0;

    e13[i] = 0;
    x13[i] = 0;
    y13[i] = 0;

    e57[i] = 0;
    x57[i] = 0;
    y57[i] = 0;

    e1[i] = 0;
    x1[i] = 0;
    y1[i] = 0;

    e3[i] = 0;
    x3[i] = 0;
    y3[i] = 0;

    e5[i] = 0;
    x5[i] = 0;
    y5[i] = 0;

    e7[i] = 0;
    x7[i] = 0;
    y7[i] = 0;
  }

  side1 = SIZE/pow(2,BOXES[mybox].level);
  side2 = side1 * 0.5;

  for ( i = 0; i < ncolleague; i++ ) {
    member = colleague[i];
    if ( BOXES[member].nchild > 0 &&
	 BOXES[member].center[0] - BOXES[mybox].center[0] > side2 ) {
      for ( j = 0; j < 8; j++ ) {
	child = BOXES[member].child[j];
	if ( child > 0 ) {
	  temp = BOXES[child].center[2] - BOXES[mybox].center[2] + side2/2;
	  temp = temp/side2;	  
	  ix = NINT(temp);
	  temp = BOXES[child].center[1] - BOXES[mybox].center[1] + side2/2;
	  temp = temp/side2;	     
	  iy = NINT(temp);

	  if ( BOXES[child].center[0] - BOXES[mybox].center[0] > side1 ) {
	    if ( ix >= -1 && ix <= 2 && iy >= -1 && iy <= 2 ) {
	      eall[*neall] = child;
	      xeall[*neall] = -ix;
	      yeall[*neall] = iy;
	      *neall += 1;	      
	    }
	  } else {
	    if ( (ix == 0 || ix == 1) && (iy == 0 || iy == 1) ) {
	      e1357[*ne1357] = child;
	      x1357[*ne1357] = -ix;
	      y1357[*ne1357] = iy;
	      *ne1357 += 1;	     
	    } else if ( ix == -1 && (iy == 0 || iy == 1) ) {
	      e13[*ne13] = child;
	      x13[*ne13] = -ix;
	      y13[*ne13] = iy;
	      *ne13 += 1;	      
	    } else if ( ix == 2 && (iy == 0 || iy == 1) ) {
	      e57[*ne57] = child;
	      x57[*ne57] = -ix;
	      y57[*ne57] = iy;
	      *ne57 += 1;
	    } else if ( iy == -1 ) {
	      if ( ix >= -1 && ix <= 1 ) {
		e1[*ne1] = child;
		x1[*ne1] = -ix;
		y1[*ne1] = iy;
		*ne1 += 1;		
	      }

	      if ( ix >= 0 && ix <= 2 ) {
		e5[*ne5] = child;
		x5[*ne5] = -ix;
		y5[*ne5] = iy;
		*ne5 += 1;		
	      }	      
	    } else if ( iy == 2 ) {	      
	      if ( ix >= -1 && ix <= 1) {
		e3[*ne3] = child;
		x3[*ne3] = -ix;
		y3[*ne3] = iy;
		*ne3 += 1;		
	      }

	      if ( ix >= 0 && ix <= 2) {
		e7[*ne7] = child;
		x7[*ne7] = -ix;
		y7[*ne7] = iy;
		*ne7 += 1;
	      }
	    }
	  }
	}
      }
    }
  }
}


void create_list_west(const int mybox, int *wall, int *nwall, int *xwall, int *ywall,
		      int *w2468, int *nw2468, int *x2468, int *y2468,
		      int *w24, int *nw24, int *x24, int *y24,
		      int *w68, int *nw68, int *x68, int *y68,
		      int *w2, int *nw2, int *x2, int *y2,
		      int *w4, int *nw4, int *x4, int *y4,
		      int *w6, int *nw6, int *x6, int *y6,
		      int *w8, int *nw8, int *x8, int *y8)
{
  int *colleague = &LIST[mybox].colleague[1];
  int ncolleague = LIST[mybox].colleague[0];
  int i, j, ix, iy, member, child;
  double side1, side2, temp;

  *nwall = 0;
  *nw2468 = 0;
  *nw24 = 0; *nw68 = 0;
  *nw2 = 0; *nw4 = 0; *nw6 = 0; *nw8 = 0;

  for ( i = 0; i < 16; i++ ) {
    wall[i] = 0;
    xwall[i] = 0;
    ywall[i] = 0;
  }

  for ( i = 0; i < 4; i++ ) {
    w2468[i] = 0;
    x2468[i] = 0;
    y2468[i] = 0;

    w24[i] = 0;
    x24[i] = 0;
    y24[i] = 0;

    w68[i] = 0;
    x68[i] = 0;
    y68[i] = 0;

    w2[i] = 0;
    x2[i] = 0;
    y2[i] = 0;

    w4[i] = 0;
    x4[i] = 0;
    y4[i] = 0;

    w6[i] = 0;
    x6[i] = 0;
    y6[i] = 0;

    w8[i] = 0;
    x8[i] = 0;
    y8[i] = 0;
  }

  side1 = SIZE/pow(2,BOXES[mybox].level);
  side2 = side1/2;

  for ( i = 0; i < ncolleague; i++ ) {
    member = colleague[i];
    if ( BOXES[member].nchild > 0 &&
	 BOXES[member].center[0] - BOXES[mybox].center[0] < -side2 ) {
      for ( j = 0; j < 8; j++ ) {
	child = BOXES[member].child[j];
	if ( child > 0 ) {
	  temp = BOXES[child].center[2] - BOXES[mybox].center[2] + side2/2;
	  temp = temp/side2;
	  ix = NINT(temp);
	  temp = BOXES[child].center[1] - BOXES[mybox].center[1] + side2/2;
	  temp = temp/side2;	 
	  iy = NINT(temp);

	  if ( BOXES[child].center[0] - BOXES[mybox].center[0] < -side1 ) {
	    if ( ix >= -1 && ix <= 2 && iy >= -1 && iy <= 2 ) {
	      wall[*nwall] = child;
	      xwall[*nwall] = -ix;
	      ywall[*nwall] = iy;
	      *nwall += 1;	      
	    }
	  } else {
	    if ( (ix == 0 || ix == 1) && (iy == 0 || iy == 1) ) {
	      w2468[*nw2468] = child;
	      x2468[*nw2468] = -ix;
	      y2468[*nw2468] = iy;
	      *nw2468 += 1;	   
	    } else if ( ix == -1 && (iy == 0 || iy == 1) ) {
	      w24[*nw24] = child;
	      x24[*nw24] = -ix;
	      y24[*nw24] = iy;
	      *nw24 += 1;	      
	    } else if ( ix == 2 && (iy == 0 || iy == 1) ) {
	      w68[*nw68] = child;
	      x68[*nw68] = -ix;
	      y68[*nw68] = iy;
	      *nw68 += 1;	      
	    } else if ( iy == -1 ) {	      
	      if ( ix >= -1 && ix <= 1 ) {
		w2[*nw2] = child;
		x2[*nw2] = -ix;
		y2[*nw2] = iy;
		*nw2 += 1;		
	      }

	      if ( ix >= 0 && ix <= 2 ) {
		w6[*nw6] = child;
		x6[*nw6] = -ix;
		y6[*nw6] = iy;
		*nw6 += 1;		
	      }

	    } else if ( iy == 2 ) {
	      if ( ix >= -1 && ix <= 1 ) {
		w4[*nw4] = child;
		x4[*nw4] = -ix;
		y4[*nw4] = iy;
		*nw4 += 1;		
	      }

	      if ( ix >= 0 && ix <= 2 ) {
		w8[*nw8] = child;
		x8[*nw8] = -ix;
		y8[*nw8] = iy;
		*nw8 += 1;		
	      }
	    }
	  }
	}
      }
    }
  }
}

 
 
