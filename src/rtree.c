/*
 *
 * rtree.c 
 *
 * Copyright (c) 2013-2023 Devon M. Powell
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */

#include "rtree.h"

// private declarations, not exposed to the user
void psi_rtree_rbox(psi_rtree_node* node, psi_rvec* rbox); 
void psi_rtree_split(psi_rtree_node* newnodes[2], psi_rtree_data* entries); 
psi_int psi_rtree_min_area(psi_rtree_node* node, psi_rvec* rbox);
psi_int psi_rtree_min_overlap(psi_rtree_node* node, psi_rvec* rbox);
static psi_real* psi_rtree_sort_helper_data; // TODO: is the helper data threadsafe?
int psi_rtree_sort_helper(const void* p1, const void* p2);
void psi_rtree_insert_fancy(psi_rtree_data newnode, psi_int level, psi_int split, psi_rtree* rtree);

void psi_rtree_insert(psi_rtree* rtree, psi_rvec* rbox, psi_int data) {
	psi_rtree_data insert_data;
	insert_data.data = data;
	memcpy(insert_data.rbox, rbox, 2*sizeof(psi_rvec));
	psi_rtree_insert_fancy(insert_data, 0, 0, rtree);
}

void psi_rtree_from_mesh(psi_rtree* rtree, psi_mesh* mesh, psi_rvec* window) {

	psi_int e, v, rtsz, elemct, tind;
	psi_int vpere = mesh->elemtype;
	psi_rvec tpos[vpere], trbox[2];

	psi_printf("Building R*-tree from mesh... ");

	// compute a good starting size and init the rtree
	rtsz = (psi_int)(0.05*mesh->nelem);
	if(rtsz < MIN_RTREE_CAP) rtsz = MIN_RTREE_CAP;
	psi_rtree_init(rtree, rtsz); // TODO: heuristic!!!

	// loop over all elements in the mesh to insert
	for(e = 0, elemct = 0; e < mesh->nelem; ++e) {
		for(v = 0; v < vpere; ++v) {
			tind = mesh->connectivity[e*vpere+v];
			tpos[v] = mesh->pos[tind];
		}
        if(!psi_aabb_periodic(tpos, trbox, window, mesh)) continue;
		psi_rtree_insert(rtree, trbox, e);
        elemct++;
	}
	psi_printf("done. Inserted %d elements.\n", elemct);

}

void psi_rtree_query_init(psi_rtree_query* qry, psi_rtree* rtree, psi_rvec* qbox) {
	qry->rtree = rtree;
	qry->rbox = qbox;
	qry->stack[0].index = rtree->root;
	qry->stack[0].level = rtree->depth;
	qry->nstack = 1;
}

psi_int psi_rtree_query_next(psi_rtree_query* qry, psi_int* ind) {

	psi_int i, c, level, index;
	psi_rtree_node* node;

	// depth-first search 
	while(qry->nstack) {
	
		// pop the stack
		--qry->nstack;
		level = qry->stack[qry->nstack].level;
		index = qry->stack[qry->nstack].index;

		// if the level is < 0, we are at a data element waiting to be returned
		// otherwise go to the next node and continue searching
		if(level < 0) {
			*ind = index;
			return 1;
		}
		else node = &qry->rtree->allnodes[index];

		// check the children
		// if this is the leaf level, push data to the buffer
		for(c = 0; c < node->nchildren; ++c) {
			for(i = 0; i < PSI_NDIM; ++i) {
				if(qry->rbox[1].xyz[i] < node->children[c].rbox[0].xyz[i]) goto next_child;
				if(qry->rbox[0].xyz[i] > node->children[c].rbox[1].xyz[i]) goto next_child;
			}
			qry->stack[qry->nstack].index = node->children[c].data;
			qry->stack[qry->nstack].level = level-1;
			qry->nstack++;

			if(qry->nstack >= STACK_MAX) {
				psi_printf("Warning! R*-tree query stack overflow!\n");
                exit(0);
			}

			next_child: continue;
		}
	}
	return 0;
}


void psi_rtree_insert_fancy(psi_rtree_data newnode, psi_int level, psi_int split, psi_rtree* rtree) {

	psi_int i, closest, c, chind, lvl, hasdata;
	psi_int sortinds[NODE_MAX+1];

	// keep track of parent history
	psi_rtree_node* pstack[rtree->depth+1];
	psi_int istack[rtree->depth+1];

	// insert algorithm
	psi_rtree_node* node;
	psi_rvec* rbox; 
	psi_rtree_node* parent;
	psi_rtree_node* splitnodes[2];
	psi_rtree_data ofdata[NODE_MAX+1];

	// for reinserting
	psi_int nreinsert, reinsertlvl;
	psi_rtree_data reinsertdata[NODE_MAX+1];
	psi_rvec boxcenter;
	psi_real axdist, cdist2s[NODE_MAX+1];

	// first, check to see if we may overflow the node space
	// reallocate if so
	if(rtree->nnodes + rtree->depth + 1 >= rtree->nodecap) {
		rtree->nodecap = (psi_int) (1.5*rtree->nodecap);
		node = (psi_rtree_node*) psi_realloc(rtree->allnodes, rtree->nodecap*sizeof(psi_rtree_node));
		if(node) rtree->allnodes = node;
		else {
			psi_printf("Tree reallocation failed.\n");
            exit(0);
			return;
		}
	}

	// descend recursively to the optimal node for inserting the new data
	// do this in-situ to build a history of parent nodes
	node = &rtree->allnodes[rtree->root];
	for(lvl = rtree->depth; lvl > level; --lvl) {
		// if the children of n are leaves, use minimum overlap criterion
		// otherwise use the minimum volume increase criterion
		if(lvl == 1) closest = psi_rtree_min_overlap(node, newnode.rbox);
		else closest = psi_rtree_min_area(node, newnode.rbox);
		pstack[lvl] = node; 
		istack[lvl] = closest;
		node = &rtree->allnodes[node->children[closest].data];
	} 

	// insert the data, splitting or reinserting as needed if overflowed
	// propagate changes upward in the tree until we reach the root
	hasdata = 1;
	nreinsert = 0;
	for(; lvl < rtree->depth; ++lvl) {
		parent = pstack[lvl+1]; // the current node's parent
		chind = istack[lvl+1]; // the index from the current node's parent to the current node
		rbox = parent->children[chind].rbox; // this node's rbox (was stored in the parent)
		if(hasdata) {
			if(node->nchildren < NODE_MAX) {
				// if the node has room, simply insert the new data
				node->children[node->nchildren++] = newnode;
				hasdata = 0;
			}	
			else {
				// else, invoke the overflow treatment
				// consolidate all overflowed children prior to reinsertion or split 
				memcpy(ofdata, node->children, NODE_MAX*sizeof(psi_rtree_data));
				ofdata[NODE_MAX] = newnode;
				if(split) {
					// we are to split this node and continue adjusting the tree
					// keep storage dense by reusing the current node
					// and sending the new one to the next loop
					splitnodes[0] = node; 
					splitnodes[1] = &rtree->allnodes[rtree->nnodes];
					psi_rtree_split(splitnodes, ofdata);
					psi_rtree_rbox(splitnodes[1], newnode.rbox);
					newnode.data = rtree->nnodes++;
					// reinsert if we cause higher levels to overflow
					split = 0;
				}
				else {
					// sort entries by their distance to the center of the node, reinsert the outliers
					for(i = 0; i < PSI_NDIM; ++i)
						boxcenter.xyz[i] = 0.5*(rbox[0].xyz[i]+rbox[1].xyz[i]);
					for(c = 0; c < NODE_MAX+1; ++c) {
						for(i = 0; i < PSI_NDIM; ++i) {
							axdist = boxcenter.xyz[i]-0.5*(node->children[c].rbox[0].xyz[i]+node->children[c].rbox[1].xyz[i]);
							cdist2s[c] += axdist*axdist;
						}
					}
					for(c = 0; c < NODE_MAX+1; ++c) sortinds[c] = c;
					psi_rtree_sort_helper_data = cdist2s; 
					qsort(sortinds, NODE_MAX+1, sizeof(psi_int), psi_rtree_sort_helper);
					// copy the closest entries back to the node
					// save the most distant ones to reinsert, in order of increasing distance
					node->nchildren = 0;
					for(c = 0; c < KEEP_ON_REINSERT; ++c)
						node->children[node->nchildren++] = ofdata[sortinds[c]];
					for(c = KEEP_ON_REINSERT; c < NODE_MAX+1; ++c) {
						reinsertdata[nreinsert++] = ofdata[sortinds[c]];
						reinsertlvl = lvl;
					}
					hasdata = 0;
				}
			}
		}
		// propagate bounding rectangles upward
		psi_rtree_rbox(node, rbox);
		node = parent;
	}
	// if there is still data when we reach the root, insert it or split
	if(hasdata) {
		if(node->nchildren < NODE_MAX) {
			node->children[node->nchildren] = newnode;
			node->nchildren++;
		}	
		else {
			// This works, please don't touch it
			memcpy(ofdata, node->children, NODE_MAX*sizeof(psi_rtree_data));
			ofdata[NODE_MAX] = newnode;
			splitnodes[0] = node; 
			splitnodes[1] = &rtree->allnodes[rtree->nnodes++];
			psi_rtree_split(splitnodes, ofdata);
			node = &rtree->allnodes[rtree->nnodes++];
			node->children[0].data = rtree->root; // the old root
			node->children[1].data = rtree->nnodes-2;
			psi_rtree_rbox(splitnodes[0], node->children[0].rbox); 
			psi_rtree_rbox(splitnodes[1], node->children[1].rbox); 
			node->nchildren = 2;
			rtree->depth++;
			rtree->root = rtree->nnodes-1;
		}
	}
	// reinsert overflowed data if we need to, specifying to split
	// on the first overflow
	for(c = 0; c < nreinsert; ++c)
		psi_rtree_insert_fancy(reinsertdata[c], reinsertlvl, 1, rtree);
}

void psi_rtree_destroy(psi_rtree* rtree) {
	rtree->root = -1;
	psi_free(rtree->allnodes);
}

void psi_rtree_init(psi_rtree* rtree, psi_int capacity) {
	rtree->nodecap = capacity; 
	rtree->allnodes = psi_malloc(capacity*sizeof(psi_rtree_node));
	if(!rtree->allnodes) {
		psi_printf("Error! Rtree allocation failed!\n");
            exit(0);
		return;
	}
	rtree->root = 0; 
	rtree->allnodes[rtree->root].nchildren = 0;
	rtree->nnodes = 1;
	rtree->depth = 0;
}

int psi_rtree_sort_helper(const void* p1, const void* p2) {
	psi_int i1 = *(int*)p1;
	psi_int i2 = *(int*)p2;
	if(psi_rtree_sort_helper_data[i1] < psi_rtree_sort_helper_data[i2]) return -1;
	if(psi_rtree_sort_helper_data[i1] == psi_rtree_sort_helper_data[i2]) return 0;
	if(psi_rtree_sort_helper_data[i1] > psi_rtree_sort_helper_data[i2]) return 1;
	return 0;
}

void psi_rtree_rbox(psi_rtree_node* node, psi_rvec* rbox) {
	// get the total bounding box of a node's children 
	psi_int i, c;
	for(i = 0; i < PSI_NDIM; ++i) {
		rbox[0].xyz[i] = 1.0e30;
		rbox[1].xyz[i] = -1.0e30;
	}
	for(c = 0; c < node->nchildren; ++c)
	for(i = 0; i < PSI_NDIM; ++i) {
		if(node->children[c].rbox[0].xyz[i] < rbox[0].xyz[i]) rbox[0].xyz[i] = node->children[c].rbox[0].xyz[i];
		if(node->children[c].rbox[1].xyz[i] > rbox[1].xyz[i]) rbox[1].xyz[i] = node->children[c].rbox[1].xyz[i];
	}
}


void psi_rtree_split(psi_rtree_node* newnodes[2], psi_rtree_data* entries) {

	psi_int lohi, splitind, group, i, j, c;
	psi_real mtmp, atmp[2];
	psi_real sortdata[NODE_MAX+1];
	psi_int sortinds[NODE_MAX+1];
	psi_rvec splitboxes[2][2]; 

	// for keeping track of the best possible split
	psi_int axis, best_axis;
	psi_real margin, overlap, area, best_margin, best_overlap, best_area;
	psi_int best_sorts[PSI_NDIM][NODE_MAX+1];
	psi_int best_split_inds[PSI_NDIM];

	psi_rtree_sort_helper_data = sortdata; 
	best_margin = 1.0e30;
	best_axis = -1;
	for(axis = 0; axis < PSI_NDIM; ++axis) {
		margin = 0.0;
		best_overlap = 1.0e30;
		best_area = 1.0e30;
		for(lohi = 0; lohi < 1; ++ lohi) {

			// sort by the low or hi box bounds
			for(c = 0; c < NODE_MAX+1; ++c) {
				sortdata[c] = entries[c].rbox[lohi].xyz[axis];
				sortinds[c] = c;
			}
			qsort(sortinds, NODE_MAX+1, sizeof(psi_int), psi_rtree_sort_helper);
	
			// divide into several distributions based on the sort
			// splitind determines the split within the sorted group
			for(splitind = NODE_MIN; splitind < (NODE_MAX-NODE_MIN+1); ++splitind) {

				// get the two bounding boxes for this split
				for(group = 0; group < 2; ++group)
					for(i = 0; i < PSI_NDIM; ++i) {
						splitboxes[group][0].xyz[i] = 1.0e30;
						splitboxes[group][1].xyz[i] = -1.0e30;
					}
				for(c = 0; c < NODE_MAX+1; ++c) {
					group = (c >= splitind);
					for(i = 0; i < PSI_NDIM; ++i) {
						if(entries[sortinds[c]].rbox[0].xyz[i] < splitboxes[group][0].xyz[i]) 
							splitboxes[group][0].xyz[i] = entries[sortinds[c]].rbox[0].xyz[i];
						if(entries[sortinds[c]].rbox[1].xyz[i] > splitboxes[group][1].xyz[i]) 
							splitboxes[group][1].xyz[i] = entries[sortinds[c]].rbox[1].xyz[i];
					}
				}

				// sum the margin-values for each box
				// TODO: This is wrong for dimensions higher than 3
				for(group = 0; group < 2; ++group)
					for(i = 0; i < PSI_NDIM; ++i) {
						mtmp = 1.0;
						for(j = 0; j < PSI_NDIM-1; ++j)
							mtmp *= splitboxes[group][1].xyz[(i+j)%PSI_NDIM]-splitboxes[group][0].xyz[(i+j)%PSI_NDIM];
						margin += 2.0*mtmp;
					}

				// get the area value
				for(group = 0; group < 2; ++group) {
					atmp[group] = 1.0;
					for(i = 0; i < PSI_NDIM; ++i)
						atmp[group] *= splitboxes[group][1].xyz[i]-splitboxes[group][0].xyz[i];
				}
				area = atmp[0] + atmp[1];

				// get the overlap between the two boxes
				overlap = 1.0;
				for(i = 0; i < PSI_NDIM; ++i) {
					if(splitboxes[0][0].xyz[i] > splitboxes[1][1].xyz[i] ||
							splitboxes[0][1].xyz[i] < splitboxes[1][0].xyz[i]) {
						overlap = 0.0;
						break;
					}
					overlap *= fmin(splitboxes[0][1].xyz[i], splitboxes[1][1].xyz[i])
								- fmax(splitboxes[0][0].xyz[i], splitboxes[1][0].xyz[i]);
				}

				if(overlap < best_overlap || (overlap == best_overlap && area < best_area)) {
					best_overlap = overlap;
					best_area = area;
					memcpy(best_sorts[axis], sortinds, (NODE_MAX+1)*sizeof(psi_int));
					best_split_inds[axis] = splitind;
				}
			}
		}
		if(margin < best_margin) {
			best_margin = margin;
			best_axis = axis;
		}
	}

	// now, we have information on the best split in hand
	// use the indices we saved to redistribute the entries
	newnodes[0]->nchildren = 0;
	newnodes[1]->nchildren = 0;
	for(c = 0; c < NODE_MAX+1; ++c) {
		group = (c >= best_split_inds[best_axis]);
		newnodes[group]->children[newnodes[group]->nchildren++] = entries[best_sorts[best_axis][c]];
	}
}

psi_int psi_rtree_min_area(psi_rtree_node* node, psi_rvec* rbox) {

	psi_int i, c, closest;
	psi_real vol, delta, ovol, minvol, mindelta;
	psi_rvec orbox[2], newbox[2];

	// query the current node's children for intersection/closeness
	// using volume of the resulting combined bounding box as a proxy for distance
	// cur_node is assumed to not be a leaf
	minvol = 1.0e30;
	mindelta = 1.0e30;
	closest = -1;
	for(c = 0; c < node->nchildren; ++c) {
		vol = 1.0;
		ovol = 1.0;
		for(i = 0; i < PSI_NDIM; ++i) {
			orbox[0].xyz[i] = node->children[c].rbox[0].xyz[i];
			orbox[1].xyz[i] = node->children[c].rbox[1].xyz[i];
			newbox[0].xyz[i] = rbox[0].xyz[i];
			newbox[1].xyz[i] = rbox[1].xyz[i];
			if(orbox[0].xyz[i] < newbox[0].xyz[i]) newbox[0].xyz[i] = orbox[0].xyz[i];
			if(orbox[1].xyz[i] > newbox[1].xyz[i]) newbox[1].xyz[i] = orbox[1].xyz[i];
			ovol *= orbox[1].xyz[i] - orbox[0].xyz[i];
			vol *= newbox[1].xyz[i] - newbox[0].xyz[i];
		}
		// minimum area enlargement criterion
		// break ties where mindelta == 0.0 using minvol
		delta = vol - ovol;
		if((delta < mindelta) || (delta == mindelta && vol < minvol)) {
			mindelta = delta;
			minvol = vol;
			closest = c;
		}
	}
	return closest;
}

psi_int psi_rtree_min_overlap(psi_rtree_node* node, psi_rvec* rbox) {

	psi_int i, c, ct, closest;
	psi_real vol, delta, ovol, minvol, mindelta;
	psi_real oldov, newov, deltaov, min, max, mindeltaov;
	psi_rvec oldbox[2], newbox[2];

	// quadratic complexity
	closest = -1;
	minvol = 1.0e30;
	mindelta = 1.0e30;
	mindeltaov = 1.0e30;
	for(c = 0; c < node->nchildren; ++c) {
		// make original and new bounding boxes, 
		// also keep track of volume and dv in case of ties
		vol = 1.0;
		ovol = 1.0;
		for(i = 0; i < PSI_NDIM; ++i) {
			oldbox[0].xyz[i] = node->children[c].rbox[0].xyz[i];
			oldbox[1].xyz[i] = node->children[c].rbox[1].xyz[i];
			newbox[0].xyz[i] = rbox[0].xyz[i];
			newbox[1].xyz[i] = rbox[1].xyz[i];
			if(oldbox[0].xyz[i] < newbox[0].xyz[i]) newbox[0].xyz[i] = oldbox[0].xyz[i];
			if(oldbox[1].xyz[i] > newbox[1].xyz[i]) newbox[1].xyz[i] = oldbox[1].xyz[i];
			ovol *= oldbox[1].xyz[i] - oldbox[0].xyz[i];
			vol *= newbox[1].xyz[i] - newbox[0].xyz[i];
		}
		delta = vol - ovol;
		// sum up the area of all overlapping other children
		deltaov = 0.0;
		for(ct = 0; ct < node->nchildren; ++ct) {
			if(c == ct) continue;
			// get the old overlap
			oldov = 1.0;
			for(i = 0; i < PSI_NDIM; ++i) {
				if(node->children[ct].rbox[0].xyz[i] > oldbox[1].xyz[i] ||
						node->children[ct].rbox[1].xyz[i] < oldbox[0].xyz[i]) {
					oldov = 0.0;
					break;
				}
				if(node->children[ct].rbox[0].xyz[i] < oldbox[0].xyz[i]) min = oldbox[0].xyz[i]; 
				else min = node->children[ct].rbox[0].xyz[i];
				if(node->children[ct].rbox[1].xyz[i] > oldbox[1].xyz[i]) max = oldbox[1].xyz[i]; 
				else max = node->children[ct].rbox[1].xyz[i];
				oldov *= max - min;
			}
			// get the new overlap
			newov = 1.0;
			for(i = 0; i < PSI_NDIM; ++i) {
				if(node->children[ct].rbox[0].xyz[i] > newbox[1].xyz[i] ||
						node->children[ct].rbox[1].xyz[i] < newbox[0].xyz[i]) {
					newov = 0.0;
					break;
				}
				if(node->children[ct].rbox[0].xyz[i] < newbox[0].xyz[i]) min = newbox[0].xyz[i]; 
				else min = node->children[ct].rbox[0].xyz[i];
				if(node->children[ct].rbox[1].xyz[i] > newbox[1].xyz[i]) max = newbox[1].xyz[i]; 
				else max = node->children[ct].rbox[1].xyz[i];
				newov *= max - min;
			}
			deltaov += newov - oldov;
		}
		// choose the node with the smalles increase in overlap, break ties with change in volume, then volume
		if((deltaov < mindeltaov) || (deltaov == mindeltaov && delta < mindelta) || (deltaov == mindeltaov && delta == mindelta && vol < minvol)) {
			mindeltaov = deltaov;
			mindelta = delta;
			minvol = vol;
			closest = c;
		}
	}
	return closest;
}

#if 0

void psi_rtree_print_recursive(FILE* stream, psi_rtree* tree, psi_rtree_node* node, psi_int lvl) {
	psi_int i, c;
	fprintf(stream, "# Node at %p, lvl = %d, %d children:\n", node, lvl, node->nchildren);
	for(c = 0; c < node->nchildren; ++c) {
		fprintf(stream, " 0x%x %d", node->children[c].data, lvl);
		for(i = 0; i < PSI_NDIM; ++i) fprintf(stream, " %f", node->children[c].rbox[0].xyz[i]);
		for(i = 0; i < PSI_NDIM; ++i) fprintf(stream, " %f", node->children[c].rbox[1].xyz[i]);
		fprintf(stream, "\n");
	}
	fprintf(stream, "# ----\n");
	if(lvl == 0) return;
	for(c = 0; c < node->nchildren; ++c)
		psi_rtree_print_recursive(stream, tree, &tree->allnodes[node->children[c].data], lvl-1);
}

void psi_rtree_print(FILE* stream, psi_rtree* tree) {
	if(!stream) stream = stdout;
	fprintf(stream, "# psi_rtree, nnodes = %d, root = %d, depth = %d\n", tree->nnodes, tree->root, tree->depth);
	psi_rtree_print_recursive(stream, tree, &tree->allnodes[tree->root], tree->depth);
}

void psi_rtree_tofile(char* fname, psi_rtree* tree) {
	FILE* f = fopen(fname, "w");
	if(!f) {
		printf("Failed to open file: %s\n", fname);
		return;
	}
	psi_rtree_print(f, tree);
	fclose(f);
}

#endif
