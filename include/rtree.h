/*
 *
 * rtree.h 
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

#ifndef _RTREE_H_
#define _RTREE_H_

#include "psi.h"
#include "mesh.h"

// heuristic values from Beckmann et al 1990
#define NODE_MAX 32
#define NODE_MIN ((4*NODE_MAX)/10)
#define KEEP_ON_REINSERT ((7*NODE_MAX)/10)
#define STACK_MAX 1024 // TODO: how big should the stack be??
#define MIN_RTREE_CAP 1024 

// a single data element (aabb and index) 
typedef struct {
	psi_rvec rbox[2];
	psi_int data;
} psi_rtree_data;

// a node containing its children
typedef struct {
	psi_short nchildren;
	psi_rtree_data children[NODE_MAX];
} psi_rtree_node;

// rtree struct containing everything we need
typedef struct {
	// tree storage
	psi_int root;
	psi_int depth;
	psi_rtree_node* allnodes;
	psi_int nodecap;
	psi_int nnodes;
} psi_rtree;

// a stack for easy query boilerplate
typedef struct {
	struct {
		psi_int index;
		psi_int level;
	} stack[STACK_MAX];
	psi_int nstack;
	psi_rvec* rbox;
	psi_rtree* rtree;
} psi_rtree_query;


void psi_rtree_from_mesh(psi_rtree* rtree, psi_mesh* mesh, psi_rvec* window); 

void psi_rtree_init(psi_rtree* rtree, psi_int capacity);

void psi_rtree_destroy(psi_rtree* rtree);

void psi_rtree_insert(psi_rtree* rtree, psi_rvec* rbox, psi_int data);

void psi_rtree_query_init(psi_rtree_query* qry, psi_rtree* rtree, psi_rvec* qbox);

psi_int psi_rtree_query_next(psi_rtree_query* qry, psi_int* ind);

void psi_rtree_print(FILE* stream, psi_rtree* tree);


#endif // _RTREE_H_
