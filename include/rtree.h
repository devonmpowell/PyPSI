#ifndef _RTREE_H_
#define _RTREE_H_

#include "psi.h"
#include "mesh.h"

// heuristic values from Beckmann et al 1990
#define NODE_MAX 32
#define NODE_MIN ((4*NODE_MAX)/10)
#define KEEP_ON_REINSERT ((7*NODE_MAX)/10)
#define STACK_MAX 512 // TODO: how big should the stack be??

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


void psi_rtree_load(psi_rtree* rtree, psi_mesh* mesh); 

void psi_rtree_init(psi_rtree* rtree, psi_int capacity);

void psi_rtree_destroy(psi_rtree* rtree);

void psi_rtree_insert(psi_rtree* rtree, psi_rvec* rbox, psi_int data);

void psi_rtree_query_init(psi_rtree_query* qry, psi_rtree* rtree, psi_rvec* qbox);

psi_int psi_rtree_query_next(psi_rtree_query* qry, psi_int* ind);

void psi_rtree_print(FILE* stream, psi_rtree* tree);


#endif // _RTREE_H_
