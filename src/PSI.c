// Python interface only
#ifdef PYMODULE 

#include <numpy/arrayobject.h>
#include "structmember.h"
#include "psi.h"
#include "grid.h"
#include "mesh.h"
#include "rtree.h"
#include "skymap.h"
#include "beamtrace.h"
#ifdef HAVE_FFTW
#include "fft.h"
#endif

#if PY_MAJOR_VERSION >= 3
#define PyInt_AsLong(x) (PyLong_AsLong((x)))
#define PyString_Check(x) (PyBytes_Check((x)))
#define PyInt_AS_LONG(x) (PyLong_AsLong((x)))
#define PyString_AsString(x) (PyUnicode_AsUTF8((x)))
#define PyInt_Check(x) (PyLong_Check((x)))
#endif

/////////////////////////////////////////////////////////////
//   Define Python type structs and conversion functions up front 
/////////////////////////////////////////////////////////////

typedef struct {
    PyObject_HEAD
    PyObject *pos, *vel, *mass;
    PyObject *connectivity;
    PyObject *boxmin, *boxmax;
} Mesh;

typedef struct {
    PyObject_HEAD
	PyObject* type; // string representation of the grid type 
	PyObject* fields; // Python dict of numpy arrays 
	PyObject* winmin, *winmax; // Python tuple of the projection window
	PyObject* n; // number of cells in each dimension
	PyObject* d; // tuple of dx in each dimension
} Grid;

typedef struct {
    PyObject_HEAD
	psi_rtree ctree; 
} RStarTree;


typedef struct {
    PyObject_HEAD
	PyObject* type; // string representation of the grid type 
	PyObject* params; // tuple of floats -- parameters for analytic metrics 
	PyObject* n; // number of cells in each dimension
	PyObject* box; // box size 
	PyObject* phi; // numpy array - potential 
	PyObject* gradphi; // numpy array - gradient of the potential 
} Metric;

static void PSI_Mesh2mesh(Mesh* mesh, psi_mesh* cmesh) {

	// parse a simple C struct from a Python object 
	psi_int ax;
	PyObject* seq;
	npy_intp* npdims;

	memset(cmesh, 0, sizeof(psi_mesh));
   
	// get the dimensionality based on shape of pos 
	// elemtype based on shape of connectivity 
	npdims = PyArray_SHAPE((PyArrayObject*)mesh->pos);
	cmesh->dim = npdims[1];
	cmesh->npart = npdims[0];
	npdims = PyArray_SHAPE((PyArrayObject*)mesh->connectivity);
	cmesh->elemtype = npdims[1];
	cmesh->nelem = npdims[0];

	// get periodicity based on whether boxmin, boxmax exist
	if(PySequence_Check(mesh->boxmin) && PySequence_Check(mesh->boxmax)){
		cmesh->periodic = 1;
		seq = PySequence_Fast(mesh->boxmin, "expected a box tuple");
		for(ax = 0; ax < cmesh->dim; ++ax)
			cmesh->box[0].xyz[ax] = PyFloat_AS_DOUBLE(PySequence_Fast_GET_ITEM(seq, ax));
		seq = PySequence_Fast(mesh->boxmax, "expected a box tuple");
		for(ax = 0; ax < cmesh->dim; ++ax)
			cmesh->box[1].xyz[ax] = PyFloat_AS_DOUBLE(PySequence_Fast_GET_ITEM(seq, ax));
	}
	else {
		cmesh->periodic = 0;
	}

	// set array pointers
	// TODO: mass...
	// TODO: increfs!!!
	cmesh->pos = PyArray_DATA((PyArrayObject*)mesh->pos);
	cmesh->vel = PyArray_DATA((PyArrayObject*)mesh->vel);
	cmesh->mass = PyArray_DATA((PyArrayObject*)mesh->mass);
	cmesh->connectivity = PyArray_DATA((PyArrayObject*)mesh->connectivity);

}

static void PSI_Grid2grid(Grid* grid, psi_grid* cgrid) {

	// parse a simple C struct from a Python object 
	psi_int ax;
	PyObject* seq;
	char* cstr;

	// clear everything first
	memset(cgrid, 0, sizeof(psi_grid));
   
	// get the enum grid type
	cstr = PyString_AsString(grid->type);
	if(strcmp(cstr, "cart") == 0)
		cgrid->type = PSI_GRID_CART;
	else if(strcmp(cstr, "hpring") == 0)
		cgrid->type = PSI_GRID_HPRING;
	else
		cgrid->type = -1;

	// get the grid dimensions, dx, etc
	cgrid->dim = PySequence_Length(grid->n);
	seq = PySequence_Fast(grid->n, "expected a grid dim tuple");
	for(ax = 0; ax < cgrid->dim; ++ax)
		cgrid->n.ijk[ax] = PyInt_AS_LONG(PySequence_Fast_GET_ITEM(seq, ax));
	
	
	if(PySequence_Check(grid->winmin)) {
		seq = PySequence_Fast(grid->winmin, "expected a grid window tuple");
		for(ax = 0; ax < cgrid->dim; ++ax)
			cgrid->window[0].xyz[ax] = PyFloat_AS_DOUBLE(PySequence_Fast_GET_ITEM(seq, ax));
		seq = PySequence_Fast(grid->winmax, "expected a grid window tuple");
		for(ax = 0; ax < cgrid->dim; ++ax)
			cgrid->window[1].xyz[ax] = PyFloat_AS_DOUBLE(PySequence_Fast_GET_ITEM(seq, ax));
	}
	if(PySequence_Check(grid->d)) {
		seq = PySequence_Fast(grid->d, "expected a grid dx tuple");
		for(ax = 0; ax < cgrid->dim; ++ax)
			cgrid->d.xyz[ax] = PyFloat_AS_DOUBLE(PySequence_Fast_GET_ITEM(seq, ax));
	}
	else if(PyFloat_Check(grid->d)) {
		cgrid->d.x = PyFloat_AsDouble(grid->d);
	}

	if(PyDict_Check(grid->fields)) {
	
		PyObject *key, *value;
		Py_ssize_t pos = 0;
		while (PyDict_Next(grid->fields, &pos, &key, &value)) {
			char* cstring = PyString_AsString(key);
			if(strcmp(cstring, "m") == 0) {
				cgrid->fields[PSI_GRID_M] = (psi_real*)PyArray_DATA((PyArrayObject*)PyDict_GetItemString(grid->fields, "m"));
			}
			else if(strcmp(cstring, "x") == 0) {
				cgrid->fields[PSI_GRID_X] = (psi_real*)PyArray_DATA((PyArrayObject*)PyDict_GetItemString(grid->fields, "x"));
			}
			else if(strcmp(cstring, "v") == 0) {
				cgrid->fields[PSI_GRID_V] = (psi_real*)PyArray_DATA((PyArrayObject*)PyDict_GetItemString(grid->fields, "v"));
			}
			else if(strcmp(cstring, "xx") == 0) {
				cgrid->fields[PSI_GRID_XX] = (psi_real*)PyArray_DATA((PyArrayObject*)PyDict_GetItemString(grid->fields, "xx"));
			}
			else if(strcmp(cstring, "xv") == 0 || strcmp(cstring, "vx") == 0) {
				cgrid->fields[PSI_GRID_XV] = (psi_real*)PyArray_DATA((PyArrayObject*)PyDict_GetItemString(grid->fields, "xv"));
			}
			else if(strcmp(cstring, "vv") == 0) {
				cgrid->fields[PSI_GRID_VV] = (psi_real*)PyArray_DATA((PyArrayObject*)PyDict_GetItemString(grid->fields, "vv"));
			}
		}
	}
}

static void PSI_Metric2metric(Metric* metric, psi_metric* cmetric) {


	// parse a simple C struct from a Python object 
	psi_int ax;
	PyObject* seq;
	char* cstr;


	// clear everything first
	memset(cmetric, 0, sizeof(psi_metric));
	cmetric->snapnum = -1;
   
	// get the enum metric type
	cstr = PyString_AsString(metric->type);
	printf("Metric type = %s\n", cstr);

	if(strcmp(cstr, "minkowski") == 0) {
		printf("Chose the Minkowski metric!\n");
		cmetric->type = PSI_METRIC_MINKOWSKI; 
	}
	else if(strcmp(cstr, "flrw") == 0) {
		printf("Chose the FLRW metric!\n");
		cmetric->type = PSI_METRIC_FLRW; 
	}
	else if(strcmp(cstr, "kerr") == 0) {
		printf("Chose the KERR metric!\n");
		cmetric->type = PSI_METRIC_KERR; 
	}
	else {
		psi_printf("Invalid metric\n");
		return NULL;
	}
}



/////////////////////////////////////////////////////////////
//    Mesh 
/////////////////////////////////////////////////////////////

static int Mesh_init(Mesh *self, PyObject *args, PyObject *kwds) {

	psi_int ax, ndim;
	psi_mesh cmesh;
	psi_dvec nside;
	PyObject *posar, *velar, *massar, *box, *n, *connar;
	char* cloader, *cfile;
	npy_intp npdims[6];
	npy_intp *dtmp; 
	static char *kwlist[] = {"loader", "filename", "posar", "velar", "massar", "connar", "box", "n", NULL};


	// TODO: init all arrays to NULL
	box = NULL;

	// parse arguments differently depending on what the grid type is
	// certain grid types must correspond to certain arg patterns
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "s|sOOOOOO", kwlist, &cloader, &cfile, &posar, &velar, &massar, &connar, &box, &n))
			return -1;

	// fill in all grid information as Python tuples
	if(strcmp(cloader, "array") == 0) {

		// TODO: check dtype!!!!

		if(!PyArray_Check(posar))
			return -1;

		// TODO: check for valid box in a dedicated function
		// TODO: more rigorous check
		if(PySequence_Check(box)) {
			self->boxmin = PySequence_GetItem(box, 0); 
			self->boxmax = PySequence_GetItem(box, 1); 
		}
		else {
			self->boxmin = Py_None; 
			self->boxmax = Py_None;
		}

		//if()
		self->connectivity = connar;
		Py_INCREF(connar);


		self->pos = posar;
		Py_INCREF(posar);
		self->vel = velar;
		Py_INCREF(velar);
		self->mass = massar;
		Py_INCREF(massar);

	}
	else if(strcmp(cloader, "block") == 0) {

		// TODO: check dtype!!!!

		setbuf(stdout, NULL);

		if(!PyArray_Check(posar))
			return -1;

		ndim = PyArray_NDIM(posar);
		dtmp = PyArray_DIMS(posar);

		psi_printf("Pos array dimensions are %d %d %d %d, ndim = %d\n", dtmp[0], dtmp[1], dtmp[2], dtmp[3], ndim);

		self->boxmin = PySequence_GetItem(box, 0); 
		self->boxmax = PySequence_GetItem(box, 1); 

		if(!PySequence_Check(n)) {
			psi_printf("Must provide particle block dimensions!\n");
			return -1;
		}
		for(ax = 0; ax < 3; ++ax)
			nside.ijk[ax] = PyInt_AsLong(PySequence_GetItem(n, ax));

		printf("Loaded nside = %d %d %d\n", nside.i, nside.j, nside.k);



		self->pos = posar;
		Py_INCREF(posar);

		// TODO: set better error handling here
		self->vel = PyArray_SimpleNew(ndim, dtmp, NPY_DOUBLE);
		self->mass = PyArray_SimpleNew(ndim-1, dtmp, NPY_DOUBLE);

		npdims[0] = nside.i*nside.j*nside.k;
		if(npdims[0] != dtmp[0]) {
			psi_printf("Number of particles does not match block dimensions!\n");
			return -1;
		}
		npdims[1] = psi_verts_per_elem(PSI_MESH_LINEAR); 
		self->connectivity = PyArray_SimpleNew(2, npdims, NPY_INT32);
		psi_int* cptr = (psi_int*) PyArray_DATA(self->connectivity); 

		psi_printf("Using the array loader.\n");

		// now build the mesh connectivity
		// trilinear elements naturally
		psi_int i, j, k, ii, jj, kk, locind, vertind;
		psi_int elemind;
		for(i = 0; i < nside.i; ++i)
		for(j = 0; j < nside.j; ++j)
		for(k = 0; k < nside.k; ++k) {
			elemind = nside.j*nside.k*i + nside.k*j + k;
			for(ii = 0; ii < 2; ++ii)
			for(jj = 0; jj < 2; ++jj)
			for(kk = 0; kk < 2; ++kk) {
				locind = 4*ii + 2*jj + kk;
				vertind = nside.j*nside.k*((i+ii)%nside.i) 
					+ nside.k*((j+jj)%nside.j) + ((k+kk)%nside.k);
				cptr[8*elemind+locind] = vertind;
			}
		}

	}
	else if(strcmp(cloader, "gadget2") == 0 || 
			strcmp(cloader, "gevolution") == 0) {

		if(strcmp(cloader, "gadget2") == 0) {
			psi_printf("Using the Gadget2 loader.\n");
		} 
		else if(strcmp(cloader, "gevolution") == 0) {
			psi_printf("Using the Gevolution loader.\n");
		} 

		// peek at the file first to make sure it's there
		// and to get header information
		if(!peek_gadget2(&cmesh, cfile))
			return -1;

		// wrap the cmesh data in a cubical Numpy array
		npdims[0] = cmesh.npart;
		npdims[1] = cmesh.dim;
		self->pos = PyArray_SimpleNew(2, npdims, NPY_DOUBLE);
		self->vel = PyArray_SimpleNew(2, npdims, NPY_DOUBLE);
		self->mass = PyArray_SimpleNew(1, npdims, NPY_DOUBLE);
		npdims[0] = cmesh.nelem;
		npdims[1] = psi_verts_per_elem(cmesh.elemtype); 
		self->connectivity = PyArray_SimpleNew(2, npdims, NPY_INT32);
		if(!self->pos || !self->vel || !self->connectivity) {
			psi_printf("new array fail!\n");
			return -1;
		} 
		if(cmesh.periodic) {
			self->boxmin = Py_BuildValue("(ddd)", cmesh.box[0].x, cmesh.box[0].y, cmesh.box[0].z);
			self->boxmax = Py_BuildValue("(ddd)", cmesh.box[1].x, cmesh.box[1].y, cmesh.box[1].z);
		}
		else {
			self->boxmin = Py_None;
			self->boxmax = Py_None;
		}

		cmesh.pos = PyArray_DATA((PyArrayObject*)self->pos);
		cmesh.vel = PyArray_DATA((PyArrayObject*)self->vel);
		cmesh.mass = PyArray_DATA((PyArrayObject*)self->mass);
		cmesh.connectivity = PyArray_DATA((PyArrayObject*)self->connectivity);

		// load the data into the numpy buffers
		if(strcmp(cloader, "gadget2") == 0) {
			if(!load_gadget2(&cmesh, cfile))
				return -1;
		} 
		else if(strcmp(cloader, "gevolution") == 0) {
			if(!load_gevolution(&cmesh, cfile))
				return -1;
		} 

	}
	else {
		PyErr_SetString(PyExc_ValueError, "Invalid mesh loader.");
		return -1;
	} 

    return 0;
}

static PyObject* Mesh_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    Mesh *self;
    self = (Mesh*)type->tp_alloc(type, 0);
    return (PyObject*)self;
}


static void Mesh_dealloc(Mesh* self) {
    Py_XDECREF(self->pos);
    Py_XDECREF(self->vel);
    Py_XDECREF(self->mass);
    Py_XDECREF(self->connectivity);
    Py_XDECREF(self->boxmin);
    Py_XDECREF(self->boxmax);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyMemberDef Mesh_members[] = {
    {"pos", T_OBJECT_EX, offsetof(Mesh, pos), 0, "pos"},
    {"vel", T_OBJECT_EX, offsetof(Mesh, vel), 0, "vel"},
    {"mass", T_OBJECT_EX, offsetof(Mesh, mass), 0, "mass"},
    {"connectivity", T_OBJECT_EX, offsetof(Mesh, connectivity), 0, "connectivity"},
    {"boxmin", T_OBJECT_EX, offsetof(Mesh, boxmin), 0, "boxmin"},
    {"boxmax", T_OBJECT_EX, offsetof(Mesh, boxmax), 0, "boxmax"},
    {NULL}  /* Sentinel */
};

static PyMethodDef Mesh_methods[] = {
    {NULL}  /* Sentinel */
};

static PyTypeObject MeshType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "PSI.Mesh",             /* tp_name */
    sizeof(Mesh),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)Mesh_dealloc, /* tp_dealloc */
    0,                         /* tp_print */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_compare */
    0,                         /* tp_repr */
    0,                         /* tp_as_number */
    0,                         /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash */
    0,                         /* tp_call */
    0,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT |
        Py_TPFLAGS_BASETYPE,   /* tp_flags */
    "Mesh",           /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    Mesh_methods,             /* tp_methods */
    Mesh_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)Mesh_init,      /* tp_init */
    0,                         /* tp_alloc */
    Mesh_new,                 /* tp_new */
};

/////////////////////////////////////////////////////////////
//     RStarTree
/////////////////////////////////////////////////////////////


static PyObject* RStarTree_query(RStarTree *self, PyObject *args, PyObject *kwds) {

	psi_int e;
	PyObject* box;
	psi_rvec qbox[2];
	npy_intp npdims[3];
	static char *kwlist[] = {"box", NULL};

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O", kwlist, &box))
		return NULL;

	// parse the query box, should be a tuple
	if(!PyArg_ParseTuple(box, "(ddd)(ddd)", &qbox[0].x, &qbox[0].y, &qbox[0].z, 
				&qbox[1].x, &qbox[1].y, &qbox[1].z))
		return NULL;

	// set up space to return all of the queried elements
	psi_int capacity = 1024; 
	psi_int nqry = 0;
	psi_int* inds_out = (psi_int*) psi_malloc(capacity*sizeof(psi_int));

	// query the rtree
	psi_rtree_query qry;
	psi_rtree_query_init(&qry, &self->ctree, qbox);
	while(psi_rtree_query_next(&qry, &e)) {

		// grow the buffer if needed
		if(nqry >= capacity) {
			capacity *= 2; 
			inds_out = (psi_int*) psi_realloc(inds_out, capacity*sizeof(psi_int));
		}
		
		// save the result in the return array
		inds_out[nqry++] = e;
	}

	// give the buffer to a Numpy array and return
	npdims[0] = nqry;
	PyObject* pyinds = PyArray_SimpleNewFromData(1, npdims, NPY_INT32, inds_out);
	PyArray_ENABLEFLAGS(pyinds, NPY_ARRAY_OWNDATA);
	return Py_BuildValue("O", pyinds);

}


static int RStarTree_init(RStarTree *self, PyObject *args, PyObject *kwds) {

	psi_int ax, dim, order, nside, npix, init_cap;
    psi_rvec tbox[2];
	npy_intp npdims[3];
	psi_mesh cmesh;
	PyObject* mesh;
	static char *kwlist[] = {"mesh", "initial_capacity", NULL};

	// parse arguments differently depending on what the grid type is
	// certain grid types must correspond to certain arg patterns
	init_cap = 1024;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|i", kwlist, &mesh, &init_cap))
			return -1;

	psi_printf("Initial capacity = %d\n", init_cap);

	PSI_Mesh2mesh(mesh, &cmesh);

	// TODO: Better user-provided mesh handling 
    for(ax = 0; ax < 3; ++ax) {
        tbox[0].xyz[ax] = -1.0e10;
        tbox[1].xyz[ax] = +1.0e10;
    }

	// TODO: pass inital capacity
	psi_rtree_from_mesh(&self->ctree, &cmesh, &tbox);

    return 0;
}

static PyObject* RStarTree_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    RStarTree *self;
    self = (RStarTree*)type->tp_alloc(type, 0);
    return (PyObject*)self;
}

static void RStarTree_dealloc(RStarTree* self) {
	psi_rtree_destroy(&self->ctree);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyMemberDef RStarTree_members[] = {
    {NULL}  /* Sentinel */
};

#if PY_MAJOR_VERSION >= 3
static PyMethodDef RStarTree_methods[] = {
	{"query", (PyCFunction)RStarTree_query, METH_KEYWORDS | METH_VARARGS,
	 "Query the RStarTree."},
    {NULL}  /* Sentinel */
};
#else
static PyMethodDef RStarTree_methods[] = {
	{"query", (PyCFunction)RStarTree_query, METH_KEYWORDS,
	 "Query the RStarTree."},
    {NULL}  /* Sentinel */
};
#endif

static PyTypeObject RStarTreeType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "PSI.RStarTree",             /* tp_name */
    sizeof(RStarTree),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)RStarTree_dealloc, /* tp_dealloc */
    0,                         /* tp_print */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_compare */
    0,                         /* tp_repr */
    0,                         /* tp_as_number */
    0,                         /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash */
    0,                         /* tp_call */
    0,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT |
        Py_TPFLAGS_BASETYPE,   /* tp_flags */
    "Basic grid",           /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    RStarTree_methods,             /* tp_methods */
    RStarTree_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)RStarTree_init,      /* tp_init */
    0,                         /* tp_alloc */
    RStarTree_new,                 /* tp_new */
};




/////////////////////////////////////////////////////////////
//     Grid
/////////////////////////////////////////////////////////////

static PyObject* Grid_getCellGeometry(Grid *self, PyObject *args, PyObject *kwds) {

	psi_rvec boundary[64], center;
	psi_real vol;
	psi_int ax, slen, v, i, allcells, bstep;
	psi_dvec grind;
	psi_grid cgrid;
	PyObject *seq=NULL, *pyind=NULL;
	PyObject *bverts=NULL, *cverts=NULL, *pyvol=NULL;
	npy_intp npdims[3];
	double *cvdata, *bvdata, *pvdata;
	static char *kwlist[] = {"cell", "bstep", NULL};

	setbuf(stdout, NULL);

	// parse args, just the PyObjetc containing the indices
	// see if it is a sequence or None to indicate all cells
	bstep = 1;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|i", kwlist, &pyind, &bstep))
		return NULL;
	allcells = (pyind == Py_None)? 1 : 0; 
	if(!allcells && !PySequence_Check(pyind)) 
		return NULL;


	// parse the grid to the simple C struct
	// make the return tuple based on the number of input elements
	// process for different indexing protocols
	PSI_Grid2grid(self, &cgrid);
	switch(cgrid.type) {
		
		case PSI_GRID_CART:

			// TODO: implement this
			slen = 16;
		
			// create the output array in the correct shape (flat for hp)
			// get direct pointers to the data
			npdims[0] = slen;
			npdims[1] = 3;
			cverts = PyArray_SimpleNew(2, npdims, NPY_DOUBLE);
			npdims[1] = 4;
			npdims[2] = 3;
			bverts = PyArray_SimpleNew(3, npdims, NPY_DOUBLE);
			pyvol = PyArray_SimpleNew(1, npdims, NPY_DOUBLE);
			break;

		case PSI_GRID_HPRING:

			// do different things if we are given a list of pixels
			// or all pixels
			// hpring is indexed by an integer
			if(allcells) {
				// get npix directly, it is stored in grid.n
				slen = cgrid.n.k; 
			}
			else {
				// check that we have a sequence of integers
				// get its length
				seq = PySequence_Fast(pyind, "not a sequence");
				slen = PySequence_Fast_GET_SIZE(seq);
				if(slen <= 0 || !PyInt_Check(PySequence_Fast_GET_ITEM(seq, 0))) 
					return NULL;
			}
	
			// create the output array in the correct shape (flat for hp)
			// get direct pointers to the data
			// iterate over all integer indices in the sequence
			// or all cells if cells=None
			npdims[0] = slen;
			npdims[1] = 3;
			cverts = PyArray_SimpleNew(2, npdims, NPY_DOUBLE);
			npdims[1] = 4*bstep;
			npdims[2] = 3;
			bverts = PyArray_SimpleNew(3, npdims, NPY_DOUBLE);
			pyvol = PyArray_SimpleNew(1, npdims, NPY_DOUBLE);
	
			cvdata = PyArray_DATA((PyArrayObject*)cverts);
			bvdata = PyArray_DATA((PyArrayObject*)bverts);
			pvdata = PyArray_DATA((PyArrayObject*)pyvol);
			for(i = 0; i < slen; ++i) {
				if(allcells) grind.i = i;
				else grind.i = PyInt_AsLong(PySequence_Fast_GET_ITEM(seq, i));
	
				// retrieve the grid geometry from the C lib
				// set the array elements 
				if(!psi_grid_get_cell_geometry(&cgrid, grind, bstep, boundary, &center, &vol)) {
					printf("bad geometry\n");
					return NULL;
				} 
				for(ax = 0; ax < 3; ++ax)
					cvdata[3*i+ax] = center.xyz[ax];
				for(v = 0; v < 4*bstep; ++v)
					for(ax = 0; ax < 3; ++ax)
						bvdata[12*bstep*i+3*v+ax] = boundary[v].xyz[ax];
				pvdata[i] = vol; 
			}

			break;


		default:
			psi_printf("Bad grid type\n");
			return NULL;
	}
	
	return Py_BuildValue("NNN", PyArray_Return((PyArrayObject*)cverts), 
			PyArray_Return((PyArrayObject*)bverts), PyArray_Return((PyArrayObject*)pyvol)); 
}


static int Grid_init(Grid *self, PyObject *args, PyObject *kwds) {

	psi_int ax, dim, order, nside, npix;
	psi_grid cgrid;
	char* type;
	npy_intp npdims[3];
	PyObject* pytype;
	PyObject* fields;
	static char *kwlist[] = {"type", "window", "n", "fields", NULL};
	static char *kwlist_hpring[] = {"type", "n", NULL};

	// default is a 32^3 unit box 
	// or healpix with nside=32
	dim = 0;
	memset(npdims, 0, sizeof(npdims));
	for(ax = 0; ax < 3; ++ax) {
		cgrid.window[0].xyz[ax] = 0.0;
		cgrid.window[1].xyz[ax] = 1.0;
		cgrid.n.ijk[ax] = 32;
	}

	// peek at the grid type before parsing args
	// parse arguments differently depending on what the grid type is
	// certain grid types must correspond to certain arg patterns
	type = PyString_AsString(PyDict_GetItemString(kwds, "type"));
	fields = NULL;
	if(strcmp(type, "cart") == 0 && 
			PyArg_ParseTupleAndKeywords(args, kwds, "S|((ddd)(ddd))(iii)O", kwlist, // for cart 
			&pytype, &cgrid.window[0].x, &cgrid.window[0].y, &cgrid.window[0].z, 
			&cgrid.window[1].x, &cgrid.window[1].y, &cgrid.window[1].z, &cgrid.n.i, &cgrid.n.j, &cgrid.n.k, &fields)) {

		// fill in all grid information as Python tuples
		self->type = pytype; 
		self->winmin = Py_BuildValue("(ddd)", cgrid.window[0].x, cgrid.window[0].y, cgrid.window[0].z); 
		self->winmax = Py_BuildValue("(ddd)", cgrid.window[1].x, cgrid.window[1].y, cgrid.window[1].z); 
		self->n = Py_BuildValue("(iii)", cgrid.n.i, cgrid.n.j, cgrid.n.k); 
		self->d = Py_BuildValue("(ddd)", (cgrid.window[1].x-cgrid.window[0].x)/cgrid.n.i, (cgrid.window[1].y-cgrid.window[0].y)/cgrid.n.j, (cgrid.window[1].z-cgrid.window[0].z)/cgrid.n.k); 
		npdims[0] = cgrid.n.i;
		npdims[1] = cgrid.n.j;
		npdims[2] = cgrid.n.k;
		dim = 3;
		Py_XDECREF(pytype);

		// make the fields dict
		// now allocate numpy storage
		self->fields = PyDict_New();
		PyDict_SetItemString(self->fields, "m", PyArray_ZEROS(dim, npdims, NPY_DOUBLE, 0));
		if(fields != NULL && fields != Py_None) {
		    PyObject* seq, *item;
		    int i, len;
		    seq = PySequence_Fast(fields, "expected a sequence");
		    len = PySequence_Size(fields);
		    for (i = 0; i < len; i++) {
		        item = PySequence_Fast_GET_ITEM(seq, i);
				if(PyString_Check(item)) {
					char* cstring = PyString_AsString(item);
					if(strcmp(cstring, "m") == 0) {;}
					else if(strcmp(cstring, "x") == 0) {
						npdims[3] = 3; 
						dim = 4;
						PyDict_SetItemString(self->fields, "x", PyArray_ZEROS(dim, npdims, NPY_DOUBLE, 0));
					}
					else if(strcmp(cstring, "v") == 0) {
						npdims[3] = 3; 
						dim = 4;
						PyDict_SetItemString(self->fields, "v", PyArray_ZEROS(dim, npdims, NPY_DOUBLE, 0));
					}
					else if(strcmp(cstring, "xx") == 0) {
						npdims[3] = 3; 
						npdims[4] = 3; 
						dim = 5;
						PyDict_SetItemString(self->fields, "xx", PyArray_ZEROS(dim, npdims, NPY_DOUBLE, 0));
					}
					else if(strcmp(cstring, "xv") == 0 || strcmp(cstring, "vx") == 0) {
						npdims[3] = 3; 
						npdims[4] = 3; 
						dim = 5;
						PyDict_SetItemString(self->fields, "xv", PyArray_ZEROS(dim, npdims, NPY_DOUBLE, 0));
					}
					else if(strcmp(cstring, "vv") == 0) {
						npdims[3] = 3; 
						npdims[4] = 3; 
						dim = 5;
						PyDict_SetItemString(self->fields, "vv", PyArray_ZEROS(dim, npdims, NPY_DOUBLE, 0));
					}
					else {
						psi_printf("Unrecognized field identifier %s\n", cstring);
					}
				}
				else {
					psi_printf("Fields must be a string!");
				}
			}
		    Py_DECREF(seq);
		
		}
	}
	else if(strcmp(type, "hpring") == 0 && 
			PyArg_ParseTupleAndKeywords(args, kwds, "S|i", kwlist_hpring, &pytype, &cgrid.n.i)) {

		// for Healpix store n as (order, nside, npix)
		self->type = pytype; 
		order = floor(log2(cgrid.n.i+0.5));
		nside = (1<<order);
		npix = 12*nside*nside;
		self->n = Py_BuildValue("(iii)", order, nside, npix); 
		self->d = Py_BuildValue("d", FOUR_PI/npix); // the nominal pixel area 
		self->winmin = Py_None;
		self->winmax = Py_None;
		npdims[0] = npix;
		dim = 1;
		Py_XDECREF(pytype);

		// make the fields dict
		// now allocate numpy storage
		self->fields = PyDict_New();
		PyDict_SetItemString(self->fields, "m", PyArray_ZEROS(dim, npdims, NPY_DOUBLE, 0));

	}
	else {
		PyErr_SetString(PyExc_ValueError, "Invalid grid type.");
		return -1;
	}

    return 0;
}

static PyObject* Grid_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    Grid *self;
    self = (Grid*)type->tp_alloc(type, 0);
    return (PyObject *)self;
}

static void Grid_dealloc(Grid* self) {
    Py_XDECREF(self->fields);
    Py_XDECREF(self->winmin);
    Py_XDECREF(self->winmax);
    Py_XDECREF(self->d);
    Py_XDECREF(self->n);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyMemberDef Grid_members[] = {
    {"type", T_OBJECT, offsetof(Grid, type), 0,
     "type"},
    {"fields", T_OBJECT, offsetof(Grid, fields), 0,
     "fields dict"},
    {"winmin", T_OBJECT, offsetof(Grid, winmin), 0,
     "projection window min"},
    {"winmax", T_OBJECT, offsetof(Grid, winmax), 0,
     "projection window max"},
    {"n", T_OBJECT, offsetof(Grid, n), 0,
     "grid size"},
    {"d", T_OBJECT, offsetof(Grid, d), 0,
     "grid dx"},
    {NULL}  /* Sentinel */
};

#if PY_MAJOR_VERSION >= 3
static PyMethodDef Grid_methods[] = {
	{"getCellGeometry", (PyCFunction)Grid_getCellGeometry, METH_KEYWORDS | METH_VARARGS,
	 "Get the vertices for the given cell."},
    {NULL}  /* Sentinel */
};
#else
static PyMethodDef Grid_methods[] = {
	{"getCellGeometry", (PyCFunction)Grid_getCellGeometry, METH_KEYWORDS,
	 "Get the vertices for the given cell."},
    {NULL}  /* Sentinel */
};
#endif

static PyTypeObject GridType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "PSI.Grid",             /* tp_name */
    sizeof(Grid),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)Grid_dealloc, /* tp_dealloc */
    0,                         /* tp_print */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_compare */
    0,                         /* tp_repr */
    0,                         /* tp_as_number */
    0,                         /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash */
    0,                         /* tp_call */
    0,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT |
        Py_TPFLAGS_BASETYPE,   /* tp_flags */
    "Basic grid",           /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    Grid_methods,             /* tp_methods */
    Grid_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)Grid_init,      /* tp_init */
    0,                         /* tp_alloc */
    Grid_new,                 /* tp_new */
};



/////////////////////////////////////////////////////////////
//     PSI module functions 
/////////////////////////////////////////////////////////////

static PyObject *PSI_skymap(PyObject *self, PyObject *args, PyObject* kwds) {

	psi_grid cgrid;
	psi_mesh cmesh;
	psi_int bstep, mode;
	Mesh* mesh;
	Grid* grid;	
	//npy_intp* npdims;
	npy_intp nverts;
	static char *kwlist[] = {"grid", "mesh", "bstep", "mode", NULL};

	// parse PyArgs
	bstep = 1;
	mode = PSI_SKYMAP_RHO_LINEAR;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OO|ii", kwlist, &grid, &mesh, &bstep, &mode))
		return NULL;

	// make C structs, check them
	PSI_Grid2grid(grid, &cgrid);
	PSI_Mesh2mesh(mesh, &cmesh);
	if(cgrid.type != PSI_GRID_HPRING || cgrid.dim != 3) 
		return NULL;
	if(cmesh.dim != cgrid.dim)
		return NULL;

	// run the skymap
	printf("----Making hpring skymap, mode = %d-----\n", mode);
	psi_skymap(&cgrid, &cmesh, bstep, mode);
	Py_RETURN_NONE;
}



static PyObject *PSI_beamtrace(PyObject *self, PyObject *args, PyObject* kwds) {

	psi_grid cgrid;
	psi_mesh cmesh;
	psi_metric cmetric;
	psi_int bstep, mode;
	Metric* metric;
	Grid* grid;	
	psi_rvec obspos, obsvel;
	psi_dvec gdim;
	npy_intp npdims[4];
	npy_intp nverts;
	static char *kwlist[] = {"grid", "metric", "obspos", "obsvel", NULL};

	setbuf(stdout, NULL);

	printf("---- GR beamtracing... -----\n");

	bstep = 1;
	mode = PSI_SKYMAP_RHO_LINEAR;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OO|(ddd)(ddd)", kwlist, &grid, &metric, &obspos.x, &obspos.y, &obspos.z, &obsvel.x, &obsvel.y, &obsvel.z))
		return NULL;

	PSI_Grid2grid(grid, &cgrid);
	if(cgrid.type != PSI_GRID_HPRING || cgrid.dim != 3) 
		return NULL;

	printf("metric2metric\n");

	PSI_Metric2metric(metric, &cmetric);

	printf("metric2metric done\n");

	printf("grid.nside = %d, metric type = %d, obspos = %f %f %f\n", cgrid.n.j, cmetric.type, obspos.x, obspos.y, obspos.z);




	npdims[0] = cgrid.n.k;
	npdims[1] = 100024; 
	npdims[2] = 4; 
	gdim.i = npdims[0];
	gdim.j = npdims[1];
	gdim.k = npdims[2];
	PyObject* rayinfo = PyArray_SimpleNew(3, npdims, NPY_DOUBLE);
	psi_real* infar = PyArray_DATA((PyArrayObject*)rayinfo);

	psi_beamtrace(&cgrid, &cmesh, bstep, &cmetric, infar, gdim);

   /* Do your stuff here. */
   //Py_RETURN_NONE;
   return rayinfo;
}




static PyObject *PSI_voxels(PyObject *self, PyObject *args, PyObject* kwds) {

	psi_int mode, maxlvl;
	psi_real reftol;
	char* modestr;
	psi_grid cgrid;
	psi_mesh cmesh;
	Mesh* mesh;
	Grid* grid;	
	npy_intp nverts;
	static char *kwlist[] = {"grid", "mesh", "mode", "refine_tolerance", "refine_max_lvl", NULL};

	// defaults
	modestr = "density";
	reftol = 1.0;
	maxlvl = 0;

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OO|sdi", kwlist, &grid, &mesh, &modestr, &reftol, &maxlvl))
		return NULL;

	// extract C pointers and such
	PSI_Grid2grid(grid, &cgrid);
	PSI_Mesh2mesh(mesh, &cmesh);
	if(cgrid.type != PSI_GRID_CART) 
		return NULL;
	if(cmesh.dim != cgrid.dim)
		return NULL;

	// parse the deposit mode
	if(strcmp(modestr, "density") == 0)
		mode = PSI_MODE_DENSITY;
	else if(strcmp(modestr, "annihilation") == 0)
		mode = PSI_MODE_ANNIHILATION;
	else {
		PyErr_SetString(PyExc_ValueError, "Invalid mode for PSI.voxels()");
		return NULL;
	}

	// call the C function
	// TODO: allow passing of a user-made r* tree
	psi_voxels(&cgrid, &cmesh, NULL, mode, reftol, maxlvl);

	Py_RETURN_NONE;
}


static PyObject *PSI_VDF(PyObject *self, PyObject *args, PyObject* kwds) {

	psi_int maxlvl;
	psi_real reftol;
	psi_rvec samppos;
	psi_mesh cmesh;
	psi_rtree* ctree;
	Mesh* mesh;
	PyObject* sPos;
	RStarTree* rst;
	npy_intp nverts;
	static char *kwlist[] = {"mesh", "sample_pos", "tree", "refine_tolerance", "refine_max_lvl", NULL};

	// defaults
	reftol = 1.0;
	maxlvl = 0;
	rst = NULL;
	ctree=NULL;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OO|O!di", kwlist, &mesh, &sPos, &RStarTreeType, &rst, &reftol, &maxlvl))
		return NULL;


	//psi_rvec samppos2;
	//if(!PyArg_ParseTuple(sPos, "(ddd)(ddd)", &samppos.x, &samppos.y, &samppos.z, &samppos2.x, &samppos2.y, &samppos2.z))
		//return NULL;
	if(!PyArg_ParseTuple(sPos, "ddd", &samppos.x, &samppos.y, &samppos.z)) {
		return NULL;
	}

	// extract C pointers and such
	PSI_Mesh2mesh(mesh, &cmesh);

	psi_real* rhoout;
	psi_rvec* velout;
	psi_int nsamp;

	if(rst)
		ctree = &rst->ctree;

	// call the C function
	// TODO: allow passing of a user-made r* tree
	psi_sample_vdf(samppos, &cmesh, ctree, 
		&rhoout, &velout, &nsamp, reftol, maxlvl);

	// make the return arrays from the c buffers
	// force numpy to own the memory
	npy_intp npdims[3];
	npdims[0] = nsamp;
	npdims[1] = 3;
	PyObject* pyrho = PyArray_SimpleNewFromData(1, npdims, NPY_DOUBLE, rhoout);
	PyObject* pyvel = PyArray_SimpleNewFromData(2, npdims, NPY_DOUBLE, velout);
	PyArray_ENABLEFLAGS(pyrho, NPY_ARRAY_OWNDATA);
	PyArray_ENABLEFLAGS(pyvel, NPY_ARRAY_OWNDATA);
	return Py_BuildValue("OO", pyrho, pyvel);
}


static PyObject *PSI_crossStreams(PyObject *self, PyObject *args, PyObject* kwds) {

	psi_int nstreams, i, j, ax, ii, jj;
	psi_real rhotot, vij, myxsn;
	psi_rvec samppos;
	psi_mesh cmesh;
	psi_rtree* ctree;
	Mesh* mesh;
	PyObject* nprho;
    PyObject* npvel; 
    PyObject* xsnfunc; 
	RStarTree* rst;
	npy_intp nverts;
	npy_intp npdims[3];
    npy_intp* arshape;
	static char *kwlist[] = {"rho", "vel", "xsnfunc", NULL};

    xsnfunc = NULL;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OO|O", kwlist, &nprho, &npvel, &xsnfunc))
		return NULL;

    // get C pointers to array data
    // TODO: check array types!!!
    // TODO: retrieve dimensions
    arshape = PyArray_SHAPE(nprho); 
    nstreams = arshape[0];
    psi_real* rho = PyArray_DATA(nprho); 
    psi_rvec* vel = PyArray_DATA(npvel); 

    // compute the total density and bulk velocity
    rhotot = 0.0;
    npdims[0] = 3;
    PyObject* pyvel = PyArray_SimpleNew(1, npdims, NPY_DOUBLE);
    psi_real* vtot = (psi_real*) PyArray_DATA(pyvel); 
    for(ax = 0; ax < 3; ++ax)
        vtot[ax] = 0.0;
    for(i = 0; i < nstreams; ++i) {
        rhotot += rho[i];
        for(ax = 0; ax < 3; ++ax)
            vtot[ax] += rho[i]*vel[i].xyz[ax]; 
    }
    for(ax = 0; ax < 3; ++ax)
        vtot[ax] /= rhotot; 

    // velocity dispersion matrix
    npdims[0] = 3;
    npdims[1] = 3;
    PyObject* pycov = PyArray_SimpleNew(2, npdims, NPY_DOUBLE);
	psi_real* cov = (psi_real*) PyArray_DATA(pycov); 
    memset(cov, 0, 9*sizeof(psi_real));
    for(i = 0; i < nstreams; ++i) {
        for(ii = 0; ii < 3; ++ii)
        for(jj = 0; jj < 3; ++jj)
            cov[3*ii+jj] += rho[i]*(vel[i].xyz[ii]-vtot[ii])*(vel[i].xyz[jj]-vtot[jj]);
    }
    for(ii = 0; ii < 3; ++ii)
    for(jj = 0; jj < 3; ++jj)
        cov[3*ii+jj] /= rhotot; 


    // lastly do the velocity-dependent annihilation rate
    psi_real annihilation_rate = 0.0;
    for(i = 0; i < nstreams; ++i) 
    for(j = 0; j < nstreams; ++j) {
        vij = sqrt((vel[i].x-vel[j].x)*(vel[i].x-vel[j].x)
                +(vel[i].y-vel[j].y)*(vel[i].y-vel[j].y)
                +(vel[i].z-vel[j].z)*(vel[i].z-vel[j].z));
        myxsn = 1.0; 
        if(xsnfunc)
            myxsn = PyFloat_AsDouble(PyObject_CallObject(xsnfunc, Py_BuildValue("(d)", vij)));
        annihilation_rate += 0.5*rho[i]*rho[j]*myxsn;
    } 

	return Py_BuildValue("idOOd", nstreams, rhotot, pyvel, pycov, annihilation_rate);
}




static PyObject *PSI_phi(PyObject *self, PyObject *args, PyObject* kwds) {

	psi_int ax;
	npy_intp npdims[3];
	psi_grid cgrid;
	psi_real Gn;
	Grid* grid;	
	static char *kwlist[] = {"grid", "Gn", NULL};

#ifdef HAVE_FFTW
	Gn = 1.0;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|d", kwlist, &grid, &Gn))
		return NULL;

	PSI_Grid2grid(grid, &cgrid);
	if(cgrid.type != PSI_GRID_CART) 
		return NULL;

	for(ax = 0; ax < 3; ++ax)
		npdims[ax] = cgrid.n.ijk[ax];

	// the return array
	PyObject* retar = PyArray_ZEROS(cgrid.dim, npdims, NPY_DOUBLE, 0);
	psi_do_phi(&cgrid, PyArray_DATA(retar), Gn);
	return retar;
#else
	PyErr_SetString(PyExc_RuntimeWarning, 
			"PSI was compiled without FFTW. PSI.phi() does nothing.");
	return NULL;
#endif

}

static PyObject *PSI_powerSpectrum(PyObject *self, PyObject *args, PyObject* kwds) {

	psi_int ax, nbins, nmin;
	psi_real lmin, lside, dk;
	npy_intp npdims[3];
	psi_grid cgrid;
	Grid* grid;	
	static char *kwlist[] = {"grid", "nbins", NULL};

#ifdef HAVE_FFTW
	nbins = -1;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|i", kwlist, &grid, &nbins))
		return NULL;

	PSI_Grid2grid(grid, &cgrid);
	if(cgrid.type != PSI_GRID_CART) 
		return NULL;

	// get the shortest physical box side for the Nyquist frequency
	// if no nbins was specified, use that box side as well
	lmin = 1.0e30;
	for(ax = 0; ax < cgrid.dim; ++ax) {
		lside = cgrid.window[1].xyz[ax]-cgrid.window[0].xyz[ax];
		if(lside < lmin) {
			lmin = lside;
			nmin = cgrid.n.ijk[ax];
		}
	}
	if(nbins < 0) nbins = nmin; 

	// the return array
	npdims[0] = nbins;
	PyObject* pspec = PyArray_ZEROS(1, npdims, NPY_DOUBLE, 0);
	PyObject* kk = PyArray_ZEROS(1, npdims, NPY_DOUBLE, 0);
	psi_do_power_spectrum(&cgrid, PyArray_DATA(pspec), PyArray_DATA(kk), TWO_PI/lmin, nbins);
	return Py_BuildValue("OO", pspec, kk);
#else
	PyErr_SetString(PyExc_RuntimeWarning, 
			"PSI was compiled without FFTW. PSI.powerSpectrum() does nothing.");
	return NULL;
#endif

}



/////////////////////////////////////////////////////////////
//     PSI module init 
/////////////////////////////////////////////////////////////

/* Docstrings */
//static char module_docstring[] = "PSI.";
//static char PSI_docstring[] =
    //"Calculate the chi-squared of some data given a model.";
#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif

#if PY_MAJOR_VERSION >= 3
static PyMethodDef module_methods[] = {
    {"skymap", (PyCFunction)PSI_skymap, METH_KEYWORDS | METH_VARARGS, "Makes a skymap"},
    {"beamtrace", (PyCFunction)PSI_beamtrace, METH_KEYWORDS | METH_VARARGS, "beamtrace"},
    {"voxels", (PyCFunction)PSI_voxels, METH_KEYWORDS | METH_VARARGS, "Voxelizes"},
    {"phi", (PyCFunction)PSI_phi, METH_KEYWORDS | METH_VARARGS, "phi"},
    {"VDF", (PyCFunction)PSI_VDF, METH_KEYWORDS | METH_VARARGS, "VDF"},
    {"crossStreams", (PyCFunction)PSI_crossStreams, METH_KEYWORDS | METH_VARARGS, "crossStreams"},
    {"powerSpectrum", (PyCFunction)PSI_powerSpectrum, METH_KEYWORDS | METH_VARARGS, "powerSpectrum"},
    {NULL, NULL, 0, NULL}
};
#else
static PyMethodDef module_methods[] = {
   	{"skymap", (PyCFunction)PSI_skymap, METH_KEYWORDS, "Makes a skymap"},
   	{"beamtrace", (PyCFunction)PSI_beamtrace, METH_KEYWORDS, "beamtrace"},
   	{"voxels", (PyCFunction)PSI_voxels, METH_KEYWORDS, "Voxelizes"},
   	{"phi", (PyCFunction)PSI_phi, METH_KEYWORDS, "phi"},
   	{"VDF", (PyCFunction)PSI_VDF, METH_KEYWORDS, "VDF"},
   	{"crossStreams", (PyCFunction)PSI_crossStreams, METH_KEYWORDS, "crossStreams"},
   	{"powerSpectrum", (PyCFunction)PSI_powerSpectrum, METH_KEYWORDS, "powerSpectrum"},
    {NULL, NULL, 0, NULL}
};
#endif

#if PY_MAJOR_VERSION >= 3
    static struct PyModuleDef PSI = {
        PyModuleDef_HEAD_INIT,
        "PSI",     /* m_name */
        "PSI, the phase space intersector",  /* m_doc */
        -1,                  /* m_size */
        module_methods    /* m_methods */
    };
#endif

#if PY_MAJOR_VERSION >= 3
PyMODINIT_FUNC PyInit_PSI(void)
#else
PyMODINIT_FUNC initPSI(void)
#endif
{
    PyObject* m;
#if PY_MAJOR_VERSION >= 3
    m = PyModule_Create(&PSI);
#else
    m = Py_InitModule3("PSI", module_methods,
          "PSI, the phase space intersector");
#endif
    if(!m) return;

	// import numpy functionality
	_import_array();

	// add the mesh type
	if(PyType_Ready(&MeshType) < 0) return;
	Py_INCREF(&MeshType);
	PyModule_AddObject(m, "Mesh", (PyObject*)&MeshType);


	// add the grid type
	if(PyType_Ready(&GridType) < 0) return;
	Py_INCREF(&GridType);
	PyModule_AddObject(m, "Grid", (PyObject*)&GridType);

	// add the rtree type
	if(PyType_Ready(&RStarTreeType) < 0) return;
	Py_INCREF(&RStarTreeType);
	PyModule_AddObject(m, "RStarTree", (PyObject*)&RStarTreeType);

#if 0
	// add the metric type
	if(PyType_Ready(&MetricType) < 0) return;
	Py_INCREF(&MetricType);
	PyModule_AddObject(m, "Metric", (PyObject*)&MetricType);
#endif

#if PY_MAJOR_VERSION >= 3
	return m;
#endif

	// TODO: add macroed constants to the Python module 
}

#endif // PYMODULE















#if 0

/////////////////////////////////////////////////////////////
//     Metric
/////////////////////////////////////////////////////////////

static int Metric_init(Metric *self, PyObject *args, PyObject *kwds) {

	printf("Init metric\n");

	psi_int ax, dim, order, nside, npix;
	psi_metric cmetric;
	char* type;
	npy_intp npdims[3];
	PyObject* pytype, *params, *filepat;
	static char *flrwkw[] = {"type", "params", "filepattern", NULL};
	static char *kerrkw[] = {"type", NULL};
	static char *minkkw[] = {"type", NULL};

	// default is 0 (Minkowski) 
	memset(npdims, 0, sizeof(npdims));
	memset(&cmetric, 0, sizeof(cmetric));

	// peek at the metric type before parsing args
	type = PyString_AsString(PyDict_GetItemString(kwds, "type"));
	//if(strcmp(type, "flrw") == 0) {
		////self->type = 
		//printf("Flrw selected\n");		
	//}

	if(strcmp(type, "minkowski") == 0 &&
			PyArg_ParseTupleAndKeywords(args, kwds, "S", minkkw, &pytype)) {
		self->type = pytype;
		printf("Minkowski selected, %s\n", PyString_AsString(pytype));		
	}
	if(strcmp(type, "kerr") == 0 &&
			PyArg_ParseTupleAndKeywords(args, kwds, "S", kerrkw, &pytype)) {
		self->type = pytype;
		printf("Kerr selected, %s\n", PyString_AsString(pytype));		
	}
	else if(strcmp(type, "flrw") == 0 &&
			PyArg_ParseTupleAndKeywords(args, kwds, "S|OS", flrwkw, &pytype, &params, &filepat)) {
		self->type = pytype;
		printf("FLRW selected, %s\n", PyString_AsString(pytype));		
	}
	else {
	
		printf("Bad metric!!!!\n");
		return -1;
	}

#if 0

	// parse arguments differently depending on what the metric type is
	// certain metric types must correspond to certain arg patterns
	if(strcmp(type, "flrw") == 0 && 
			PyArg_ParseTupleAndKeywords(args, kwds, "S|((ddd)(ddd))(iii)", kwlist, // for cart 
			&pytype, &cmetric.window[0].x, &cmetric.window[0].y, &cmetric.window[0].z, 
			&cmetric.window[1].x, &cmetric.window[1].y, &cmetric.window[1].z, &cmetric.n.i, &cmetric.n.j, &cmetric.n.k)) {

		// fill in all metric information as Python tuples
		self->type = pytype; 
		self->winmin = Py_BuildValue("(ddd)", cmetric.window[0].x, cmetric.window[0].y, cmetric.window[0].z); 
		self->winmax = Py_BuildValue("(ddd)", cmetric.window[1].x, cmetric.window[1].y, cmetric.window[1].z); 
		self->n = Py_BuildValue("(iii)", cmetric.n.i, cmetric.n.j, cmetric.n.k); 
		self->d = Py_BuildValue("(ddd)", (cmetric.window[1].x-cmetric.window[0].x)/cmetric.n.i, (cmetric.window[1].y-cmetric.window[0].y)/cmetric.n.j, (cmetric.window[1].z-cmetric.window[0].z)/cmetric.n.k); 
		npdims[0] = cmetric.n.i;
		npdims[1] = cmetric.n.j;
		npdims[2] = cmetric.n.k;
		dim = 3;
		Py_XDECREF(pytype);
	}
	else if(strcmp(type, "hpring") == 0 && 
			PyArg_ParseTupleAndKeywords(args, kwds, "S|i", kwlist_hpring, &pytype, &cmetric.n.i)) {

		// for Healpix store n as (order, nside, npix)
		self->type = pytype; 
		order = floor(log2(cmetric.n.i+0.5));
		nside = (1<<order);
		npix = 12*nside*nside;
		self->n = Py_BuildValue("(iii)", order, nside, npix); 
		self->d = Py_BuildValue("d", FOUR_PI/npix); // the nominal pixel area 
		self->winmin = Py_None;
		self->winmax = Py_None;
		npdims[0] = npix;
		dim = 1;
		Py_XDECREF(pytype);
	}
	else {
		return -1;
	}

	// make the fields dict
	// now allocate numpy storage
	self->fields = PyDict_New();
	PyDict_SetItemString(self->fields, "m", PyArray_ZEROS(dim, npdims, NPY_DOUBLE, 0));

#endif
    return 0;
}

static PyObject* Metric_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    Metric *self;
    self = (Metric*)type->tp_alloc(type, 0);
    return (PyObject *)self;
}

static void Metric_dealloc(Metric* self) {
    Py_XDECREF(self->type);
    Py_XDECREF(self->params);
    Py_XDECREF(self->n);
    Py_XDECREF(self->box);
    Py_XDECREF(self->phi);
    Py_XDECREF(self->gradphi);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyMemberDef Metric_members[] = {
    {"type", T_OBJECT, offsetof(Metric, type), 0,
     "type"},
    {"params", T_OBJECT, offsetof(Metric, params), 0,
     "params"},
    {"n", T_OBJECT, offsetof(Metric, n), 0,
     "n"},
    {"box", T_OBJECT, offsetof(Metric, box), 0,
     "box"},
    {"phi", T_OBJECT, offsetof(Metric, phi), 0,
     "potential"},
    {"gradphi", T_OBJECT, offsetof(Metric, gradphi), 0,
     "potential gradient"},
    {NULL}  /* Sentinel */
};

static PyMethodDef Metric_methods[] = {
	//{"getCellGeometry", (PyCFunction)Metric_getCellGeometry, METH_KEYWORDS,
	 //"Get the vertices for the given cell."},
    {NULL}  /* Sentinel */
};

static PyTypeObject MetricType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "PSI.Metric",             /* tp_name */
    sizeof(Metric),             /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)Metric_dealloc, /* tp_dealloc */
    0,                         /* tp_print */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_compare */
    0,                         /* tp_repr */
    0,                         /* tp_as_number */
    0,                         /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash */
    0,                         /* tp_call */
    0,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT |
        Py_TPFLAGS_BASETYPE,   /* tp_flags */
    "Basic metric",           /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    Metric_methods,             /* tp_methods */
    Metric_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)Metric_init,      /* tp_init */
    0,                         /* tp_alloc */
    Metric_new,                 /* tp_new */
};



#endif


