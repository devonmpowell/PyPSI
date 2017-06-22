// Python interface only
#ifdef PYMODULE 

#include <numpy/arrayobject.h>
#include "structmember.h"
#include "psi.h"
#include "grid.h"
#include "mesh.h"
#include "skymap.h"

#define MAKE_PY_NONE Py_BuildValue("")

/////////////////////////////////////////////////////////////
//   Define Python type structs and conversion functions up front 
/////////////////////////////////////////////////////////////

typedef struct {
    PyObject_HEAD
    PyObject *type; 
    PyObject *filename; 
    PyArrayObject *pos, *vel, *mass;
} Mesh;

typedef struct {
    PyObject_HEAD
	PyObject* type; // string representation of the grid type 
	PyObject* fields; // Python dict of numpy arrays 
	PyObject* winmin, *winmax; // Python tuple of the projection window
	PyObject* n; // number of cells in each dimension
	PyObject* d; // tuple of dx in each dimension
} Grid;

static void PSI_Grid2grid(Grid* grid, psi_grid* cgrid) {

	// parse a simple C struct from a Python object 
	psi_int ax;
	PyObject* seq;
	cgrid->type = PyString_AsString(grid->type);
	seq = PySequence_Fast(grid->n, "expected a grid dim tuple");
	for(ax = 0; ax < 3; ++ax)
		cgrid->n.ijk[ax] = PyInt_AS_LONG(PySequence_Fast_GET_ITEM(seq, ax));
	if(PySequence_Check(grid->winmin)) {
		seq = PySequence_Fast(grid->winmin, "expected a grid window tuple");
		for(ax = 0; ax < 3; ++ax)
			cgrid->winmin.xyz[ax] = PyFloat_AS_DOUBLE(PySequence_Fast_GET_ITEM(seq, ax));
		seq = PySequence_Fast(grid->winmax, "expected a grid window tuple");
		for(ax = 0; ax < 3; ++ax)
			cgrid->winmax.xyz[ax] = PyFloat_AS_DOUBLE(PySequence_Fast_GET_ITEM(seq, ax));
	}
	if(PySequence_Check(grid->d)) {
		seq = PySequence_Fast(grid->d, "expected a grid dx tuple");
		for(ax = 0; ax < 3; ++ax)
			cgrid->d.xyz[ax] = PyFloat_AS_DOUBLE(PySequence_Fast_GET_ITEM(seq, ax));
	}
	else if(PyFloat_Check(grid->d)) {
		cgrid->d.x = PyFloat_AsDouble(grid->d);
	}

	// temporary, just to get the mass array
	cgrid->fields[0] = (psi_real*)PyArray_DATA((PyArrayObject*)PyDict_GetItemString(grid->fields, "m"));
}


/////////////////////////////////////////////////////////////
//    Mesh 
/////////////////////////////////////////////////////////////

static int Mesh_init(Mesh *self, PyObject *args, PyObject *kwds) {

	psi_int nside;
	psi_mesh cmesh;
	npy_intp npdims[6];
	PyObject* pyfile, *pytype;
	static char *kwlist[] = {"filename", "type", NULL};

	// parse arguments differently depending on what the grid type is
	// certain grid types must correspond to certain arg patterns
	pyfile = NULL;
	pytype = NULL;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "SS", kwlist, &pyfile, &pytype))
			return -1;

	// fill in all grid information as Python tuples
	self->type = pytype;
	self->filename = pyfile;
	cmesh.type = PyString_AsString(pytype);
	if(strcmp(cmesh.type, "gadget2") == 0) {

		// peek at the file first to make sure it's there
		// and to get header information
		if(!peek_gadget2(&cmesh, PyString_AsString(pyfile)))
			return -1;

		// wrap the cmesh data in a cubical Numpy array
		nside = floor(pow(cmesh.npart+0.5, ONE_THIRD));
		npdims[0] = nside;
		npdims[1] = nside;
		npdims[2] = nside;
		npdims[3] = 3;
		self->pos = (PyArrayObject*)PyArray_SimpleNew(4, npdims, NPY_DOUBLE);
		self->vel = (PyArrayObject*)PyArray_SimpleNew(4, npdims, NPY_DOUBLE);
		if(!self->pos || !self->vel) {
			printf("new array fail!\n");
			return -1;
		} 
		cmesh.pos = PyArray_DATA(self->pos);
		cmesh.vel = PyArray_DATA(self->vel);

		// load the data into the numpy buffers
		if(!load_gadget2(&cmesh, PyString_AsString(pyfile)))
			return -1;
	}
	else {
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
    Py_XDECREF(self->type);
    Py_XDECREF(self->filename);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyMemberDef Mesh_members[] = {
    {"type", T_OBJECT_EX, offsetof(Mesh, type), 0, "type"},
    {"filename", T_OBJECT_EX, offsetof(Mesh, filename), 0, "filename"},
    {"pos", T_OBJECT_EX, offsetof(Mesh, pos), 0, "pos"},
    {"vel", T_OBJECT_EX, offsetof(Mesh, vel), 0, "vel"},
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
		return MAKE_PY_NONE;
	allcells = (pyind == Py_None)? 1 : 0; 
	if(!allcells && !PySequence_Check(pyind)) 
		return MAKE_PY_NONE;


	// parse the grid to the simple C struct
	// make the return tuple based on the number of input elements
	// process for different indexing protocols
	PSI_Grid2grid(self, &cgrid);
	if(strcmp(cgrid.type, "cart3d") == 0) {
	
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

		

	}
	else if(strcmp(cgrid.type, "hpring") == 0) {

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
				return MAKE_PY_NONE;
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
				return MAKE_PY_NONE;
			} 
			for(ax = 0; ax < 3; ++ax)
				cvdata[3*i+ax] = center.xyz[ax];
			for(v = 0; v < 4*bstep; ++v)
				for(ax = 0; ax < 3; ++ax)
					bvdata[12*bstep*i+3*v+ax] = boundary[v].xyz[ax];
			pvdata[i] = vol; 
		}
	}
	else {
		// TODO: error handling
		printf("Bad grid type\n");
		return MAKE_PY_NONE;
	}
	
	return Py_BuildValue("(NNN)", PyArray_Return((PyArrayObject*)cverts), 
			PyArray_Return((PyArrayObject*)bverts), PyArray_Return((PyArrayObject*)pyvol)); 
}


static int Grid_init(Grid *self, PyObject *args, PyObject *kwds) {

	psi_int ax, dim, order, nside, npix;
	psi_grid cgrid;
	char* type;
	npy_intp npdims[3];
	PyObject* pytype;
	static char *kwlist[] = {"type", "window", "n", NULL};
	static char *kwlist_hpring[] = {"type", "n", NULL};

	// default is a 32^3 unit box 
	// or healpix with nside=32
	dim = 0;
	memset(npdims, 0, sizeof(npdims));
	for(ax = 0; ax < 3; ++ax) {
		cgrid.winmin.xyz[ax] = 0.0;
		cgrid.winmax.xyz[ax] = 1.0;
		cgrid.n.ijk[ax] = 32;
	}

	// peek at the grid type before parsing args
	// parse arguments differently depending on what the grid type is
	// certain grid types must correspond to certain arg patterns
	type = PyString_AsString(PyDict_GetItemString(kwds, "type"));
	if(strcmp(type, "cart3d") == 0 && 
			PyArg_ParseTupleAndKeywords(args, kwds, "S|((ddd)(ddd))(iii)", kwlist, // for cart3d 
			&pytype, &cgrid.winmin.x, &cgrid.winmin.y, &cgrid.winmin.z, 
			&cgrid.winmax.x, &cgrid.winmax.y, &cgrid.winmax.z, &cgrid.n.i, &cgrid.n.j, &cgrid.n.k)) {

		// fill in all grid information as Python tuples
		self->type = pytype; 
		self->winmin = Py_BuildValue("(ddd)", cgrid.winmin.x, cgrid.winmin.y, cgrid.winmin.z); 
		self->winmax = Py_BuildValue("(ddd)", cgrid.winmax.x, cgrid.winmax.y, cgrid.winmax.z); 
		self->n = Py_BuildValue("(iii)", cgrid.n.i, cgrid.n.j, cgrid.n.k); 
		self->d = Py_BuildValue("(ddd)", (cgrid.winmax.x-cgrid.winmin.x)/cgrid.n.i, (cgrid.winmax.y-cgrid.winmin.y)/cgrid.n.j, (cgrid.winmax.z-cgrid.winmin.z)/cgrid.n.k); 
		npdims[0] = cgrid.n.i;
		npdims[1] = cgrid.n.j;
		npdims[2] = cgrid.n.k;
		dim = 3;
		Py_XDECREF(pytype);
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
		self->winmin = MAKE_PY_NONE;
		self->winmax = MAKE_PY_NONE;
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

static PyMethodDef Grid_methods[] = {
	{"getCellGeometry", (PyCFunction)Grid_getCellGeometry, METH_KEYWORDS,
	 "Get the vertices for the given cell."},
    {NULL}  /* Sentinel */
};

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
	psi_int bstep;
	Mesh* mesh;
	Grid* grid;	
	//npy_intp* npdims;
	npy_intp nverts;
	static char *kwlist[] = {"grid", "mesh", "bstep", NULL};

	bstep = 1;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OO|i", kwlist, &grid, &mesh, &bstep))
		return MAKE_PY_NONE;

	printf("Got grid verts...\n");
	nverts = PyArray_SIZE((PyArrayObject*)mesh->pos);
	if(nverts%3) return MAKE_PY_NONE;
	nverts /= 3;
	printf("npart = %d\n", (int)nverts);
	//for(p = 0; p < )
	

	PSI_Grid2grid(grid, &cgrid);
	if(strcmp(cgrid.type, "hpring") != 0) 
		Py_RETURN_NONE;

	printf("----Making hpring skymap-----\n");


	psi_skymap(&cgrid, bstep);



   /* Do your stuff here. */
   Py_RETURN_NONE;

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

static PyMethodDef module_methods[] = {
   	{"skymap", (PyCFunction)PSI_skymap, METH_KEYWORDS, "Makes a skymap"},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initPSI(void) {
    PyObject* m;
    if (PyType_Ready(&GridType) < 0)
        return;

    m = Py_InitModule3("PSI", module_methods,
	       "Example module that creates an extension type.");
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


	// TODO: add macroed constants to the module
}

#endif // PYMODULE
