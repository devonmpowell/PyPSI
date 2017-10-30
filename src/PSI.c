// Python interface only
#ifdef PYMODULE 

#include <numpy/arrayobject.h>
#include "structmember.h"
#include "psi.h"
#include "grid.h"
#include "mesh.h"
#include "skymap.h"
#include "beamtrace.h"

#define MAKE_PY_NONE Py_BuildValue("")


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

	// temporary, just to get the mass array
	cgrid->fields[0] = (psi_real*)PyArray_DATA((PyArrayObject*)PyDict_GetItemString(grid->fields, "m"));
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
		Py_RETURN_NONE;
	}
}



/////////////////////////////////////////////////////////////
//    Mesh 
/////////////////////////////////////////////////////////////

static int Mesh_init(Mesh *self, PyObject *args, PyObject *kwds) {

	psi_int ax, ndim;
	psi_mesh cmesh;
	PyObject *posar, *velar, *massar, *box, *n;
	char* cloader, *cfile;
	npy_intp npdims[6];
	npy_intp *dtmp; 
	static char *kwlist[] = {"loader", "filename", "posar", "velar", "massar", "box", "n", NULL};

	// parse arguments differently depending on what the grid type is
	// certain grid types must correspond to certain arg patterns
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "s|sOOOOO", kwlist, &cloader, &cfile, &posar, &velar, &massar, &box, &n))
			return -1;

	// fill in all grid information as Python tuples
	if(strcmp(cloader, "array") == 0) {

		// TODO: check dtype!!!!

		setbuf(stdout, NULL);

		if(!PyArray_Check(posar))
			return -1;

		ndim = PyArray_NDIM(posar);
		dtmp = PyArray_DIMS(posar);

		psi_printf("Pos array dimensions are %d %d %d %d, ndim = %d\n", dtmp[0], dtmp[1], dtmp[2], dtmp[3], ndim);

		self->boxmin = PySequence_GetItem(box, 0); 
		self->boxmax = PySequence_GetItem(box, 1); 

		self->pos = posar;
		Py_INCREF(posar);

		// TODO: set better error handling here
		self->vel = PyArray_SimpleNew(ndim, dtmp, NPY_DOUBLE);
		self->mass = PyArray_SimpleNew(ndim-1, dtmp, NPY_DOUBLE);


		psi_int nside = 64; 
		//psi_int nside = 512; 
		npdims[0] = nside*nside*nside;
		npdims[1] = psi_verts_per_elem(PSI_MESH_LINEAR); 
		self->connectivity = PyArray_SimpleNew(2, npdims, NPY_INT32);
		psi_int* cptr = (psi_int*) PyArray_DATA(self->connectivity); 

		psi_printf("Using the array loader.\n");

		// now build the mesh connectivity
		// trilinear elements naturally
		psi_int i, j, k, ii, jj, kk, locind, vertind;
		psi_int elemind;
		for(i = 0; i < nside; ++i)
		for(j = 0; j < nside; ++j)
		for(k = 0; k < nside; ++k) {
			elemind = nside*nside*i + nside*j + k;
			for(ii = 0; ii < 2; ++ii)
			for(jj = 0; jj < 2; ++jj)
			for(kk = 0; kk < 2; ++kk) {
				locind = 4*ii + 2*jj + kk;
				vertind = nside*nside*((i+ii)%nside) 
					+ nside*((j+jj)%nside) + ((k+kk)%nside);
				cptr[8*elemind+locind] = vertind;
			}
		}

	}
	else if(strcmp(cloader, "gadget2") == 0) {

		psi_printf("Using the Gadget2 loader.\n");

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
			self->boxmin = MAKE_PY_NONE;
			self->boxmax = MAKE_PY_NONE;
		}

		cmesh.pos = PyArray_DATA((PyArrayObject*)self->pos);
		cmesh.vel = PyArray_DATA((PyArrayObject*)self->vel);
		cmesh.mass = PyArray_DATA((PyArrayObject*)self->mass);
		cmesh.connectivity = PyArray_DATA((PyArrayObject*)self->connectivity);

		// load the data into the numpy buffers
		if(!load_gadget2(&cmesh, cfile))
			return -1;

	}
	else if(strcmp(cloader, "hacky_test") == 0) {

		psi_printf("Using the hacked test vertices.\n");

		// set up the cmesh...
		cmesh.npart = 8; 
		cmesh.nelem = 2; 
		cmesh.periodic = 1;
		cmesh.elemtype = PSI_MESH_SIMPLEX;
		cmesh.dim = 3; 
		for(ax = 0; ax < 3; ++ax) {
			cmesh.box[0].xyz[ax] = 0.0;
			cmesh.box[1].xyz[ax] = 40.0; 
		}

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
			self->boxmin = MAKE_PY_NONE;
			self->boxmax = MAKE_PY_NONE;
		}

		cmesh.pos = PyArray_DATA((PyArrayObject*)self->pos);
		cmesh.vel = PyArray_DATA((PyArrayObject*)self->vel);
		cmesh.mass = PyArray_DATA((PyArrayObject*)self->mass);
		cmesh.connectivity = PyArray_DATA((PyArrayObject*)self->connectivity);
	
		//  fill in the vertex pos, mass, vel...
		cmesh.mass[0] = 1.0;
		cmesh.mass[1] = 1.0;
		cmesh.mass[2] = 1.0;
		cmesh.mass[3] = 1.0;
		cmesh.mass[4] = 1.0;
		cmesh.mass[5] = 1.0;
		cmesh.mass[6] = 1.0;
		cmesh.mass[7] = 1.0;
		cmesh.pos[0].x = 17.0;
		cmesh.pos[0].y = 17.0;
		cmesh.pos[0].z = 17.0;
		cmesh.pos[1].x = 29.0;
		cmesh.pos[1].y = 17.0;
		cmesh.pos[1].z = 17.0;
		cmesh.pos[2].x = 17.0;
		cmesh.pos[2].y = 29.0;
		cmesh.pos[2].z = 17.0;
		cmesh.pos[3].x = 17.0;
		cmesh.pos[3].y = 17.0;
		cmesh.pos[3].z = 29.0;
		cmesh.pos[4].x = 25.0;
		cmesh.pos[4].y = 25.0;
		cmesh.pos[4].z = 25.0;
		cmesh.pos[5].x = 9.0;
		cmesh.pos[5].y = 25.0;
		cmesh.pos[5].z = 25.0;
		cmesh.pos[6].x = 25.0;
		cmesh.pos[6].y = 9.0;
		cmesh.pos[6].z = 25.0;
		cmesh.pos[7].x = 25.0;
		cmesh.pos[7].y = 25.0;
		cmesh.pos[7].z = 9.0;
		cmesh.connectivity[0] = 0;
		cmesh.connectivity[1] = 1;
		cmesh.connectivity[2] = 2;
		cmesh.connectivity[3] = 3;
		cmesh.connectivity[4] = 4;
		cmesh.connectivity[5] = 5;
		cmesh.connectivity[6] = 6;
		cmesh.connectivity[7] = 7;

	}

	else {
		psi_printf("Bad loader: %s\n", cloader);
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

			break;


		default:
			psi_printf("Bad grid type\n");
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
		cgrid.window[0].xyz[ax] = 0.0;
		cgrid.window[1].xyz[ax] = 1.0;
		cgrid.n.ijk[ax] = 32;
	}

	// peek at the grid type before parsing args
	// parse arguments differently depending on what the grid type is
	// certain grid types must correspond to certain arg patterns
	type = PyString_AsString(PyDict_GetItemString(kwds, "type"));
	if(strcmp(type, "cart") == 0 && 
			PyArg_ParseTupleAndKeywords(args, kwds, "S|((ddd)(ddd))(iii)", kwlist, // for cart 
			&pytype, &cgrid.window[0].x, &cgrid.window[0].y, &cgrid.window[0].z, 
			&cgrid.window[1].x, &cgrid.window[1].y, &cgrid.window[1].z, &cgrid.n.i, &cgrid.n.j, &cgrid.n.k)) {

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
		return MAKE_PY_NONE;

	// make C structs, check them
	PSI_Grid2grid(grid, &cgrid);
	PSI_Mesh2mesh(mesh, &cmesh);
	if(cgrid.type != PSI_GRID_HPRING || cgrid.dim != 3) 
		Py_RETURN_NONE;
	if(cmesh.dim != cgrid.dim)
		Py_RETURN_NONE;

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
		return MAKE_PY_NONE;

	PSI_Grid2grid(grid, &cgrid);
	if(cgrid.type != PSI_GRID_HPRING || cgrid.dim != 3) 
		Py_RETURN_NONE;

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

	psi_grid cgrid;
	psi_mesh cmesh;
	Mesh* mesh;
	Grid* grid;	
	//npy_intp* npdims;
	npy_intp nverts;
	static char *kwlist[] = {"grid", "mesh", NULL};

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OO", kwlist, &grid, &mesh))
		return MAKE_PY_NONE;

	PSI_Grid2grid(grid, &cgrid);
	PSI_Mesh2mesh(mesh, &cmesh);
	if(cgrid.type != PSI_GRID_CART) 
		Py_RETURN_NONE;
	if(cmesh.dim != cgrid.dim)
		Py_RETURN_NONE;

	psi_voxels(&cgrid, &cmesh);

   Py_RETURN_NONE;
}

static PyObject *PSI_phi(PyObject *self, PyObject *args, PyObject* kwds) {

	psi_grid cgrid;
	psi_real Gn;
	Grid* grid;	

	static char *kwlist[] = {"grid", "Gn", NULL};

	Gn = 1.0;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|d", kwlist, &grid, &Gn))
		return MAKE_PY_NONE;


	PSI_Grid2grid(grid, &cgrid);
	if(cgrid.type != PSI_GRID_CART) 
		Py_RETURN_NONE;

	psi_int ax;
	npy_intp npdims[3];
	for(ax = 0; ax < 3; ++ax)
		npdims[ax] = cgrid.n.ijk[ax];

	// the return array
	PyObject* retar = PyArray_ZEROS(cgrid.dim, npdims, NPY_DOUBLE, 0);
	//psi_do_phi(&cgrid, PyArray_DATA(retar), Gn);
	return retar;
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
   	{"voxels", (PyCFunction)PSI_voxels, METH_KEYWORDS, "Voxelizes"},
   	{"phi", (PyCFunction)PSI_phi, METH_KEYWORDS, "phi"},
   	{"beamtrace", (PyCFunction)PSI_beamtrace, METH_KEYWORDS, "beamtrace"},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initPSI(void) {
    PyObject* m;
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

	// add the metric type
	if(PyType_Ready(&MetricType) < 0) return;
	Py_INCREF(&MetricType);
	PyModule_AddObject(m, "Metric", (PyObject*)&MetricType);


	// TODO: add macroed constants to the Python module 
}

#endif // PYMODULE
