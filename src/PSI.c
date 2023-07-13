// Python interface only
#ifdef PYMODULE 

#include "psi.h"
#include <numpy/arrayobject.h>
#include "structmember.h"
#include "grid.h"
#include "mesh.h"
#include "rtree.h"
#include "skymap.h"
#include "beamtrace.h"
#ifdef HAVE_FFTW
#include "fft.h"
#endif

// TODO: type-checking and error handling all around!
// TODO: ref counts and memory leaks!
// TODO: make grid and mesh immutable!!!


/////////////////////////////////////////////////////////////
//     Grid
/////////////////////////////////////////////////////////////


typedef struct {
    PyObject_HEAD
	PyObject* type; // string representation of the grid type 
	PyObject* fields; // Python dict of numpy arrays 
	PyObject* bounds; // Python tuple of the projection window
	PyObject* n; // number of cells in each dimension
	PyObject* d; // tuple of dx in each dimension
	psi_grid cgrid; // c values and pointers	
} Grid;

static int Grid_init(Grid *self, PyObject *args, PyObject *kwds) {

	psi_int ax, dim, i, len, order, nside, npix;
	const char *ctype, *cstring;
	npy_intp npdims[6];
	PyObject *seq, *tmp, *arr;
	PyObject *type_in = NULL, *bounds_in = NULL, *n_in = NULL, *fields_in = NULL;
	static char *kwlist[] = {"type", "n", "bounds", "fields", NULL};

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "|UOOO", kwlist, &type_in, &n_in, &bounds_in, &fields_in)) 
		return -1;

	if(type_in) {
		ctype = PyUnicode_AsUTF8(type_in);
		if(strcmp(ctype, "cart") == 0) {
			self->cgrid.type = PSI_GRID_CART;
		}
		else if(strcmp(ctype, "hpring") == 0) {
			self->cgrid.type = PSI_GRID_HPRING;
		}
		else {
			self->cgrid.type = -1; 
			PyErr_Format(PyExc_ValueError, "Invalid grid type: %s", ctype);
			return -1;
		}
	    tmp = self->type;
    	self->type = type_in;
    	Py_INCREF(type_in);
    	Py_DECREF(tmp);
	}

	if(self->cgrid.type == PSI_GRID_CART) {

		if(n_in && !PyArg_ParseTuple(n_in, "iii", &self->cgrid.n.i, &self->cgrid.n.j, &self->cgrid.n.k)) {
			PyErr_SetString(PyExc_TypeError, "n must be a tuple containing 3 integers");
			return -1;
		}
		tmp = self->n;
	    self->n = Py_BuildValue("(iii)", self->cgrid.n.i, self->cgrid.n.j, self->cgrid.n.k);
	    Py_DECREF(tmp);

		if(bounds_in && !PyArg_ParseTuple(bounds_in, "(ddd)(ddd)", 
				&self->cgrid.window[0].x, &self->cgrid.window[0].y, &self->cgrid.window[0].z, 
				&self->cgrid.window[1].x, &self->cgrid.window[1].y, &self->cgrid.window[1].z)) {
			PyErr_SetString(PyExc_TypeError, "bounds must be a tuple containing two tuples of 3 values each");
			return -1;
		}
		tmp = self->bounds;
	    self->bounds = Py_BuildValue("((ddd)(ddd))", 
			self->cgrid.window[0].x, self->cgrid.window[0].y, self->cgrid.window[0].z, 
			self->cgrid.window[1].x, self->cgrid.window[1].y, self->cgrid.window[1].z);
	    Py_DECREF(tmp);

		for(ax = 0; ax < self->cgrid.dim; ++ax) {
			self->cgrid.d.xyz[ax] = (self->cgrid.window[1].xyz[ax]-self->cgrid.window[0].xyz[ax])/self->cgrid.n.ijk[ax];
		}
		tmp = self->d;
	    self->d = Py_BuildValue("(ddd)", self->cgrid.d.x, self->cgrid.d.y, self->cgrid.d.z);
	    Py_DECREF(tmp);

		dim = self->cgrid.dim;
		for(ax = 0; ax < self->cgrid.dim; ++ax) {
			npdims[ax] = self->cgrid.n.ijk[ax];
		}

		// always do m
		arr = PyArray_ZEROS(dim, npdims, NPY_FLOAT64, 0);
		if(!arr) return -1;
		self->cgrid.fields[PSI_GRID_M] = (psi_real*)PyArray_DATA((PyArrayObject*)arr);
		if(PyDict_SetItemString(self->fields, "m", arr) < 0) 
			return -1;

		if(fields_in != NULL && fields_in != Py_None) {
		    seq = PySequence_Fast(fields_in, "Expected a sequence for fields.");
			if(!seq) return -1;
		    len = PySequence_Size(fields_in);
		    for (i = 0; i < len; i++) {
		        tmp = PySequence_Fast_GET_ITEM(seq, i);
				if(PyUnicode_Check(tmp)) {
					cstring = PyUnicode_AsUTF8(tmp);
					if(strcmp(cstring, "m") == 0) {;}
					else if(strcmp(cstring, "x") == 0) {
						npdims[3] = 3; 
						dim = 4;
						arr = PyArray_ZEROS(dim, npdims, NPY_FLOAT64, 0);
						if(!arr) return -1;
						self->cgrid.fields[PSI_GRID_X] = (psi_real*)PyArray_DATA((PyArrayObject*)arr);
						if(PyDict_SetItemString(self->fields, "x", arr) < 0) 
							return -1;
					}
					else if(strcmp(cstring, "v") == 0) {
						npdims[3] = 3; 
						dim = 4;
						arr = PyArray_ZEROS(dim, npdims, NPY_FLOAT64, 0);
						if(!arr) return -1;
						self->cgrid.fields[PSI_GRID_V] = (psi_real*)PyArray_DATA((PyArrayObject*)arr);
						if(PyDict_SetItemString(self->fields, "v", arr) < 0) 
							return -1;
					}
					else if(strcmp(cstring, "xx") == 0) {
						npdims[3] = 3;
						npdims[4] = 3; 
						dim = 5;
						arr = PyArray_ZEROS(dim, npdims, NPY_FLOAT64, 0);
						if(!arr) return -1;
						self->cgrid.fields[PSI_GRID_XX] = (psi_real*)PyArray_DATA((PyArrayObject*)arr);
						if(PyDict_SetItemString(self->fields, "xx", arr) < 0) 
							return -1;
					}
					else if(strcmp(cstring, "xv") == 0 || strcmp(cstring, "vx") == 0) {
						npdims[3] = 3; 
						npdims[4] = 3; 
						dim = 5;
						arr = PyArray_ZEROS(dim, npdims, NPY_FLOAT64, 0);
						if(!arr) return -1;
						self->cgrid.fields[PSI_GRID_XV] = (psi_real*)PyArray_DATA((PyArrayObject*)arr);
						if(PyDict_SetItemString(self->fields, "xv", arr) < 0) 
							return -1;
					}
					else if(strcmp(cstring, "vv") == 0) {
						npdims[3] = 3; 
						npdims[4] = 3; 
						dim = 5;
						arr = PyArray_ZEROS(dim, npdims, NPY_FLOAT64, 0);
						if(!arr) return -1;
						self->cgrid.fields[PSI_GRID_VV] = (psi_real*)PyArray_DATA((PyArrayObject*)arr);
						if(PyDict_SetItemString(self->fields, "vv", arr) < 0) 
							return -1;
					}
					else {
						PyErr_Format(PyExc_ValueError, "Unrecognized field identifier: %s", cstring);
						return -1;
					}
				}
				else {
					PyErr_SetString(PyExc_TypeError, "Field identifiers must be strings");
					return -1;
				}
			}
			Py_DECREF(seq);

			// check if we have all the 1st-order moments we need for the 2nd-order ones
			if(self->cgrid.fields[PSI_GRID_XX] && !self->cgrid.fields[PSI_GRID_X]) {
				PyErr_SetString(PyExc_ValueError, "Field \'xx\' requires you to also provide field \'x\'.");
				return -1;
			}
			if(self->cgrid.fields[PSI_GRID_XV] && !(self->cgrid.fields[PSI_GRID_V] && self->cgrid.fields[PSI_GRID_X] )){
				PyErr_SetString(PyExc_ValueError, "Field \'xv\' requires you to also provide fields \'x\' and \'v\'.");
				return -1;
			}
			if(self->cgrid.fields[PSI_GRID_VV] && !self->cgrid.fields[PSI_GRID_V]) {
				PyErr_SetString(PyExc_ValueError, "Field \'vv\' requires you to also provide field \'v\'.");
				return -1;
			}
		}
	}
	else if(self->cgrid.type == PSI_GRID_HPRING) {

		PyErr_SetString(PyExc_ValueError, "Grid type hpring is not yet implemented");
		return -1;
	
#if 0
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
		PyDict_SetItemString(self->fields, "m", PyArray_ZEROS(dim, npdims, NPY_FLOAT64, 0));

#endif
	
	}
	else {
		PyErr_SetString(PyExc_ValueError, "Invalid grid type");
		return -1;
	}

    return 0;
}

static PyObject* Grid_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {

	// set up default c struct
	psi_int ax;
	psi_grid cgrid;
	memset(&cgrid, 0, sizeof(psi_grid));
	char *ctype = "cart";
	cgrid.type = PSI_GRID_CART;
	cgrid.dim = 3;
	for(ax = 0; ax < cgrid.dim; ++ax) {
		cgrid.n.ijk[ax] = 32;
		cgrid.window[0].xyz[ax] = 0.0;
		cgrid.window[1].xyz[ax] = 1.0;
		cgrid.d.xyz[ax] = (cgrid.window[1].xyz[ax]-cgrid.window[0].xyz[ax])/cgrid.n.ijk[ax];
	}

	// copy c struct, initialize Python values to None 
	Grid *self;
	self = (Grid*)type->tp_alloc(type, 0);
	if (self != NULL) {
		self->cgrid = cgrid;
		self->fields = PyDict_New();
		self->type = PyUnicode_FromString(ctype);
		self->bounds = Py_None; 
		Py_INCREF(Py_None);
		self->n = Py_None; 
		Py_INCREF(Py_None);
		self->d = Py_None; 
		Py_INCREF(Py_None);
	}
	return (PyObject *)self;
}

static void Grid_dealloc(Grid* self) {
    Py_XDECREF(self->fields);
    Py_XDECREF(self->bounds);
    Py_XDECREF(self->n);
    Py_XDECREF(self->d);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyMemberDef Grid_members[] = {
    {"type", T_OBJECT, offsetof(Grid, type), 0,
     "type"},
    {"fields", T_OBJECT, offsetof(Grid, fields), 0,
     "fields dict"},
    {"bounds", T_OBJECT, offsetof(Grid, bounds), 0,
     "lower-left and upper-right corners of the grid"},
    {"n", T_OBJECT, offsetof(Grid, n), 0,
     "number of gridpoints"},
    {"d", T_OBJECT, offsetof(Grid, d), 0,
     "grid dx"},
    {NULL}  /* Sentinel */
};

static PyMethodDef Grid_methods[] = {
	//{"getCellGeometry", (PyCFunction)Grid_getCellGeometry, METH_KEYWORDS | METH_VARARGS,
	 //"Get the vertices for the given cell."},
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
//     Mesh
/////////////////////////////////////////////////////////////


typedef struct {
    PyObject_HEAD
    PyObject *pos, *vel, *mass;
    PyObject *connectivity;
    PyObject *wrapbox;
	psi_mesh cmesh;
} Mesh;


static int Mesh_init(Mesh *self, PyObject *args, PyObject *kwds) {

	psi_int ax, ndim, i;
	psi_dvec nside, nelem;
	char* cloader, *cfile;
	npy_intp npdims[6];
	npy_intp *dtmp; 
	PyObject *seq, *tmp, *arr;
	PyObject *posar = NULL, *velar = NULL, *massar = NULL, *box = NULL, *connar = NULL;
	static char *kwlist[] = {"loader", "filename", "pos", "vel", "mass", "connectivity", "periodic_box", NULL};

	// parse arguments differently depending on what the grid type is
	// certain grid types must correspond to certain arg patterns
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "s|sOOOOO", kwlist, &cloader, &cfile, &posar, &velar, &massar, &connar, &box))
			return -1;

	// check if we have a periodic wrap box
	if(strcmp(cloader, "array") == 0 || strcmp(cloader, "block") == 0) {
		if(box) {
			if(!PyArg_ParseTuple(box, "(ddd)(ddd)", 
					&self->cmesh.box[0].x, &self->cmesh.box[0].y, &self->cmesh.box[0].z, 
					&self->cmesh.box[1].x, &self->cmesh.box[1].y, &self->cmesh.box[1].z)) {
				PyErr_SetString(PyExc_TypeError, "periodic_box must be a tuple containing two tuples of 3 values each");
				return -1;
			}
			self->cmesh.periodic = 1;
		}
		else {
			self->cmesh.periodic = 0;
		}
	}

	// fill in all grid information as Python tuples
	if(strcmp(cloader, "array") == 0) {

		if(PyArray_Check(posar)) {

			// get the dimensionality based on shape of pos 
			if(PyArray_TYPE(posar) != NPY_FLOAT64) {
				PyErr_SetString(PyExc_TypeError, "pos must be an array of doubles");
				return -1;
			}
			ndim = PyArray_NDIM((PyArrayObject*)posar);
			dtmp = PyArray_DIMS((PyArrayObject*)posar);
			if(ndim != 2 || dtmp[1] != 3) {
				PyErr_SetString(PyExc_TypeError, "pos must be a 2D array with shape (npart, 3)");
				return -1;
			}
			self->cmesh.dim = dtmp[1];
			self->cmesh.npart = dtmp[0];
			tmp = self->pos;
			self->pos = posar;
			Py_INCREF(posar);
			Py_DECREF(tmp);
			self->cmesh.pos = (psi_rvec*)PyArray_DATA((PyArrayObject*)self->pos);
		}
		else {
			PyErr_SetString(PyExc_TypeError, "pos must be a numpy array");
			return -1;
		}
		if(PyArray_Check(connar)) {

			// elemtype based on shape of connectivity 
			if(PyArray_TYPE(connar) != NPY_INT32) {
				PyErr_SetString(PyExc_TypeError, "connectivity must be an array of ints");
				return -1;
			}
			ndim = PyArray_NDIM((PyArrayObject*)connar);
			dtmp = PyArray_DIMS((PyArrayObject*)connar);
			if(ndim != 2) {
				PyErr_SetString(PyExc_TypeError, "connectivity must be a 2D array with shape (nelem, nverts_per_elem)");
				return -1;
			}
			if(dtmp[1] != PSI_MESH_SIMPLEX && dtmp[1] != PSI_MESH_LINEAR && dtmp[1] != PSI_MESH_QUADRATIC) {
				PyErr_SetString(PyExc_TypeError, "connectivity must be a 2D array where nverts_per_elem is 4, 8, or 27.");
				return -1;
			}
			self->cmesh.elemtype = dtmp[1];
			self->cmesh.nelem = dtmp[0];
			tmp = self->connectivity;
			self->connectivity = connar;
			Py_INCREF(connar);
			Py_DECREF(tmp);
			self->cmesh.connectivity = (psi_int*)PyArray_DATA((PyArrayObject*)self->connectivity);
		}
		else {
			PyErr_SetString(PyExc_TypeError, "connectivity must be a numpy array");
			return -1;
		}

		// mass and velocity are not required, they will created if not provided 
		if(velar == NULL) {
			npdims[0] = self->cmesh.npart;
			npdims[1] = self->cmesh.dim;
			self->vel = PyArray_SimpleNew(2, npdims, NPY_FLOAT64);
			if(!self->vel) return -1;
			self->cmesh.vel = (psi_real*)PyArray_DATA((PyArrayObject*)self->vel);
			for(i = 0; i < self->cmesh.npart; ++i) {
				self->cmesh.vel[i].x = 0.0;
				self->cmesh.vel[i].y = 0.0;
				self->cmesh.vel[i].z = 0.0;
			}
		}
		else if(PyArray_Check(velar)) {
			if(PyArray_TYPE(velar) != NPY_FLOAT64) {
				PyErr_SetString(PyExc_TypeError, "if provided, vel must be an array of doubles");
				return -1;
			}
			ndim = PyArray_NDIM((PyArrayObject*)velar);
			dtmp = PyArray_DIMS((PyArrayObject*)velar);
			if(ndim != 2 || dtmp[1] != self->cmesh.dim || dtmp[0] != self->cmesh.npart) {
				PyErr_SetString(PyExc_TypeError, "if provided, vel must have the same shape as pos");
				return -1;
			}
			tmp = self->vel;
			self->vel = velar;
			Py_INCREF(velar);
			Py_DECREF(tmp);
			self->cmesh.vel = (psi_rvec*)PyArray_DATA((PyArrayObject*)self->vel);
		}
		else {
			PyErr_SetString(PyExc_TypeError, "if provided, vel must be a numpy array");
			return -1;
		}

		if(massar == NULL) {
			npdims[0] = self->cmesh.nelem;
			self->mass = PyArray_SimpleNew(1, npdims, NPY_FLOAT64);
			if(!self->mass) return -1;
			self->cmesh.mass = (psi_real*)PyArray_DATA((PyArrayObject*)self->mass);
			for(i = 0; i < self->cmesh.nelem; ++i)
				self->cmesh.mass[i] = 1.0;
		}
		else if(PyArray_Check(massar)) {
			if(PyArray_TYPE(massar) != NPY_FLOAT64) {
				PyErr_SetString(PyExc_TypeError, "if provided, mass must be an array of doubles");
				return -1;
			}
			ndim = PyArray_NDIM((PyArrayObject*)massar);
			dtmp = PyArray_DIMS((PyArrayObject*)massar);
			if(ndim != 1 || dtmp[0] != self->cmesh.nelem) {
				PyErr_SetString(PyExc_TypeError, "if provided, mass must be a 1D array of length nelem");
				return -1;
			}
			tmp = self->mass;
			self->mass = massar;
			Py_INCREF(massar);
			Py_DECREF(tmp);
			self->cmesh.mass = (psi_real*)PyArray_DATA((PyArrayObject*)self->mass);
		}
		else {
			PyErr_SetString(PyExc_TypeError, "if provided, mass must be a numpy array");
			return -1;
		}
	}
	else if(strcmp(cloader, "block") == 0) {

		if(PyArray_Check(posar)) {

			// get the dimensionality based on shape of pos 
			if(PyArray_TYPE(posar) != NPY_FLOAT64) {
				PyErr_SetString(PyExc_TypeError, "pos must be an array of doubles");
				return -1;
			}
			ndim = PyArray_NDIM((PyArrayObject*)posar);
			dtmp = PyArray_DIMS((PyArrayObject*)posar);
			if(ndim != 4 || dtmp[3] != 3) {
				PyErr_SetString(PyExc_TypeError, "pos must be a 4D array with shape (nx, ny, nz, 3)");
				return -1;
			}
			self->cmesh.dim = dtmp[3];
			self->cmesh.npart = dtmp[0]*dtmp[1]*dtmp[2];
			nside.i = dtmp[0];
			nside.j = dtmp[1];
			nside.k = dtmp[2];
			tmp = self->pos;
			self->pos = posar;
			Py_INCREF(posar);
			Py_DECREF(tmp);
			self->cmesh.pos = (psi_rvec*)PyArray_DATA((PyArrayObject*)self->pos);
		}
		else {
			PyErr_SetString(PyExc_TypeError, "pos must be a numpy array");
			return -1;
		}

		// set number of elements depending whether we wrap periodically
		if(self->cmesh.periodic) {
			for(ax = 0; ax < 3; ++ax)
				nelem.ijk[ax] = nside.ijk[ax];
		}
		else {
			for(ax = 0; ax < 3; ++ax)
				nelem.ijk[ax] = nside.ijk[ax] - 1;
		}

		// trilinear elements naturally
		self->cmesh.nelem = nelem.i*nelem.j*nelem.k;
		self->cmesh.elemtype = PSI_MESH_LINEAR;

		// mass and velocity are not required, they will created if not provided 
		if(velar == NULL) {
			npdims[0] = nside.i;
			npdims[1] = nside.j;
			npdims[2] = nside.k;
			npdims[3] = self->cmesh.dim;
			self->vel = PyArray_SimpleNew(4, npdims, NPY_FLOAT64);
			if(!self->vel) return -1;
			self->cmesh.vel = (psi_real*)PyArray_DATA((PyArrayObject*)self->vel);
			for(i = 0; i < self->cmesh.npart; ++i) {
				self->cmesh.vel[i].x = 0.0;
				self->cmesh.vel[i].y = 0.0;
				self->cmesh.vel[i].z = 0.0;
			}
		}
		else if(PyArray_Check(velar)) {
			if(PyArray_TYPE(velar) != NPY_FLOAT64) {
				PyErr_SetString(PyExc_TypeError, "if provided, vel must be an array of doubles");
				return -1;
			}
			ndim = PyArray_NDIM((PyArrayObject*)velar);
			dtmp = PyArray_DIMS((PyArrayObject*)velar);
			if(ndim != 4 || dtmp[3] != self->cmesh.dim || 
					dtmp[0] != nside.i || dtmp[1] != nside.j || dtmp[2] != nside.k) {
				PyErr_SetString(PyExc_TypeError, "if provided, vel must have the same shape as pos");
				return -1;
			}
			tmp = self->vel;
			self->vel = velar;
			Py_INCREF(velar);
			Py_DECREF(tmp);
			self->cmesh.vel = (psi_rvec*)PyArray_DATA((PyArrayObject*)self->vel);
		}
		else {
			PyErr_SetString(PyExc_TypeError, "if provided, vel must be a numpy array");
			return -1;
		}

		if(massar == NULL) {
			npdims[0] = nelem.i;
			npdims[1] = nelem.j;
			npdims[2] = nelem.k;
			self->mass = PyArray_SimpleNew(3, npdims, NPY_FLOAT64);
			if(!self->mass) return -1;
			self->cmesh.mass = (psi_real*)PyArray_DATA((PyArrayObject*)self->mass);
			for(i = 0; i < self->cmesh.nelem; ++i)
				self->cmesh.mass[i] = 1.0;
		}
		else if(PyArray_Check(massar)) {
			if(PyArray_TYPE(massar) != NPY_FLOAT64) {
				PyErr_SetString(PyExc_TypeError, "if provided, mass must be an array of doubles");
				return -1;
			}
			ndim = PyArray_NDIM((PyArrayObject*)massar);
			dtmp = PyArray_DIMS((PyArrayObject*)massar);
			if(ndim != 3 || dtmp[0] != nelem.i || dtmp[1] != nelem.j || dtmp[2] != nelem.k) {
				PyErr_SetString(PyExc_TypeError, "if provided, mass must be a 3D array of shape (nelemx, nelemy, nelemz)");
				return -1;
			}
			tmp = self->mass;
			self->mass = massar;
			Py_INCREF(massar);
			Py_DECREF(tmp);
			self->cmesh.mass = (psi_real*)PyArray_DATA((PyArrayObject*)self->mass);
		}
		else {
			PyErr_SetString(PyExc_TypeError, "if provided, mass must be a numpy array");
			return -1;
		}

		// now build the mesh connectivity
		npdims[0] = self->cmesh.nelem;
		npdims[1] = self->cmesh.elemtype;
		self->connectivity = PyArray_SimpleNew(2, npdims, NPY_INT32);
		if(!self->connectivity) return -1;
		self->cmesh.connectivity = (psi_int*) PyArray_DATA((PyArrayObject*)self->connectivity); 
		psi_int i, j, k, ii, jj, kk, locind, vertind;
		psi_int elemind;
		for(i = 0; i < nelem.i; ++i)
		for(j = 0; j < nelem.j; ++j)
		for(k = 0; k < nelem.k; ++k) {
			elemind = nelem.j*nelem.k*i + nelem.k*j + k;
			for(ii = 0; ii < 2; ++ii)
			for(jj = 0; jj < 2; ++jj)
			for(kk = 0; kk < 2; ++kk) {
				locind = 4*ii + 2*jj + kk;
				vertind = nside.j*nside.k*((i+ii)%nside.i) 
					+ nside.k*((j+jj)%nside.j) + ((k+kk)%nside.k);
				self->cmesh.connectivity[8*elemind+locind] = vertind;
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
		if(!peek_gadget2(&self->cmesh, cfile)) {
			PyErr_Format(PyExc_ValueError, "Bad file name: %s", cfile);
			return -1;
		}

		// wrap the self->cmesh data in a cubical Numpy array
		npdims[0] = self->cmesh.npart;
		npdims[1] = self->cmesh.dim;
		self->pos = PyArray_SimpleNew(2, npdims, NPY_FLOAT64);
		if(!self->pos) return -1;
		self->vel = PyArray_SimpleNew(2, npdims, NPY_FLOAT64);
		if(!self->vel) return -1;
		self->mass = PyArray_SimpleNew(1, npdims, NPY_FLOAT64);
		if(!self->mass) return -1;
		npdims[0] = self->cmesh.nelem;
		npdims[1] = self->cmesh.elemtype; 
		self->connectivity = PyArray_SimpleNew(2, npdims, NPY_INT32);
		if(!self->connectivity) return -1;
		self->cmesh.pos = PyArray_DATA((PyArrayObject*)self->pos);
		self->cmesh.vel = PyArray_DATA((PyArrayObject*)self->vel);
		self->cmesh.mass = PyArray_DATA((PyArrayObject*)self->mass);
		self->cmesh.connectivity = PyArray_DATA((PyArrayObject*)self->connectivity);

		// load the data into the numpy buffers
		if(strcmp(cloader, "gadget2") == 0) {
			if(!load_gadget2(&self->cmesh, cfile)) {
				PyErr_Format(PyExc_ValueError, "Could not read %s", cfile);
				return -1;
			}
		} 
		else if(strcmp(cloader, "gevolution") == 0) {
			if(!load_gevolution(&self->cmesh, cfile)) {
				PyErr_Format(PyExc_ValueError, "Could not read %s", cfile);
				return -1;
			}
		} 
	}
	else {
		PyErr_SetString(PyExc_ValueError, "Invalid mesh loader.");
		return -1;
	} 

	if(self->cmesh.periodic) {
		tmp = self->wrapbox;
		self->wrapbox = Py_BuildValue("((ddd)(ddd))", 
			self->cmesh.box[0].x, self->cmesh.box[0].y, self->cmesh.box[0].z, 
			self->cmesh.box[1].x, self->cmesh.box[1].y, self->cmesh.box[1].z);
		Py_DECREF(tmp);
	}

    return 0;
}

static PyObject* Mesh_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    Mesh *self;
    self = (Mesh*)type->tp_alloc(type, 0);
	if (self != NULL) {
		memset(&self->cmesh, 0, sizeof(psi_mesh));
		self->pos = Py_None; 
		Py_INCREF(Py_None);
		self->vel = Py_None; 
		Py_INCREF(Py_None);
		self->mass = Py_None; 
		Py_INCREF(Py_None);
		self->connectivity = Py_None; 
		Py_INCREF(Py_None);
		self->wrapbox = Py_None; 
		Py_INCREF(Py_None);
	}
    return (PyObject*)self;
}


static void Mesh_dealloc(Mesh* self) {
    Py_XDECREF(self->pos);
    Py_XDECREF(self->vel);
    Py_XDECREF(self->mass);
    Py_XDECREF(self->connectivity);
    Py_XDECREF(self->wrapbox);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyMemberDef Mesh_members[] = {
    {"pos", T_OBJECT_EX, offsetof(Mesh, pos), 0, "pos"},
    {"vel", T_OBJECT_EX, offsetof(Mesh, vel), 0, "vel"},
    {"mass", T_OBJECT_EX, offsetof(Mesh, mass), 0, "mass"},
    {"connectivity", T_OBJECT_EX, offsetof(Mesh, connectivity), 0, "connectivity"},
    {"wrapbox", T_OBJECT_EX, offsetof(Mesh, wrapbox), 0, "wrapbox"},
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

typedef struct {
    PyObject_HEAD
	psi_rtree ctree; 
} RStarTree;

static PyObject* RStarTree_query(RStarTree *self, PyObject *args, PyObject *kwds) {

	psi_int e;
	PyObject* box;
	psi_rvec qbox[2];
	npy_intp npdims[6];
	static char *kwlist[] = {"box", NULL};

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O", kwlist, &box))
		return NULL;

	// parse the query box, should be a tuple
	if(!PyArg_ParseTuple(box, "(ddd)(ddd)", &qbox[0].x, &qbox[0].y, &qbox[0].z, 
				&qbox[1].x, &qbox[1].y, &qbox[1].z)) {
		PyErr_SetString(PyExc_TypeError, "Query box must be a tuple containing two tuples of 3 values each");
		return NULL;
	}

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
	if(!pyinds) return NULL;
	PyArray_ENABLEFLAGS((PyArrayObject*)pyinds, NPY_ARRAY_OWNDATA);
	return pyinds;

}


static int RStarTree_init(RStarTree *self, PyObject *args, PyObject *kwds) {

	psi_int ax, dim, order, nside, npix, init_cap;
    psi_rvec tbox[2];
	npy_intp npdims[6];
	Mesh* mesh = NULL;
	static char *kwlist[] = {"mesh", "initial_capacity", NULL};

	// parse arguments differently depending on what the grid type is
	// certain grid types must correspond to certain arg patterns
	init_cap = 1024;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O!|i", kwlist, &MeshType, &mesh, &init_cap))
			return -1;

	// if mesh is periodic, use its box, otherwise use an arbitrarily large box
	if(mesh->cmesh.periodic) {
    	for(ax = 0; ax < 3; ++ax) {
    	    tbox[0].xyz[ax] = mesh->cmesh.box[0].xyz[ax];
    	    tbox[1].xyz[ax] = mesh->cmesh.box[1].xyz[ax];
    	}
	}
	else {
    	for(ax = 0; ax < 3; ++ax) {
    	    tbox[0].xyz[ax] = -1.0e30;
    	    tbox[1].xyz[ax] = +1.0e30;
    	}
	}

	// TODO: pass inital capacity
	psi_rtree_from_mesh(&self->ctree, &mesh->cmesh, &tbox);

    return 0;
}

static PyObject* RStarTree_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    RStarTree *self;
    self = (RStarTree*)type->tp_alloc(type, 0);
	if(self) {
		memset(&self->ctree, 0, sizeof(psi_rtree));
	}
    return (PyObject*)self;
}

static void RStarTree_dealloc(RStarTree* self) {
	psi_rtree_destroy(&self->ctree);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyMemberDef RStarTree_members[] = {
    {NULL}  /* Sentinel */
};

static PyMethodDef RStarTree_methods[] = {
	{"query", (PyCFunction)RStarTree_query, METH_KEYWORDS | METH_VARARGS,
	 "Query the RStarTree."},
    {NULL}  /* Sentinel */
};

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
//     PSI module functions 
/////////////////////////////////////////////////////////////



static PyObject *PSI_voxels(PyObject *self, PyObject *args, PyObject* kwds) {

	psi_int mode, maxlvl;
	psi_real reftol;
	char* modestr;
	Mesh* mesh;
	Grid* grid;	
	npy_intp nverts;
	static char *kwlist[] = {"grid", "mesh", "mode", "refine_tolerance", "refine_max_lvl", NULL};

	// defaults
	modestr = "density";
	reftol = 1.0;
	maxlvl = 0;

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!|sdi", kwlist, &GridType, &grid, &MeshType, &mesh, &modestr, &reftol, &maxlvl)) {
		return NULL;
	}

	// extract C pointers and such
	if(grid->cgrid.type != PSI_GRID_CART) {
		PyErr_SetString(PyExc_ValueError, "Only cartesian grids are supported by PSI.voxels()");
		return NULL;
	}
	
	if(mesh->cmesh.dim != grid->cgrid.dim) {
		PyErr_SetString(PyExc_ValueError, "Input mesh and grid must have the same number of dimensions");
		return NULL;
	}

	// parse the deposit mode
	if(strcmp(modestr, "density") == 0)
		mode = PSI_MODE_DENSITY;
	else if(strcmp(modestr, "annihilation") == 0)
		mode = PSI_MODE_ANNIHILATION;
	else {
		PyErr_Format(PyExc_ValueError, "Invalid mode for PSI.voxels(): %s", modestr);
		return NULL;
	}

	// call the C function
	// TODO: allow passing of a user-made r* tree
	psi_voxels(&grid->cgrid, &mesh->cmesh, NULL, mode, reftol, maxlvl);

	Py_RETURN_NONE;
}



static PyObject *PSI_VDF(PyObject *self, PyObject *args, PyObject* kwds) {

	psi_int maxlvl, nsamp;
	psi_real reftol;
	psi_rvec samppos;
	psi_rtree* ctree;
	psi_real* rhoout;
	psi_rvec* velout;
	Mesh* mesh = NULL;
	RStarTree* rst = NULL;
	npy_intp npdims[3];
	npy_intp nverts;
	static char *kwlist[] = {"mesh", "sample_pos", "tree", "refine_tolerance", "refine_max_lvl", NULL};

	reftol = 1.0;
	maxlvl = 0;

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O!(ddd)|O!di", kwlist, &MeshType, &mesh, &samppos.x, &samppos.y, &samppos.z, 
				&RStarTreeType, &rst, &reftol, &maxlvl))
		return NULL;

	// call the C function
	if(rst) ctree = &rst->ctree;
	psi_sample_vdf(samppos, &mesh->cmesh, ctree, 
		&rhoout, &velout, &nsamp, reftol, maxlvl);

	// make the return arrays from the c buffers
	// force numpy to own the memory
	npdims[0] = nsamp;
	npdims[1] = 3;
	PyObject* pyrho = PyArray_SimpleNewFromData(1, npdims, NPY_FLOAT64, rhoout);
	if(!pyrho) return NULL;
	PyObject* pyvel = PyArray_SimpleNewFromData(2, npdims, NPY_FLOAT64, velout);
	if(!pyvel) return NULL;
	PyArray_ENABLEFLAGS((PyArrayObject*)pyrho, NPY_ARRAY_OWNDATA);
	PyArray_ENABLEFLAGS((PyArrayObject*)pyvel, NPY_ARRAY_OWNDATA);
	return Py_BuildValue("OO", pyrho, pyvel);
}




/////////////////////////////////////////////////////////////
//     PSI module init 
/////////////////////////////////////////////////////////////

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif

static PyMethodDef module_methods[] = {
    //{"skymap", (PyCFunction)PSI_skymap, METH_KEYWORDS | METH_VARARGS, "Makes a skymap"},
    //{"beamtrace", (PyCFunction)PSI_beamtrace, METH_KEYWORDS | METH_VARARGS, "beamtrace"},
    {"voxels", (PyCFunction)PSI_voxels, METH_KEYWORDS | METH_VARARGS, "Voxelizes"},
    //{"phi", (PyCFunction)PSI_phi, METH_KEYWORDS | METH_VARARGS, "phi"},
    {"VDF", (PyCFunction)PSI_VDF, METH_KEYWORDS | METH_VARARGS, "VDF"},
    //{"crossStreams", (PyCFunction)PSI_crossStreams, METH_KEYWORDS | METH_VARARGS, "crossStreams"},
    //{"powerSpectrum", (PyCFunction)PSI_powerSpectrum, METH_KEYWORDS | METH_VARARGS, "powerSpectrum"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef PSI = {
    PyModuleDef_HEAD_INIT,
    "PSI",     /* m_name */
    "PSI, the phase space intersector",  /* m_doc */
    -1,                  /* m_size */
    module_methods    /* m_methods */
};

PyMODINIT_FUNC PyInit_PSI(void)
{
    PyObject* m;
    m = PyModule_Create(&PSI);
    if(!m) return NULL;

	// import numpy functionality
	_import_array();

	// add the mesh type
	if(PyType_Ready(&MeshType) < 0) return NULL;
	Py_INCREF(&MeshType);
	PyModule_AddObject(m, "Mesh", (PyObject*)&MeshType);


	// add the grid type
	if(PyType_Ready(&GridType) < 0) return NULL;
	Py_INCREF(&GridType);
	PyModule_AddObject(m, "Grid", (PyObject*)&GridType);

	// add the rtree type
	if(PyType_Ready(&RStarTreeType) < 0) return NULL;
	Py_INCREF(&RStarTreeType);
	PyModule_AddObject(m, "RStarTree", (PyObject*)&RStarTreeType);

	return m;
}

#endif // PYMODULE




#if 0


static PyObject *PSI_skymap(PyObject *self, PyObject *args, PyObject* kwds) {

	psi_int bstep, mode;
	Mesh* mesh;
	Grid* grid;	
	//npy_intp* npdims;
	npy_intp nverts;
	static char *kwlist[] = {"grid", "mesh", "bstep", "mode", NULL};

	// parse PyArgs
	bstep = 1;
	mode = PSI_SKYMAP_RHO_LINEAR;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!|ii", kwlist, &GridType, &grid, &MeshType, &mesh, &bstep, &mode))
		return NULL;

	// make C structs, check them
	if(grid->cgrid.type != PSI_GRID_HPRING || grid->cgrid.dim != 3) 
		return NULL;
	if(mesh->cmesh.dim != grid->cgrid.dim)
		return NULL;

	// run the skymap
	psi_printf("----Making hpring skymap, mode = %d-----\n", mode);
	psi_skymap(&grid->cgrid, &mesh->cmesh, bstep, mode);
	Py_RETURN_NONE;
}



static PyObject *PSI_phi(PyObject *self, PyObject *args, PyObject* kwds) {

	psi_int ax;
	npy_intp npdims[6];
	psi_grid cgrid;
	psi_real Gn;
	Grid* grid;	
	static char *kwlist[] = {"grid", "Gn", NULL};

#ifdef HAVE_FFTW
	Gn = 1.0;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|d", kwlist, &grid, &Gn))
		return NULL;

	cgrid = grid->cgrid;
	if(cgrid.type != PSI_GRID_CART) 
		return NULL;

	for(ax = 0; ax < 3; ++ax)
		npdims[ax] = cgrid.n.ijk[ax];

	// the return array
	PyObject* retar = PyArray_ZEROS(cgrid.dim, npdims, NPY_FLOAT64, 0);
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

	cgrid = grid->cgrid;
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
	PyObject* pspec = PyArray_ZEROS(1, npdims, NPY_FLOAT64, 0);
	PyObject* kk = PyArray_ZEROS(1, npdims, NPY_FLOAT64, 0);
	psi_do_power_spectrum(&cgrid, PyArray_DATA(pspec), PyArray_DATA(kk), TWO_PI/lmin, nbins);
	return Py_BuildValue("OO", pspec, kk);
#else
	PyErr_SetString(PyExc_RuntimeWarning, 
			"PSI was compiled without FFTW. PSI.powerSpectrum() does nothing.");
	return NULL;
#endif

}



static PyObject* Grid_getCellGeometry(Grid *self, PyObject *args, PyObject *kwds) {

	psi_rvec boundary[64], center;
	psi_real vol;
	psi_int ax, slen, v, i, allcells, bstep;
	psi_dvec grind;
	psi_grid cgrid;
	PyObject *seq=NULL, *pyind=NULL;
	PyObject *bverts=NULL, *cverts=NULL, *pyvol=NULL;
	npy_intp npdims[6];
	double *cvdata, *bvdata, *pvdata;
	static char *kwlist[] = {"cell", "bstep", NULL};

	// TODO: update this function

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
			cverts = PyArray_SimpleNew(2, npdims, NPY_FLOAT64);
			npdims[1] = 4;
			npdims[2] = 3;
			bverts = PyArray_SimpleNew(3, npdims, NPY_FLOAT64);
			pyvol = PyArray_SimpleNew(1, npdims, NPY_FLOAT64);
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
				if(slen <= 0 || !PyLong_Check(PySequence_Fast_GET_ITEM(seq, 0))) 
					return NULL;
			}
	
			// create the output array in the correct shape (flat for hp)
			// get direct pointers to the data
			// iterate over all integer indices in the sequence
			// or all cells if cells=None
			npdims[0] = slen;
			npdims[1] = 3;
			cverts = PyArray_SimpleNew(2, npdims, NPY_FLOAT64);
			npdims[1] = 4*bstep;
			npdims[2] = 3;
			bverts = PyArray_SimpleNew(3, npdims, NPY_FLOAT64);
			pyvol = PyArray_SimpleNew(1, npdims, NPY_FLOAT64);
	
			cvdata = PyArray_DATA((PyArrayObject*)cverts);
			bvdata = PyArray_DATA((PyArrayObject*)bverts);
			pvdata = PyArray_DATA((PyArrayObject*)pyvol);
			for(i = 0; i < slen; ++i) {
				if(allcells) grind.i = i;
				else grind.i = PyLong_AsLong(PySequence_Fast_GET_ITEM(seq, i));
	
				// retrieve the grid geometry from the C lib
				// set the array elements 
				if(!psi_grid_get_cell_geometry(&cgrid, grind, bstep, boundary, &center, &vol)) {
					psi_printf("bad geometry\n");
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
	Py_RETURN_NONE;
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
	npy_intp npdims[6];
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
    PyObject* pyvel = PyArray_SimpleNew(1, npdims, NPY_FLOAT64);
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
    PyObject* pycov = PyArray_SimpleNew(2, npdims, NPY_FLOAT64);
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



/////////////////////////////////////////////////////////////
//     Metric
/////////////////////////////////////////////////////////////


typedef struct {
    PyObject_HEAD
	PyObject* type; // string representation of the grid type 
	PyObject* params; // tuple of floats -- parameters for analytic metrics 
	PyObject* n; // number of cells in each dimension
	PyObject* box; // box size 
	PyObject* phi; // numpy array - potential 
	PyObject* gradphi; // numpy array - gradient of the potential 
} Metric;

static void PSI_Metric2metric(Metric* metric, psi_metric* cmetric) {


	// parse a simple C struct from a Python object 
	psi_int ax;
	PyObject* seq;
	char* cstr;


	// clear everything first
	memset(cmetric, 0, sizeof(psi_metric));
	cmetric->snapnum = -1;
   
	// get the enum metric type
	cstr = PyUnicode_AsUTF8(metric->type);
	psi_printf("Metric type = %s\n", cstr);

	if(strcmp(cstr, "minkowski") == 0) {
		psi_printf("Chose the Minkowski metric!\n");
		cmetric->type = PSI_METRIC_MINKOWSKI; 
	}
	else if(strcmp(cstr, "flrw") == 0) {
		psi_printf("Chose the FLRW metric!\n");
		cmetric->type = PSI_METRIC_FLRW; 
	}
	else if(strcmp(cstr, "kerr") == 0) {
		psi_printf("Chose the KERR metric!\n");
		cmetric->type = PSI_METRIC_KERR; 
	}
	else {
		psi_printf("Invalid metric\n");
		return NULL;
	}
}



static int Metric_init(Metric *self, PyObject *args, PyObject *kwds) {

	psi_printf("Init metric\n");

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
	type = PyUnicode_AsUTF8(PyDict_GetItemString(kwds, "type"));
	//if(strcmp(type, "flrw") == 0) {
		////self->type = 
		//psi_printf("Flrw selected\n");		
	//}

	if(strcmp(type, "minkowski") == 0 &&
			PyArg_ParseTupleAndKeywords(args, kwds, "S", minkkw, &pytype)) {
		self->type = pytype;
		psi_printf("Minkowski selected, %s\n", PyUnicode_AsUTF8(pytype));		
	}
	if(strcmp(type, "kerr") == 0 &&
			PyArg_ParseTupleAndKeywords(args, kwds, "S", kerrkw, &pytype)) {
		self->type = pytype;
		psi_printf("Kerr selected, %s\n", PyUnicode_AsUTF8(pytype));		
	}
	else if(strcmp(type, "flrw") == 0 &&
			PyArg_ParseTupleAndKeywords(args, kwds, "S|OS", flrwkw, &pytype, &params, &filepat)) {
		self->type = pytype;
		psi_printf("FLRW selected, %s\n", PyUnicode_AsUTF8(pytype));		
	}
	else {
	
		psi_printf("Bad metric!!!!\n");
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
	PyDict_SetItemString(self->fields, "m", PyArray_ZEROS(dim, npdims, NPY_FLOAT64, 0));

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


static PyObject *PSI_beamtrace(PyObject *self, PyObject *args, PyObject* kwds) {

	psi_grid cgrid;
	psi_mesh cmesh;
	psi_metric cmetric;
	psi_int bstep, mode;
	Metric* metric;
	Grid* grid;	
	psi_rvec obspos, obsvel;
	psi_dvec gdim;
	npy_intp npdims[6];
	npy_intp nverts;
	static char *kwlist[] = {"grid", "metric", "obspos", "obsvel", NULL};

	setbuf(stdout, NULL);

	psi_printf("---- GR beamtracing... -----\n");

	bstep = 1;
	mode = PSI_SKYMAP_RHO_LINEAR;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OO|(ddd)(ddd)", kwlist, &grid, &metric, &obspos.x, &obspos.y, &obspos.z, &obsvel.x, &obsvel.y, &obsvel.z))
		return NULL;

	cgrid = grid->cgrid;
	if(cgrid.type != PSI_GRID_HPRING || cgrid.dim != 3) 
		return NULL;

	psi_printf("metric2metric\n");

	PSI_Metric2metric(metric, &cmetric);

	psi_printf("metric2metric done\n");

	psi_printf("grid.nside = %d, metric type = %d, obspos = %f %f %f\n", cgrid.n.j, cmetric.type, obspos.x, obspos.y, obspos.z);




	npdims[0] = cgrid.n.k;
	npdims[1] = 100024; 
	npdims[2] = 4; 
	gdim.i = npdims[0];
	gdim.j = npdims[1];
	gdim.k = npdims[2];
	PyObject* rayinfo = PyArray_SimpleNew(3, npdims, NPY_FLOAT64);
	psi_real* infar = PyArray_DATA((PyArrayObject*)rayinfo);

	psi_beamtrace(&cgrid, &cmesh, bstep, &cmetric, infar, gdim);

   /* Do your stuff here. */
   //Py_RETURN_NONE;
   return rayinfo;
}




#endif


