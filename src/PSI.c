/*
 *
 * PSI.c 
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

// Python interface only
#ifdef PYMODULE 

#include "psi.h"
#include <numpy/arrayobject.h>
#include "structmember.h"
#include "grid.h"
#include "mesh.h"
#include "rtree.h"

// forward declarations for functions in cpsi.c

void psi_voxels(psi_grid* grid, psi_mesh* mesh, psi_rtree* rtree, psi_int mode, psi_real reftol, psi_int max_ref_lvl);
void psi_sample_vdf(psi_rvec samppos, psi_mesh* mesh, psi_rtree* rtree, 
		psi_real** rhoout, psi_rvec** velout, psi_int* nsamp , psi_real reftol, psi_int max_ref_lvl);

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

	psi_int ax, dim, i, len;
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
	PyObject *tmp = NULL, *posar = NULL, *velar = NULL, *massar = NULL, *box = NULL, *connar = NULL;
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
			if(PyArray_TYPE((PyArrayObject*)posar) != NPY_FLOAT64) {
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
			if(PyArray_TYPE((PyArrayObject*)connar) != NPY_INT32) {
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
			self->cmesh.vel = (psi_rvec*)PyArray_DATA((PyArrayObject*)self->vel);
			for(i = 0; i < self->cmesh.npart; ++i) {
				self->cmesh.vel[i].x = 0.0;
				self->cmesh.vel[i].y = 0.0;
				self->cmesh.vel[i].z = 0.0;
			}
		}
		else if(PyArray_Check(velar)) {
			if(PyArray_TYPE((PyArrayObject*)velar) != NPY_FLOAT64) {
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
			if(PyArray_TYPE((PyArrayObject*)massar) != NPY_FLOAT64) {
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
			if(PyArray_TYPE((PyArrayObject*)posar) != NPY_FLOAT64) {
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
			self->cmesh.vel = (psi_rvec*)PyArray_DATA((PyArrayObject*)self->vel);
			for(i = 0; i < self->cmesh.npart; ++i) {
				self->cmesh.vel[i].x = 0.0;
				self->cmesh.vel[i].y = 0.0;
				self->cmesh.vel[i].z = 0.0;
			}
		}
		else if(PyArray_Check(velar)) {
			if(PyArray_TYPE((PyArrayObject*)velar) != NPY_FLOAT64) {
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
			if(PyArray_TYPE((PyArrayObject*)massar) != NPY_FLOAT64) {
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
		self->cmesh.pos = (psi_rvec*)PyArray_DATA((PyArrayObject*)self->pos);
		self->cmesh.vel = (psi_rvec*)PyArray_DATA((PyArrayObject*)self->vel);
		self->cmesh.mass = (psi_real*)PyArray_DATA((PyArrayObject*)self->mass);
		self->cmesh.connectivity = (psi_int*)PyArray_DATA((PyArrayObject*)self->connectivity);

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

	psi_int ax, init_cap;
    psi_rvec tbox[2];
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
	psi_rtree_from_mesh(&self->ctree, &mesh->cmesh, &tbox[0]);

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
	psi_real* rhoout;
	psi_rvec* velout;
	Mesh* mesh = NULL;
	RStarTree* rst = NULL;
	psi_rtree* ctree = NULL;
	npy_intp npdims[3];
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
    {"voxels", (PyCFunction)PSI_voxels, METH_KEYWORDS | METH_VARARGS, "Voxelizes"},
    {"VDF", (PyCFunction)PSI_VDF, METH_KEYWORDS | METH_VARARGS, "VDF"},
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


