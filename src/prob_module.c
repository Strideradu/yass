//
// Created by Nan on 7/24/2017.
//

/*  Example of wrapping cos function from math.h with the Python-C-API. */

#include <Python.h>
#include <math.h>
#include "proba.h"

/*  wrapped cosine function */
static PyObject* stat_bound_func(PyObject* self, PyObject* args)
{
    double pI;
    long int L;
    double alpha;

    double value;
    double answer;

    /*  parse the input, from python float to c double */
    if (!PyArg_ParseTuple(args, "dld", &pI, &L, &alpha))
        return NULL;
    /* if the above function returns -1, an appropriate Python exception will
     * have been set, and the function simply returns NULL
     */

    /* call cos from libm */
    answer = statistical_bound_of_randomwalk2(pI, L, alpha);

    /*  construct the output from cos, from c double to python float */
    return Py_BuildValue("l", answer);
}

/*  define functions in module */
static PyMethodDef ProbMethods[] =
        {
                {"stat_bound_func", stat_bound_func, METH_VARARGS, "statistical bound"},
                {NULL, NULL, 0, NULL}
        };

/* module initialization */
PyMODINIT_FUNC

initprob_module(void) {
    (void) Py_InitModule("prob_module", ProbMethods);
}