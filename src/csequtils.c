// Implements fast C versions of some of the functions in the sequtils
// module for the mapmuts package. Written by Jesse Bloom.
//
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Documentation string for MeanQValueFromString function
static char MeanQValueFromString_docs[] = "Fast C function for computing average Q value from string.\n\nExactly mimics the MeanQValueFromString in the mapmuts.sequtils module.\nThe calling and return values are the same as for that function.\nThe calling variable is the Python string qstring.\nThe return variable is a Python float giving the mean Q-value.";

// Implementation of MeanQValueFromString function
static PyObject *MeanQValueFromString(PyObject *self, PyObject *args) {
    // Calling variables: qstring
    char *qstring;
    size_t i, lenqstring, qsum;
    if (! PyArg_ParseTuple(args, "s", &qstring)) {
        PyErr_SetString(PyExc_TypeError, "Invalid calling arguments to MeanQValueFromString.");
        return NULL;
    }
    lenqstring = strlen(qstring);
    if (! lenqstring) {
        PyErr_SetString(PyExc_ValueError, "Empty qstring");
        return NULL;
    }
    qsum = 0;
    for (i = 0; i < lenqstring; i++)
        qsum += ((long) qstring[i]) - 33;
    return PyFloat_FromDouble(((double) qsum) / ((double) lenqstring));
}

// Documentation string for ReverseComplement function
static char ReverseComplement_docs[] = "Fast C function for DNA sequence reverse complements.\n\nTakes a single calling argument, which must be a Python string.\nThis string should be composed exclusively of DNA nucleotide codes:\nA, T, C, G, a, t, c, g, N, n.\nThe returned variable is a new string of the same length in which the\nnucleotides are reverse-complemented. Case is preserved\nand N/n reverse complements to N/n.";

// Implementation of ReverseComplement
static PyObject *ReverseComplement(PyObject *self, PyObject *args) {
    // Calling variables: seq
    char *seq;
    size_t seqlen, i, seqlen1;
    PyObject *pyrc;
    // Parse the arguments.  
    if (! PyArg_ParseTuple(args, "s", &seq)) {
        PyErr_SetString(PyExc_TypeError, "Invalid calling arguments to ReverseComplement.");
        return NULL;
    }
    seqlen = strlen(seq);
    seqlen1 = seqlen - 1;
    char *rc = PyMem_New(char, seqlen + 1);
    if (rc == NULL) {
        return PyErr_NoMemory();
    }
    for (i = 0; i < seqlen; i++) {
        switch (seq[i]) {
            case 'A' : rc[seqlen1 - i] = 'T';
                       break;
            case 'a' : rc[seqlen1 - i] = 't';
                       break;
            case 'C' : rc[seqlen1 - i] = 'G';
                       break;
            case 'c' : rc[seqlen1 - i] = 'g';
                       break;
            case 'G' : rc[seqlen1 - i] = 'C';
                       break;
            case 'g' : rc[seqlen1 - i] = 'c';
                       break;
            case 'T' : rc[seqlen1 - i] = 'A';
                       break;
            case 't' : rc[seqlen1 - i] = 'a';
                       break;
            case 'N' : rc[seqlen1 - i] = 'N';
                       break;
            case 'n' : rc[seqlen1 - i] = 'n';
                       break;
            default : PyErr_SetString(PyExc_ValueError, "Invalid nucleotide code.");
                      return NULL;
        }
    }
    rc[seqlen] = '\0'; // string termination character
    pyrc = PyString_FromString(rc);
    PyMem_Del(rc);
    return pyrc;
}


// Module documentation string
static char csequtils_docs[] = "Fast implementations of some functions from mapmuts.sequtils.\n\nReverseComplement in this module mimics some of the same function from sequtils.\n\nMeanQValueFromString in this module mimics the same function from sequtils.";

// The module functions
static PyMethodDef csequtils_funcs[] = {
    {"ReverseComplement", (PyCFunction) ReverseComplement, METH_VARARGS, ReverseComplement_docs},
    {"MeanQValueFromString", (PyCFunction) MeanQValueFromString, METH_VARARGS, MeanQValueFromString_docs},
    {NULL}
};

// Initialize the module
void initcsequtils(void) {
    Py_InitModule3("csequtils", csequtils_funcs, csequtils_docs);
}
