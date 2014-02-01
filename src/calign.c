// Implements fast C versions of some of the functions in the align
// module for the mapmuts package. Written by Jesse Bloom.
//
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

// Documentation string for ParsePairedAlignment function
static char ParsePairedAlignment_docs[] = "Fast C implementation of align.ParsePairedAlignment in the mapmuts package.\n\nCalling arguments the same as for align.ParsePairedAlignment.\nReturn tuple is the same as for align.ParsePairedAlignment.";

// Implementation of ParsePairedAlignment function
static PyObject *ParsePairedAlignment(PyObject *self, PyObject *args) {
    // Calling variables: aligntup, r1exclude, r2exclude
    PyObject *aligntup, *r1exclude, *r2exclude;
    long gstart, gend, r1start, r1end, r2start, lena, ig, ir1, ir2, i;
    char *gseq, *r1seq, *r2seq;
    int r1isrc;
    char nt1, nt2;
    // Parse the arguments
    if (! PyArg_ParseTuple(args, "O!O!O!", &PyTuple_Type, &aligntup, &PySet_Type, &r1exclude, &PySet_Type, &r2exclude)) {
        PyErr_SetString(PyExc_TypeError, "Invalid calling arguments to ParsePairedAlignment.");
        return NULL;
    }
    gstart = PyInt_AS_LONG(PyTuple_GET_ITEM(aligntup, 0));
    gend = PyInt_AS_LONG(PyTuple_GET_ITEM(aligntup, 1));
    gseq = PyString_AS_STRING(PyTuple_GET_ITEM(aligntup, 2));
    r1start = PyInt_AS_LONG(PyTuple_GET_ITEM(aligntup, 3));
    r1end = PyInt_AS_LONG(PyTuple_GET_ITEM(aligntup, 4));
    r1seq = PyString_AS_STRING(PyTuple_GET_ITEM(aligntup, 5));
    r2start = PyInt_AS_LONG(PyTuple_GET_ITEM(aligntup, 6));
    r2seq = PyString_AS_STRING(PyTuple_GET_ITEM(aligntup, 8));
    lena = gend - gstart + 1;
    if (r1start > r1end)
        r1isrc = 1;
    else
        r1isrc = 0;
    char *ntidentities = PyMem_New(char, (size_t) lena + 1);
    ig = gstart;
    ir1 = r1start;
    ir2 = r2start;
    for (i = 0; i < lena; i++) {
        nt1 = r1seq[i];
        nt2 = r2seq[i];
        if ((nt1 == 'N') || (PySet_Contains(r1exclude, PyInt_FromLong(ir1)))) {
            if ((nt2 == 'N') || (PySet_Contains(r2exclude, PyInt_FromLong(ir2))))
                ntidentities[i] = 'N';
            else {
                if (nt2 == '.')
                    ntidentities[i] = tolower(gseq[i]);
                else
                    ntidentities[i] = tolower(nt2);
            }
        }
        else {
            if ((nt2 == 'N') || (PySet_Contains(r2exclude, PyInt_FromLong(ir2)))) {
                if ((nt1 == 'N') || (PySet_Contains(r1exclude, PyInt_FromLong(ir1))))
                    ntidentities[i] = 'N';
                else {
                    if (nt1 == '.')
                        ntidentities[i] = tolower(gseq[i]);
                    else
                        ntidentities[i] = tolower(nt1);
                }
            }
            else {
                if (nt1 == nt2) {
                    if (nt1 == '.')
                        ntidentities[i] = gseq[i];
                    else
                        ntidentities[i] = nt1;
                }
                else
                    ntidentities[i] = 'N';
            }
        }
        ig++;
        if (r1isrc) {
            ir1--;
            ir2++;
        }
        else {
            ir1++;
            ir2--;
        }
    }
    ntidentities[lena] = '\0';
    PyObject *pyntidentities = PyString_FromString(ntidentities);
    PyMem_Del(ntidentities);
    return Py_BuildValue("(lN)", gstart, pyntidentities);
}

// Documentation string for ReadToGeneMismatches function
static char ReadToGeneMismatches_docs[] = "Fast C implementation of align.ReadToGeneMismatches in the mapmuts package.\n\nCalling arguments the same as for align.ReadToGeneMismatches.\nUpper and lower case letters are treated differently.\nThe return list is the same as for align.ReadToGeneMismatches.";

// Implementation of ReadToGeneMismatches function
static PyObject *ReadToGeneMismatches(PyObject *self, PyObject *args) {
    // Calling variables g, r, gi, maxmm, indexback, n_is_mismatch
    char *g, *r;
    long gi, maxmm, imax, i, imms;
    size_t lenr;
    int ib, nismm;
    PyObject *indexback, *pymms, *n_is_mismatch;
    // Parse the arguments
    if (! PyArg_ParseTuple(args, "ssllOO", &g, &r, &gi, &maxmm, &indexback, &n_is_mismatch)) {
        PyErr_SetString(PyExc_TypeError, "Invalid calling arguments to ReadToGeneMismatches.");
        return NULL;
    }
    ib = PyObject_IsTrue(indexback);
    nismm = PyObject_IsTrue(n_is_mismatch);
    lenr = strlen(r);
    imax = gi + lenr;
    if (imax > strlen(g)) 
        return PyBool_FromLong((long) 0);
    long *mms = PyMem_New(long, maxmm);
    if (mms == NULL) 
        return PyErr_NoMemory();
    imms = 0;
    for (i = gi; i < imax; i++) {
        if (r[i - gi] != g[i]) {
            if ((nismm) || ((r[i - gi] != 'N') && (g[i] != 'N') && (r[i - gi] != 'n') && (g[i] != 'n'))) {
                if (imms >= maxmm) {
                    PyMem_Free(mms);
                    return PyBool_FromLong((long) 0);
                }
                if (ib) {
                    mms[imms] = lenr - i + gi;
                    imms++;
                }
                else {
                    mms[imms] = i + 1 - gi;
                    imms++;
                }
            }
        }
    }
    pymms = PyList_New(imms);
    for (i = 0; i < imms; i++) 
        PyList_SET_ITEM(pymms, i, PyInt_FromLong(mms[i]));
    PyMem_Free(mms);
    return pymms;
}

// Documentation string for AddDots function
static char AddDots_docs[] = "Fast C implementation of align.AddDots in the mapmuts package.\n\nCalling arguments the same as for align.AddDots.\nUpper and lower case letters are treated differently.\nThe return string is the dotted version of s2, as for align.AddDots.";

// Implementation of AddDots
static PyObject *AddDots(PyObject *self, PyObject *args) {
    // Calling variables s1, s2
    PyObject *pys2dotted;
    char *s1, *s2;
    size_t slength, i;
    // Parse the arguments
    if (! PyArg_ParseTuple(args, "ss", &s1, &s2)) {
        PyErr_SetString(PyExc_TypeError, "Invalid calling arguments to AddDots.");
        return NULL;
    }
    slength = strlen(s1);
    if (slength != strlen(s2)) {
        PyErr_SetString(PyExc_ValueError, "AddDots: s1 and s2 differ in length.");
        return NULL;
    }
    char *s2dotted = PyMem_New(char, slength + 1);
    if (s2dotted == NULL) {
        return PyErr_NoMemory();
    }
    for (i = 0; i < slength; i++) {
        if (s1[i] == s2[i]) 
            s2dotted[i] = '.';
        else
            s2dotted[i] = s2[i];
    }
    s2dotted[slength] = '\0'; // string termination character
    pys2dotted = PyString_FromString(s2dotted);
    PyMem_Del(s2dotted);
    return pys2dotted;
}


// Documentation string for AlignReadToGene function
static char AlignReadToGene_docs[] = "Fast C implementation of align.AlignReadToGene in the mapmuts package.\n\nCalling arguments the same except that there is no upcase argument.\nAny conversion of reads to upper case must be done prior to calling.\nUpper and lower case letters are treated differently by this function.\n\nThe return tuple is the same as for align.AlignReadToGene";

// Implementation of AlignReadToGene
static PyObject *AlignReadToGene(PyObject *self, PyObject *args) {
    // Calling variables r, gene, gene_rc, maxm
    char *r, *gene, *gene_rc, *g;
    char c1, c2;
    long maxm, lenr, lengene, i, j, nm, isrc, brokeloop, imax;
    PyObject *n_is_mismatch;
    int nismm;
    // Parse the arguments
    if (! PyArg_ParseTuple(args, "ssslO", &r, &gene, &gene_rc, &maxm, &n_is_mismatch)) {
        PyErr_SetString(PyExc_TypeError, "Invalid calling arguments to AlignReadToGene.");
        return NULL;
    }
    nismm = PyObject_IsTrue(n_is_mismatch);
    lenr = strlen(r);
    lengene = strlen(gene);
    if (lengene != strlen(gene_rc)) {
        PyErr_SetString(PyExc_ValueError, "AlignReadToGene: gene and gene_rc differ in length");
        return NULL;
    }
    imax = lengene - lenr + 1;
    for (isrc = 0; isrc <= 1; isrc++) {
        if (isrc == 0)
            g = gene;
        else
            g = gene_rc;

        for (i = 0; i < imax; i++) {
            nm = 0;
            brokeloop = 0;
            for (j = 0; j < lenr; j++) {
                c1 = g[i + j];
                c2 = r[j];
                if (c1 != c2) {
                    if (nismm || (c1 != 'N' && c1 != 'n' && c2 != 'N' && c2 != 'n')) {
                        nm++;
                        if (nm > maxm) {
                            brokeloop = 1;
                            break;
                        }
                    }
                }
            }
            if (! brokeloop) {
                return Py_BuildValue("(lNl)", i, PyBool_FromLong(isrc), nm);
            }
        }
    }
    return PyBool_FromLong((long) 0);
}



// Documentation string for AlignReads function
static char AlignReads_docs[] = "Fast C implementation of align.AlignReads in the mapmuts package.\n\nCalling arguments the same except that there is no upcase argument.\nAny conversion of reads to upper case must be done prior to calling.\nUpper and lower case letters are treated differently by this function.\n\nThe return tuple is the same as for align.AlignReads";

// Implementation of AlignReads
static PyObject *AlignReads(PyObject *self, PyObject *args) {
    // Calling variables: r1, r2, a1, a2, minoverlap, maxrm, maxa1m, maxa2m, n_is_mismatch
    char *r1, *r2, *a1, *a2;
    char c1, c2;
    size_t lenr1, lenr2, lena1, lena2;
    long minoverlap, maxrm, maxa1m, maxa2m;
    long i2max, i1, i2, nrm, j, nma1, nma2, k, i1max, jmax, kmax, r1a, r2a, overlap, x1, x2;
    int brokeloop, nismm;
    PyObject *n_is_mismatch;
    // Parse the arguments.  
    if (! PyArg_ParseTuple( args, "ssssllllO", &r1, &r2, &a1, &a2, &minoverlap, &maxrm, &maxa1m, &maxa2m, &n_is_mismatch)) {
        PyErr_SetString(PyExc_TypeError, "Invalid calling arguments to AlignReads.");
        return NULL;
    }
    nismm = PyObject_IsTrue(n_is_mismatch);
    lenr1 = strlen(r1);
    lenr2 = strlen(r2);
    lena1 = strlen(a1);
    lena2 = strlen(a2);
    i1max = lenr1 - minoverlap + 1;
    i2max = lenr2 - minoverlap + 1;
    for (i1 = 0; i1 < i1max; i1++) {
        for (i2 = 0; i2 < i2max; i2++) {
            nrm = 0;
            jmax = lenr1 - i1;
            x1 = lenr2 - i2;
            if (x1 < jmax) {
                jmax = x1;
            }
            brokeloop = 0;
            for (j = 0; j < jmax; j++) {
                c1 = r1[i1 + j];
                c2 = r2[i2 + j];
                if (c1 != c2) {
                    if (nismm || ((c1 != 'N') && (c1 != 'n') && (c2 != 'N') && (c2 != 'n'))) {
                        nrm++;
                        if (nrm > maxrm) {
                            brokeloop = 1;
                            break;
                        }
                    }
                }
            }
            if (! brokeloop) {
                r1a = i1 + j;
                r2a = lenr2 - i2;
                overlap = j;
                nma1 = 0;
                kmax = lenr1 - r1a;
                if (lena1 < kmax) {
                    kmax = lena1;
                }
                for (k = 0; k < kmax; k++) {
                    c1 = r1[overlap + k];
                    c2 = a1[k];
                    if (c1 != c2) {
                        if (nismm || ((c1 != 'N') && (c1 != 'n') && (c2 != 'N') && (c2 != 'n'))) {
                            nma1++;
                            if (nma1 > maxa1m) {
                                brokeloop = 1;
                                break;
                            }
                        }
                    }
                }
                if (! brokeloop) {
                    nma2 = 0;
                    kmax = lenr2 - r2a;
                    if (lena2 < kmax) {
                        kmax = lena2;
                    }
                    x1 = lenr2 - overlap - 1;
                    x2 = lena2 - 1;
                    for (k = 0; k < kmax; k++) {
                        c1 = r2[x1 - k];
                        c2 = a2[x2 - k];
                        if (c1 != c2) {
                            if (nismm || ((c1 != 'N') && (c1 != 'n') && (c2 != 'N') && (c2 != 'n'))) {
                                nma2++;
                                if (nma2 > maxa2m) {
                                    brokeloop = 1;
                                    break;
                                }
                            }
                        }
                    }
                    if (! brokeloop) {
                        return Py_BuildValue("(NNlll)", Py_BuildValue("(ll)", i1, i1 + overlap), Py_BuildValue("(ll)", i2, i2 + overlap), nrm, r1a, r2a);
                    }
                }
            }
        }
        i2max = 1;
    }
    return PyBool_FromLong((long) 0);
}


// Module documentation string
static char calign_docs[] = "Fast implementations of some methods from align in mapmuts package.\n\nThe AlignReads function from this module mimics the same function from align.\n\nThe AlignReadToGene function from this module mimics the same function from align.\n\nThe AddDots function from this module mimics the same function from align.\n\nThe ReadToGeneMismatches function from this module mimics the same function from align.\n\nThe ParsePairedAlignment function from this module mimics the same function from align.";

// The module functions
static PyMethodDef calign_funcs[] = {
    {"AlignReads", (PyCFunction) AlignReads, METH_VARARGS, AlignReads_docs},
    {"AlignReadToGene", (PyCFunction) AlignReadToGene, METH_VARARGS, AlignReadToGene_docs},
    {"AddDots", (PyCFunction) AddDots, METH_VARARGS, AddDots_docs},
    {"ReadToGeneMismatches", (PyCFunction) ReadToGeneMismatches, METH_VARARGS, ReadToGeneMismatches_docs},
    {"ParsePairedAlignment", (PyCFunction) ParsePairedAlignment, METH_VARARGS, ParsePairedAlignment_docs},
    {NULL}
};

// Initialize the module
void initcalign(void) {
    Py_InitModule3("calign", calign_funcs, calign_docs);
}
