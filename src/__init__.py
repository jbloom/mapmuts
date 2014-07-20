"""Package for mapping mutations for analysis of deep-sequencing libraries.

Python modules
----------------
The package consists of the following Python modules:

* main

* align

* io

* latex

* sequtils

* plot

* simulate

* stats

* bayesian

* dssp

* weblogo


C-extensions
-------------
Some of the key functions in the Python modules are also coded in
C-extensions for better performance. These C-extensions can be accessed
directly through the Python functions, using options such as
`use_calign` when you call the function, and so you typically would not need
to address these C-extensions directly. The extensions are:

* calign

* csequtils


Acknowledgements
--------------------

``mapmuts`` is written by Jesse Bloom.
"""

__version__ = '1.1'
