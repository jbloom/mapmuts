==============================
Installation
==============================

`mapmuts`_ is written in `Python`_. In addition to the pure `Python`_ code, it also has a few C-extensions which require a C compiler (such as `gcc`_). 



Dependencies
--------------
Although the core of `mapmuts`_  package can be run just using the standard `Python`_ libraries, there are lots of useful features that require the following external packages:

* `mapmuts`_ has been tested on Linux and Mac OS X.

* `mapmuts`_ requires `Python`_. It has been tested with the following `Python`_ versions:

    - 2.6.8
    
    - 2.7.3

* The graphing performed by `mapmuts`_ requires `matplotlib`_. It is known to work with the following `matplotlib`_ versions:

    - 1.2.1

    - 1.3.1

* Automatic compilation of plots into PDF summaries requires `pdflatex`_. 

* `numpy`_ is required for some operations. `mapmuts`_ has been tested with the following `numpy`_ versions:

    - 1.6.1

    - 1.7.1

* The MCMC Bayesian inference of the amino-acid preferences (or enrichment ratios) requires `pymc`_. `mapmuts`_ has been tested with `pymc`_ version 2.3. It may **not** work with other versions of `pymc`_.

* A few operations require `scipy`_. The following versions of `scipy`_ have been tested:

    - 0.9.0
    
    - 0.12.0

* Creation of sequence logos requires installation of `weblogo`_ in a location that puts it in the current search path so that it can be run from the shell. This package has been tested with `weblogo`_ version 3.4. It is not known if it would work with other versions. Note that if `weblogo`_ crashes with a Ghostscript error, you may want to see `this bug report <https://code.google.com/p/weblogo/issues/detail?id=36>`_. Some of the functionality also requires installation of the `Python`_ `weblogo`_ library (``weblogolib``) in a location that can be imported into `Python`_.

* Overlay of RSA and SS information on the sequence logos requires the `pyPdf`_ package to be installed. This package has only been tested with `pyPdf`_ version 1.13 -- it is unknown whether it would work with ``pyPdf2``.

* The `rpy2` package is required for the :ref:`mapmuts_entropycomparison.py` script. This package has been tested with version 2.3.9.

* The `R`_ statistical package is required for the :ref:`mapmuts_entropycomparison.py` script. This package has been tested with version 3.0.3.


Installation 
-----------------------

To install `mapmuts`_, first download the source ZIP repository `on GitHub`_. After unzipping the file, run the following commands::

    cd mapmuts
    python setup.py build
    python setup.py test
    python setup.py install

The last command might need to be replaced with::

    sudo python setup.py install
    
if you want to install the package globally do not have privileges to the default global installation directory. Alternatively, if you just want to install the package locally for the current user, you can use::
    
    python setup.py install --user

The ``test`` command runs a variety of tests to check that the program appears to be working properly. This is optional, and the tests may take a few minutes to run. Running these tests is recommended.

These commands install the `Python`_ modules and also install several scripts, which provide the most convenient high-level interface into `mapmuts`_ package. There is also an ``./examples/`` subdirectory in the main `mapmuts`_ directory that contains some examples. Looking at these examples may be helpful for understanding the typical workflow.

.. include:: weblinks.txt
