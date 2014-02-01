To build the documentation from the *.rst files, make sure that you have Sphinx (http://sphinx-doc.org/tutorial.html) installed, and then run::

    make html

The documentation will be built in ./_build/html

*****************************

Note that the conf.py file has been modified to autodocument examples from their README files. To do this, put the example in its own subdirectory in ../examples/ with a README.rst file, and then also list it in conf.py.
