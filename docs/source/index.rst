.. FE2D_DomainDecop documentation master file, created by
   sphinx-quickstart on Fri Jun 18 01:42:35 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to FE2D_DomainDecop's documentation!
============================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. topic:: Overview
   
   This is the documentation page for the FORTRAN code developed by Jianbo Long in support of research in domain decomposition methods at Prof. Ronald Haynes's group. For any questions, please contact me at jl7037@mun.ca.

   Currently, the FORTRAN code only acts as a "forward solver" in a domain decomposition algorithm, i.e., the code only computes the solution to a given partial differential equation (PDE) over the global domain or a local subdomain of the problem. For computing the global solution, an initial boundary condition is required; whereas for a subdomain solution, the boundary values of the subdomain are extracted from the global solution and the initial boundary values.
   
   The FORTRAN code uses a finite element method and solves only 2-D PDEs. The finite element method implemented here uses unstructured triangular meshes in the spatial discretization of the problem domain. In this code, only Maxwell's equations for 2-D magnetotelluric (MT) problems are solved.
   
.. contents:: 
    :depth: 3
    
Get started
===========

To run the code, follow these steps:

**Step 1**:
   Compile the code using the given Makefile. Make sure the required libraries are installed first. For details about how to compile the code, see `Compilation`_ below.
   
**Step 2**:
   Generate the mesh for the given 2-D problem. To do this, use the external program `Triangle <https://www.cs.cmu.edu/~quake/triangle.html>`_ to generate the mesh files, which are provided to the code as part of input files.
   
**Step 3**:
   Specify the modelling parameters (e.g., MT frequency). Use the templates of input files to make sure all required parameters are provided in these files. The details of input files are discussed in `Input files`_. Note that the code requires that **one special file** is passed to the (compiled) program as an argument during the running; that input file itself contains information of all other input files:
   
   .. code-block::

      >>> FE2D -iFWD_input.dat
   
   In the above, "FE2D" is the name of the compiled FORTRAN program, "FWD_input.dat" is the name of the special input file, and "-i" is one of the flags of the program.
   
Compilation
-----------

Input files
-----------

There are a few input files required for running the forward solvers. 

File FWD_input.dat
^^^^^^^^^^^^^^^^^^

This is the special file passed to the program when running the code. A template is given as:
.. code-block::
   # --put all input files for the forward modelling in one place. This file is the direct input file for modelling program.
   # -- "Name Value" format; "value" has the file name which can include the path info.

   input_path          "/media/jack/NewVolume/3DEM_Jianbo_test_data/MT2D_3layer/"    # path info for all file defined here (if they are in one place)  
   input_modelling     "modelling_parameter.in"       # general modelling parameters
   input_EM_generic    "EM_generic_parameter.in"
   input_mesh          "mesh_parameter.in"
   input_MT            "MT_parameter.in"      # provide this if doing MT modelling

Output files
------------

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
