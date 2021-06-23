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
   
   The FORTRAN code uses a finite element method and solves only 2-D PDEs. The finite element method implemented here uses unstructured triangular meshes in the spatial discretization of the problem domain.
   
.. contents:: 
    :depth: 3
    
Get started
===========

**Step 1**:
   Compile the code using the given Makefile. Make sure the required libraries are installed first. For details about how to compile the code, see `Compilation`_.
   
Compilation
-----------

Input files
-----------



Output files
------------

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
