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
   
   The FORTRAN code uses a finite element method and solves only 2-D PDEs. The finite element method implemented here uses unstructured triangular meshes in the spatial discretization of the problem domain. In the current code, only Maxwell's equations for 2-D magnetotelluric (MT) problems are solved.
   
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
   Specify the modelling parameters (e.g., MT frequency). Use the templates of input files to make sure all required parameters are provided in these files. The details of input files are discussed in :ref:`label-input-files`. Note that the code requires that **one special file** is passed to the (compiled) program as an argument during the running; that input file itself contains information of all other input files:
   
   .. code-block::

      >>> FE2D -iFWD_input.dat
   
   In the above, "FE2D" is the name of the compiled FORTRAN program, "FWD_input.dat" is the name of the special input file, and "-i" is one of the flags of the program and is used to specify the name of the input file.
   
Compilation
-----------


.. _label-mesh-generation:

Mesh generation
---------------

The required mesh files for the code are .node, .ele, .edge and .neigh files generated from `Triangle <https://www.cs.cmu.edu/~quake/triangle.html>`_. The file names and their locations are defined in :ref:`label-mesh-parameter`. For the requirements of the FORTRAN code, each region of the problem domain needs to be assigned an integer marker when generating the mesh files; this is typically done in a .poly file that is required for generating the mesh. For details of regional marker requirements in this code, see :ref:`label-mesh-parameter`.


.. _label-input-files:

Input files
-----------

There are a few input files required for running the forward solver. 

``FWD_input.dat``
^^^^^^^^^^^^^^^^^^^

This is the special file passed to the program when running the code. A template of the file is given as:

   .. literalinclude:: /MT2D_3layer (example 1)/FWD_input.dat

As we can see, a general rule is that any texts after "#" are considered as comments. The information in the file is listed in a "**Name Value**" format: the "**Name**" is the parameter name recognized by the code and is not allowed to change; and "**Value**" is the value of that parameter. There should be at least a space between "**Name**" and "**Value**".

This special file in the above example has five parameters (all parameter names are case-sensitive):

   * ``input_path``: the path for all other input files defined here, for example, the file ``modelling_parameter.in`` defined as the value to the parameter   ``input_modelling`` should be found using this path.
   * ``input_modelling``: its value is a file name; that file contains general modelling parameters and their values.
   * ``input_EM_generic``: its value is a file name; that file contains general parameters related to EM problems (MT is a special case of EM).
   * ``input_mesh``: its value is a file name; that file contains information of mesh files used by the code.
   .. * ``input_MT``: its value is a file name; that file contains all necessary modelling parameters that are specific to MT.

All paths and file names should be double or single quoted.

``modelling_parameter.in``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An example template of the file is:

   .. literalinclude:: /MT2D_3layer (example 1)/modelling_parameter.in

For MT problems, always choose ``DataType`` as "CM" (i.e, complex-valued data type). For the matrix equation solver, two Krylov iterative methods, GMRES and BCGSTAB (Bi-CG stabilized), are provided here. If the iterative solver parameter ``Iter_solver`` is "f" (i.e., false), then a direct solver will be used.

``EM_generic_parameter.in``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An example template is:

   .. literalinclude:: /MT2D_3layer (example 1)/EM_generic_parameter.in

In this file, the regular frequency (in Hz) and a conductivity list file name need to be defined. The file name can include an absolute or relative path. For the content of the conductivity list file, see :ref:`label-cond-file`.



.. _label-mesh-parameter:

``mesh_parameter.in``
^^^^^^^^^^^^^^^^^^^^^^^

An example template is:

   .. literalinclude:: /MT2D_3layer (example 1)/mesh_parameter.in

The definitions of the above parameters are:
   * ``nregion``: integer parameter; number of homogeneous regions of conductivity in the mesh.
   * ``regionMarkFile``: the file name of lists of regional markers in the mesh. Each homogenous region is assigned an integer marker in order to distinguish itself from other regions. The order of listing the markers can be specified by the user when generating the mesh. See :ref:`label-marker-file` for an example of this file.
   * ``meshfilepath``: string parameter; the path of all mesh files (see :ref:`label-mesh-generation`).
   * ``basefilename`` string parameter; the "base" part of the file names of mesh files. By default, all mesh files generated using the above-mentioned external program should share this part in their names. For example, if the node file is named "foo.node", then ``basefilename`` should be "foo". 


.. _label-cond-file:

``list_regional_cond.txt``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An example template is:

   .. literalinclude:: /MT2D_3layer (example 1)/list_regional_cond.txt
   
In this file, the conductivities (i.e., material properties, unit: S/m) of different uniform regions are defined. The first line must begin with "L", which is followed by the number of actual regions in the mesh/problem domain. In this example, there are 3 regions in the mesh. Subsequent lines list the conductivity of each region. The order of this list must be the same order that is used to define different regions (see :ref:`label-mesh-parameter`) during meshing.



.. _label-marker-file:

``list_regionMark.txt``
^^^^^^^^^^^^^^^^^^^^^^^^^^^

An example template is:

   .. literalinclude:: /MT2D_3layer (example 1)/list_regionMark.txt
   
   
Output files
------------

Indices and tables (coming up)
==============================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
