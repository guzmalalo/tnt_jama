

Welcome to TNT documentation!
=======================================

The Template Numerical Toolkit (TNT) is a collection of interfaces and reference implementations of numerical objects useful for scientific computing in C++. The toolkit defines interfaces for basic data structures, such as multidimensional arrays and sparse matrices, commonly used in numerical applications. The goal of this package is to provide reusable software components that address many of the portability and maintenance problems with C++ codes.

TNT provides a distinction between **interfaces** and **implementations** of TNT components. For example, there is a TNT interface for two-dimensional arrays which describes how individual elements are accessed and how certain information, such as the array dimensions, can be used in algorithms; however, there can be several implementations of such an interface: one that uses expression templates, or one that uses BLAS kernels, or another that is instrumented to provide debugging information. By specifying only the interface, applications codes may utilize such algorithms, while giving library developers the greatest flexibility in employing optimization or portability strategies.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   api/library_root


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
