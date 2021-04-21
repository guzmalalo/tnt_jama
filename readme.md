## Template Numerical Toolkit (TNT)

This is a fork of the Template Numerical Toolkit [(TNT)](https://math.nist.gov/tnt/index.html)

The Template Numerical Toolkit is a collection of interfaces and reference implementations of numerical objects useful for scientific computing in C++. The toolkit defines interfaces for basic data structures, such as multidimensional arrays and sparse matrices, commonly used in numerical applications. The goal of this package is to provide reusable software components that address many of the portability and maintenance problems with C++ codes.

TNT provides a distinction between interfaces and implementations of TNT components. For example, there is a TNT interface for two-dimensional arrays which describes how individual elements are accessed and how certain information, such as the array dimensions, can be used in algorithms; however, there can be several implementations of such an interface: one that uses expression templates, or one that uses BLAS kernels, or another that is instrumented to provide debugging information. By specifying only the interface, applications codes may utilize such algorithms, while giving library developers the greatest flexibility in employing optimization or portability strategies.

<p xmlns:dct="http://purl.org/dc/terms/">
<a rel="license" href="http://creativecommons.org/publicdomain/mark/1.0/">
<img src="https://licensebuttons.net/p/mark/1.0/80x15.png"
     style="border-style: none;" alt="Public Domain Mark" />
</a>
<br />
This work (<span property="dct:title">Template Numerical Toolkit</span>, by <a href="https://math.nist.gov/~RPozo/" rel="dct:creator"><span property="dct:title">Roldan Pozo: Original Author)</span></a>, and modified by <a href="https://www.linkedin.com/in/eguzmanmal/" rel="dct:creator"><span property="dct:title"> Eduardo Guzman: Build, Test, Doc, CI, Hosting, extended interfaces, Maintenance</span></a>), identified by <a href="https://gitlab.com/egm_foss/tnt" rel="dct:publisher">https://gitlab.com/egm_foss/tnt</a>, is free of known copyright restrictions.
</p>
