![Arpra](./doc/arpra_logo.svg)

Copyright 2016-2022 James Paul Turner.

This file is part of the Arpra library.

The Arpra library is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The Arpra library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the Arpra library. If not, see <http://www.gnu.org/licenses/>.

For any copyright year range specified as YYYY-ZZZZ in this package,
note that the range specifies every single year in that closed interval.


Introduction
============

Arpra is a C library for (Ar)bitrary-(p)recision (r)ange (a)nalysis
of IEEE-754 floating-point computations, based on GNU MPFR. The main
use-case of Arpra is to maintain computed upper and lower bounds of
numerical error for all variables, at all times, throughout a
computation. Arpra uses mixed trimmed interval/affine arithmetic with
deviation term reduction to accomplish this.

Affine arithmetic is a variant of interval arithmetic which accounts for
variable correlations. As such, it does not suffer from the so-called
'dependency problem', where intervals grow overly large due to lack of
consideration for variable correlations. The problem is described further
at <https://en.wikipedia.org/wiki/Interval_arithmetic#Dependency_problem>.
By combining the results of interval arithmetic and affine arithmetic,
one avoids both the dependency problem of interval arithmetic, and the
nonlinear function overshoot/undershoot problem of affine arithmetic.

Arpra implements affine arithmetic using a GNU MPFR backend. MPFR is an
arbitrary-precision floating-point library, meaning the floating-point
MPFR variables can be of arbitrary precision. For more information,
refer to the MPFR project website at <http://www.mpfr.org/>. By
implementing affine arithmetic with an arbitrary-precision backend,
one is able, for example, to test how a change in numerical precision
or integration scheme affects local and global error during a long
numerical simulation, without the interval 'explosion' problem regular
interval arithmetic suffers from.

For further information on the implementation and features of Arpra,
such as range trimming and deviation term reduction, refer to the
original published article:

Turner, J. P., & Nowotny, T. (2021). Arpra: An Arbitrary Precision
Range Analysis Library. Frontiers in Neuroinformatics, 30.

https://doi.org/10.3389/fninf.2021.632729


Quickstart
==========

Arpra follows the familiar GNU/Linux software building paradigm. The
typical installation procedure consists of the following.

If installing from the Git source repository (i.e. not a dist tarball),
the configure script and other auxillary files need to be generated
by running the following command in the root of the repository:

    autoreconf -i -Wall

This generates the configuration files from the configure.ac and the
Makefile.am files (note that the GNU Autotools must be installed in
order to run autoreconf). Next run the configure, build and install
commands:

    ./configure
    make
    sudo make install

All installed Arpra files can be cleanly uninstalled from the system by
running the following command:

    sudo make uninstall

A suite of test programs can be executed with the following command:

    make check


Contribute
==========

All contributions (e.g. bug reports, feature requests, expert knowledge,
source code and documentation contributions) are gratefully received via
the issue tracker <https://github.com/arpra-project/arpra/issues> or
pull request.

The source code repository for Arpra is hosted at GitHub. Clone it using:

    git clone https://github.com/arpra-project/arpra
