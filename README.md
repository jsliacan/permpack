# permpack

In many ways a "barebones" implementation of Razborov's flag algebras framework for permutations in
[SageMath](http://sagemath.org). The usage resembles [Flagmatic's](http://flagmatic.org) (written by Emil Vaughan). I worked with Flagmatic before (see [Flagmatic](https://github.com/jsliacan/flagmatic.git)
repository on my Github).

set up
------

Assuming you have downloaded or cloned this repository, navigate to its `permpack/pkg` subdirectory.

    $ cd path/to/permpack/pkg/
  
Then run Sage

    $ sage

and import everything once in Sage.

    sage: from permpack.all import *
  
For later stages you might need to install CSDP solver to work with
your Sage.

    $ sage -i csdp

You can also find a CSDP 6.1.1 Mac binaries in `permpack/utils` subdirectory.

notes
-----
1. This package assumes that the user knows how to use it: no input checking and only basic docs.

2. Tested on

  * Mac OSX El Capitan with Sage-6.4 and Sage-7.3beta3
  * Ubuntu 14.04 with Sage-7.2

3. When using Sage-6.x, x<10, on El Capitan, you need to disable *System Integrity Protection* by pressing `cmd+R` on boot, then launching Terminal and typing `csrutil disable`. Reboot. To enable, repeat the process and type `csrutil enable`. You should know what this does before doing it.

4. Permpack stores flag products for problems that it encountered before. Then loads them from the file next time. This is sometimes useful (faster for big problems) and sometimes annoying (creates 60MB files). 
