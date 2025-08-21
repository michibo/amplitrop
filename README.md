# amplitrop

This repository contains accompanying computer code to the paper ['Tropicalized quantum field theory and global tropical sampling'](https://arxiv.org/abs/2508.14263) by Michael Borinsky.

Installation
--

To download, compile and run the code the [git](https://git-scm.com/) program and [gcc](https://gcc.gnu.org/), a C++ compiler are needed. Your locally installed g++ compiler should support [OpenMP](https://www.openmp.org/).

To download the repository, run

``git clone https://github.com/michibo/amplitrop.git``

This code internally uses the packages [Eigen](https://eigen.tuxfamily.org/), [cxxopts](https://github.com/jarro2783/cxxopts), [xsum](https://github.com/yafshar/xsum), and [JSON](https://github.com/nlohmann/json).

To compile the code, run

``make``

Usage
--

To evaluate the massive phi^3 theory 3-point function in D=3 at 10 loops, run

``./amplitrop -k3 -D3 -n3 -L10 -N1000000``

The parameter -k sets the k in phi^k, -D sets the dimension, -n sets the multiplicity, -L sets the loop number, and -N fixes the number of sample points to evaluate.

Analogously, to evaluate the primitive contribution to the phi^4 theory beta function at 10 loops, run

``./amplitrop -k4 -D4 -n4 --prim -L10 -N1000000``
