# amplitrop


Installation
--

To download, compile and run the code, git and a C++ compiler are needed.

To download the repository run

``git clone https://github.com/michibo/amplitrop.git``

This code uses the submodules, [Eigen](https://eigen.tuxfamily.org/), [cxxopts](https://github.com/jarro2783/cxxopts), [xsum](https://github.com/yafshar/xsum), and [JSON](https://github.com/nlohmann/json). To download these submodules run

``cd amplitrop && git submodule update --init --recursive``

To compile the code run

``make``

To evaluate the massive phi^3 theory 3-point function in D=3 at 10 loops run

``./amplitrop -k3 -D3 -n3 -L10 -N1000000``

The parameter -k sets the k in phi^k, -D sets the dimension, -n sets the multiplicity, -L sets the loop number, and -N fixes the number of sample points to evaluate.

