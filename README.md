# ilemt-dual-rate-ukf

## Build steps for WSL on Ubuntu
1. Clone the repository
2. Install g++, cmake:
```
sudo apt-get install g++
sudo apt-get install cmake
```
3. Install CMake Tools extension for VSCode
4. Run CMake: Configure using CMake Tools
5. Install Eigen by downloading the eigen folder from http://eigen.tuxfamily.org/index.php?title=Main_Page#Download and moving the Eigen folder to `/usr/local/include`
6. Download and install zlib 
```
sudo apt install zlib1g-dev
```
7. Download, build, and install matio with zlib:
```
git clone git://git.code.sf.net/p/matio/matio
cd matio
git submodule update --init  # for datasets used in unit tests
./autogen.sh
./configure --with-zlib=/usr/include
make
make check
sudo make install
```
If `./autogen.sh` fails you may have to run `sudo apt-get install libtool m4 automake`

8. Run CMake: Build using CMake Tools
You may have to set `export LD_LIBRARY_PATH=/usr/local/lib`

9. Run the program:
```
build/ilemt-dual-rate-ukf
```

## Usage
You must use the `trace_to_2d.m` script in robmacl/ilemt to output new traces to a `*_couplings.mat` file that this implementation can read. Place the file in the directory you run the executable from. Change the `std::string trace` variable in `main()` of `dual_rate_ukf.cpp` to the name of the trace (without "_couplings.mat"). The biases and poses will be written to .mat files that can be opened by MATLAB for visualization.
