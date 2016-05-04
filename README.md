#Soth test#

## Installation ##
1) Download compile and install the soth package,(any folder will do, this is _NOT_ a catkin package)

```
git clone --recursive  git@github.com:stack-of-tasks/soth.git
cd soth
mkdir build
cd build/
cmake ..
sudo make istall
```

2) Dowload and compile the expressionGraph (catkin)

https://github.com/eaertbel/expressiongraph

3) and this package (catkin)

compile the source,
and try to run `test1` or `test2`.
if you get 
```bash
$ ./build/soth_test/test2 
./build/soth_test/test2: error while loading shared libraries: libsoth.so.2.0.1-9-g07d3: cannot open shared object file: No such file or directory
```
you probably have to 

```
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
```




