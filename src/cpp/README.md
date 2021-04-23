# paction_cpp (SPRUCE)

## Compilation instructions

### Dependencies

* [CMake](http://www.cmake.org/) (>=3.1)
* [Boost](http://www.boost.org) (>= 1.38)
* [LEMON](http://lemon.cs.elte.hu/trac/lemon) graph library (>= 1.3)

### Compilation

    $ mkdir build
    $ cd build
    $ cmake ..
    $ make

### Usage

```
./enumerate <CNA tree file> <CNA proportions file> <SNV tree file> <SNV proportions file> 0 <output prefix>
```
