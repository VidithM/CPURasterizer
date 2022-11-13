# 3D CPU rasterizer written for CSCE 441 at Texas A&M
## Features
* Can load variety of .obj files (found in /resources)
* Texture mapping and coloring (switch between modes using 't')
* Anti-aliasing (using nearest neighbor, bilinear interpolation, and mipmapping methods)

## How to build
* CMake is used to build the project. Install [here](https://cmake.org/).
* [GLFW](https://www.glfw.org/), [GLM](https://glm.g-truc.net/0.9.9/index.html), and [GLEW](http://glew.sourceforge.net/) are dependencies. Click on the respective links to download.
* After installing the above dependencies, set global environment variables ```GLFW_INCLUDE_DIR```, ```GLFW_DIR```, and ```GLEW_DIR``` to the respective install locations.
* Download this project and create a new directory in it called ```build```. 
* In ```build```, run CMake using the ```CMakeLists.txt``` in the parent directory.
* The result will depend on if you are using Windows or a Unix system. In the latter case, there should be an executable called ```CPURasterizer``` in ```build```, that you can run. 
