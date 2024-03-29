# Ionizer [![License](https://img.shields.io/github/license/sinaatalay/Ionizer.svg)](https://github.com/sinaatalay/Ionizer/blob/main/LICENSE)
An ion thruster simulation tool. It's a work in progress. Please come to revisit later.

It's a tricky simulation because the applicability of the continuum assumption is questionable since the operating pressures are pretty low. Therefore, the molecular theory has to be employed.

**Ionizer** is being written in C++, and its goal is to solve the 3D Poisson's equation to calculate the electrostatic field inside the thruster and move the ions accordingly with the well-established molecular approach, PIC-DSMC.

Currently, it solves 3D [Poisson's equation](https://en.wikipedia.org/wiki/Poisson%27s_equation) for the potential field, which is how it looks:

<p align="center">
  <img src="https://github.com/sinaatalay/Ionizer/blob/main/figures/preview.png?raw=true">
</p>

## Installation

>These installation instructions assume that [Windows](https://www.microsoft.com/en-us/windows/), [Git](https://git-scm.com/), [Visual Studio 2022](https://visualstudio.microsoft.com/vs/), and [Vulkan SDK](https://vulkan.lunarg.com/) are installed.

 1. Clone the repository with 
  ```
  git clone --recursive https://github.com/sinaatalay/Ionizer
  ```
  If the repository was cloned non-recursively previously, use 
  ```
  git submodule update --init --recursive
  ``` 
  to clone the necessary submodules.

 2. Run [scripts/SetupVS2022.bat](https://github.com/sinaatalay/Ionizer/blob/master/scripts/SetupVS2022.bat) file to generate Visual Studio 2022 solution and project files.
