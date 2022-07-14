# Ionizer [![License](https://img.shields.io/github/license/sinaatalay/Ionizer.svg)](https://github.com/sinaatalay/Ionizer/blob/main/LICENSE)
An ion thruster simulation tool. It's a work in progress. Please come to revisit later.

Currently, **Ionizer** solves 3D [Poisson's equation](https://en.wikipedia.org/wiki/Poisson%27s_equation) for the potential field, and this is how it looks:

<p align="center">
  <img src="https://github.com/sinaatalay/Ionizer/blob/main/figures/preview.png?raw=true">
</p>

## Installation

>These installation instructions assume that [Windows](https://www.microsoft.com/en-us/windows/), [Git](https://git-scm.com/), [Visual Studio 2022](https://devcenter.heroku.com/articles/heroku-cli), and [Vulkan SDK](https://vulkan.lunarg.com/) are installed.

 1. Clone the repository with `git clone --recursive https://github.com/sinaatalay/Ionizer`. If the repository was cloned non-recursively previously, use `git submodule update --init --recursive` to clone the necessary submodules.

 2. Run [SetupVS2022.bat](https://github.com/sinaatalay/Ionizer/blob/master/scripts/SetupVS2022.bat) file found in `scripts` folder to generate Visual Studio solution and project files.
