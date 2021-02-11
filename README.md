## GAPFix Fractal

This is a tool for creating Mandelbrot and Julia sets, which differs from
most of the numerous other such tools in two points. The admittedly somewhat
unwieldy name stands for GPU Arbitrary Precision Fixed Point Fractals.
Which is to say, it runs on CUDA-capable GPUs, and it uses arbitrary
precision fixed point math, compiled on the fly to any desired precision.
The goal is to be able to zoom in deeper than when using normal 64 bit
floating point math, while still getting decent speed due to the power of
modern GPUs.

![screenshot](screens/screenshot.png)

At the moment, this is still very new and fairly experimental code. It's
slowly becoming a real application, but expect it to be rough around the
edges in this early stage.

## Features

GAPFix Fractal has the following main features:
- Up to 512 bits of precision currently
- Arbitrary precision kernels generated on the fly as needed
- Julia set previews
- Supersampling
- A few different formulas
- DEM black and white images for the standard formula
- Storing parameters with previews and using them for navigation

## Usage

Most of the user interface should be self-explanatory or easily
understandable for anyone who has played with fractal software before.

Left-click recenters and zooms. Control-click selects parameters for
Julia sets. The "Mixture" formula allows choosing an additional
parameter q, currently only from a few presets in the menu.

The "Store" button saves the current position for later, with a small
preview image.

## Requirements

There are some prerequisites for building and running this program:
- An nvidia GPU, let's say GTX 1060 and up.
- The Qt library for the GUI
- Linux (it should build and run on Windows once the build system is
  adapted to find CUDA there)
- A fairly large amount of memory, especially when playing with large
  resolutions and supersampling factors.

## Limitations

Among the things currently known to need improvement are:
- The DEM algorithm works better than in the initial release (thanks Claude!),
  but still has issues rendering Julia sets near the origin point.
- There are also overflow issues at higher power variants that lead
  to some banding in color gradients which are supposed to be smooth.
- The "Spider" formula doesn't produce quite the same images seen
  elsewhere on the web. Probably a math error somewhere.
- Color palettes are still somewhat haphazard and subject to change.

## Gallery

I've uploaded a [sample gallery](https://photos.app.goo.gl/fFZyEvVNFHzDMaFu5)
of images created with this program.

## Compiling

On Linux, make a build subdirectory, enter it, and run
```sh
  qmake ../src/frac.pro PREFIX=/where/you/want/to/install
```
followed by make and make install.

## License

GAPFix Fractal is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

GAPFix Fractal is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with GAPFix Fractal.  If not, see <http://www.gnu.org/licenses/>.

### Additional permission under GNU GPL version 3 section 7

_The source files of GAPFix Fractal have the following additional
permission, as allowed under GNU GPL version 3 section 7:_

If you modify this Program, or any covered work, by linking or
combining it with NVIDIA Corporation's libraries from the NVIDIA CUDA
Toolkit a (or a modified version of those libraries), containing parts
covered by the terms of the respective license agreement, the
licensors of this Program grant you additional permission to convey
the resulting work.
