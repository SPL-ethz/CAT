      ___           ___                   
     /\__\         /\  \                  
    /:/  /        /::\  \         ___     
   /:/  /        /:/\:\  \       /\__\    
  /:/  /  ___   /:/ /::\  \     /:/  /    
 /:/__/  /\__\ /:/_/:/\:\__\   /:/__/     
 \:\  \ /:/  / \:\/:/  \/__/  /::\  \     
  \:\  /:/  /   \::/__/      /:/\:\  \    
   \:\/:/  /     \:\  \      \/__\:\  \   
    \::/  /       \:\__\          \:\__\  
     \/__/         \/__/           \/__/  

Crystallization Analysis Toolbox
a SPLIfA project

brought to you by:

Martin Iggland, Separation Process Laboratory, miggland@ipe.mavt.ethz.ch
Dave Ochsenbein, Automatic Control Laboratory, ochsenbein@control.ee.ethz.ch

CAT is an open-source software designed to solve population balance equations as they typically arise in particulate processes and to analyze the results. 

Numerical Methods currently supported
- Moving Pivot
- Central Difference
- High Resolution

Features supported by all solvers
- Nucleation (homogeneous/heterogeneous)
- Growth (size dependent/independent)
- Dissolution (only size independent verified)
- Antisolvent and Temperature profiles in form of anonymous functions or piecewise-linear functions
- Arbitrary grid sizing

Features that are currently planned for the future
- Ostwald ripening
- Lattice-Boltzman method
- 1D agglomeration

** Installation

See Install.m for installation instructions

** Getting started

After installation, type
>> help CAT
or have a look at the demos to get an idea of how to use the toolbox.

** Bug reports

Please submit bugs via the Github page:
https://github.com/miggland/PBEToolbox



This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 3 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

