## Application : RegularizationModel: Regularizes the non-linear convection term 

### Author:
- Rajib Roy
- University of Wyoming
- rroy@uwyo.edu, roy.rajib@live.com

### Description
The non-linear convection term is explicitly estimated. THe velocity (U) and mass-flux (phi) are filtered using a LES type filter specified in the __regularization__ sub-dictionary in the system/fvSolution case file. In this class, the regularization is 4th order accurate; implemented as:

    Uf = filter(U);     UPrime = (U - Uf)
    phif = filter(phi); phiPrime = (phi - phif)
    C(phi, U) = fvc::div(phi, U)
    C4(phi, U) = C(phif, Uf) + filter(C(phif, UPrime)) + filter(C(phiPrime, Uf))

Sample filter description:

    regularization
    {
        filter    simple;
    }


    regularization
    {
        filter    polyLaplace;
        d1        0.16666667;
        d2        0.00416667;
    }

Filters of the 2nd and 6th order would be added in the future release, if needed. More about the regularization filters are available in 

    Verstappen, R. (2008).
    On restraining the production of small scales of motion in a turbulent channel flow.
    Computers & Fluids, 37(7), 887–897. https://doi.org/10.1016/J.COMPFLUID.2007.01.013


### Disclaiimer:

This application is built based on [OpenFOAM version-6](https://openfoam.org/release/6/). Please read the _About OpenFOAM_ section to learn more on OpenFOAM.

The application is free to use. The author neither provide any warranty nor shall be liable for any damage incurred from this application.



#### About OpenFOAM

OpenFOAM is the leading free, open source software for computational fluid dynamics (CFD), owned by the OpenFOAM Foundation and distributed exclusively under the [General Public Licence (GPL)](http://www.gnu.org/copyleft/gpl.html). The GPL gives users the freedom to modify and redistribute the software and a guarantee of continued free use, within the terms of the licence. To learn more visit [https://openfoam.org/](https://openfoam.org/)
