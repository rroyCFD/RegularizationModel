## Application : RegularizationModel: Regularizes the non-linear convection term 

### Author:
- Rajib Roy
- University of Wyoming
- rroy@uwyo.edu, roy.rajib@live.com

### Description
the non-linear convection term is explicitly estimated. THe velocity (U) and mass-flux (phi) are filtered using a LES type filter specified in the __regularization__ sub-dictionary in the system/fvSolution case file. In this class, the regularization is 4th order accurate; implemented as:

    Uf = filter(U);     UPrime = (U - Uf)
    phif       = filter(phi); phiPrime = (phi - phif)
    C(phi, U)  = fvc::div(phi, U)
    C2(phi, U) = filter(C(phif, Uf))
    C4(phi, U) = C(phif, Uf) + filter(C(phif, UPrime)) + filter(C(phiPrime, Uf))
    C6(phi, U) = C(phif, Uf) + C(phif, UPrime))+ C(phiPrime, Uf) + filter(C(phiPrime, UPrime))

    A6(phi, U) =  C(phi, U) - residual(C(phiPrime, UPrime))

Sample regularization dictionary in fvSolution file:

    regularization
    {
        // Verstappen regularization
        regOrder  C6; // C4; // C2;
        filter    polyLaplace;

        // A6: compact formulation of C6
        regOrder  A6;
        filter    polyLaplaceResidual;

        // epsilon 2
        d1        0.16666667;
        d2        0.00416667;

        // epsilon 3
        d1        0.375;
        d2        0.0375;
    }

    Verstappen, R. (2008).
    On restraining the production of small scales of motion in a turbulent channel flow.
    Computers & Fluids, 37(7), 887–897. https://doi.org/10.1016/J.COMPFLUID.2007.01.013


### Disclaiimer:

This application is built based on [OpenFOAM version-6](https://openfoam.org/release/6/). Please read the _About OpenFOAM_ section to learn more on OpenFOAM.

The application is free to use. The author neither provide any warranty nor shall be liable for any damage incurred from this application.



#### About OpenFOAM

OpenFOAM is the leading free, open source software for computational fluid dynamics (CFD), owned by the OpenFOAM Foundation and distributed exclusively under the [General Public Licence (GPL)](http://www.gnu.org/copyleft/gpl.html). The GPL gives users the freedom to modify and redistribute the software and a guarantee of continued free use, within the terms of the licence. To learn more visit [https://openfoam.org/](https://openfoam.org/)
