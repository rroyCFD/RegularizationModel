/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::RegularizationModel

Description:
    Regularization of the non-linear convection term

SourceFiles
    RegularizationModel.C
\*---------------------------------------------------------------------------*/
#ifndef RegularizationModel_H
#define RegularizationModel_H

#include "fvCFD.H"
#include "IOdictionary.H"
#include "LESfilter.H"

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class RegularizationModel Declaration
\*---------------------------------------------------------------------------*/

class RegularizationModel
{
    // Private data

        //- Runtime pointer
        const Time& runTime_;;

        //- Mesh pointer
        const fvMesh& mesh_;

         // regularization On or projection?
        bool regOn_;

        // make extrapolated fields divergence free?
        bool extpFilterFieldDivFree_;

        // Adam-bashforth extrapolation coefficient
        scalar k_;

        // regularization dict: contaning regularization and filter types
        const dictionary regDict_;

        // regularization order specifier; default C4
        word regOrder_;

        // LES filter for U and phi
        autoPtr<LESfilter> filterPtr_;
        LESfilter& filter_;

        // const dictionary residualDict_;
        word residualOrder_;
        // autoPtr<LESfilter> resPtr_;
        LESfilter& residual_;


    // Pointer to solver fields and scalars
        // pressure-correction field pp pointer
        volScalarField& pp_;

        // Reference pressure parameters
        const label&  pRefCell_;
        const scalar& pRefValue_;

        // number of non-orthogonal correctors
        const label nNonOrthCorr_;

    // Private member functions
        void setDivergenceFree(surfaceScalarField& ssf_ );

        void setDivergenceFree(
            surfaceScalarField& ssf_,
            volVectorField& vvf_
        );

        inline volVectorField convOperator(
            const surfaceScalarField& ssf_,
            const volVectorField& vvf_
        );

        // void calcContinuityError (surfaceScalarField& ssf_);

public:

    //- Declare name of the class and its debug switch
    ClassName("RegularizationModel");


    // Constructors
    RegularizationModel
    (
        volScalarField& pp,
        const label&  pRefCell,
        const scalar& pRefValue,
        const label nNonOrthCorr
    );


    // Destructor
    virtual ~RegularizationModel();


    // Public Member functions

        // Set regulariozation on/off
        void setRegOn(bool regOn)
        {
            regOn_ = regOn;
        }

        void calcContinuityError (const surfaceScalarField& ssf_);
        void calcContinuityError (const volVectorField& vvf_);

        void setFilterFieldsDivFree(bool extpFilterFieldDivFree)
        {
            extpFilterFieldDivFree_ = extpFilterFieldDivFree;
        }

        volVectorField getConvectionTerm(
                                        const surfaceScalarField& phie_,
                                        const volVectorField& Ue_
                                        );

        surfaceScalarField surfaceFieldFilter(const surfaceScalarField& ssf_);

        volVectorField convectionResidual(
                                        const surfaceScalarField& phie_,
                                        const volVectorField& Ue_
                                        );

        surfaceScalarField surfaceFieldResidual(const surfaceScalarField& ssf_);

        // volVectorField SGSfilter(volVectorField& vvf)
        // {
        //     return filter_(vvf);
        // }

        // volVectorField SGSresidual(volVectorField& vvf)
        // {
        //     return residual_(vvf);
        // }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif