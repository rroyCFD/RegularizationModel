/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "RegularizationModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(RegularizationModel, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

surfaceScalarField Foam::RegularizationModel::surfaceFieldFilter
(const surfaceScalarField& ssf_)
{
    volVectorField vvf_(ssf_.name(), fvc::reconstruct(ssf_));
    vvf_ = filter_(vvf_);

    return(fvc::flux(vvf_));
}

// explicit convection function
inline volVectorField Foam::RegularizationModel::convOperator
(
    const surfaceScalarField& ssf_,
    const volVectorField& vvf_
)
{
    // cross product of the convection operation
    if(extpFilterFieldDivFree_)
    {
       return (fvc::div(ssf_, vvf_, "div(phi,U)"));
    }
    else
    {
        return (fvc::div(ssf_, vvf_, "div(phi,U)") - fvc::div(ssf_)*vvf_);
    }
}


void Foam::RegularizationModel::calcContinuityError (surfaceScalarField& ssf_)
{
    volScalarField contErr(fvc::div(ssf_));

    scalar sumLocalContErr(
        runTime_.deltaTValue()*
        mag(contErr)().weightedAverage(mesh_.V()).value());

    scalar globalContErr(
        runTime_.deltaTValue()*
        contErr.weightedAverage(mesh_.V()).value());

    Info<< "time step continuity errors : sum local = " << sumLocalContErr
        << ", global = " << globalContErr << endl;
}

// solve a poison pressure equation to make filtered flux divergence free
void Foam::RegularizationModel::setDivergenceFree
(
    surfaceScalarField& ssf_,
    volVectorField& vvf_
)
{
    Info << "continuity error before correction" << endl;
    calcContinuityError(ssf_);

    dimensionedScalar dt = runTime_.deltaT();

    // initialize pp to zero
    pp_.primitiveFieldRef() = 0.0;

    for (label nonOrth = 0; nonOrth <= nNonOrthCorr_; nonOrth++)
    {
        fvScalarMatrix ppEqn
        (
          fvm::laplacian(pp_, "laplacian(p)") == (fvc::div(ssf_)*((k_ + 0.5)/dt) )
        );
        // no correction to reference cell pressure
        ppEqn.setReference(pRefCell_, 0);

        if (nonOrth < nNonOrthCorr_)
        {
            ppEqn.solve(mesh_.solver(pp_.name()));
        }
        else
        {
            ppEqn.solve(mesh_.solver(pp_.name()+"Final")); // pp_.select(1)

            // Info << "Making " << ssf_.name() << " divergence free" << endl;
            ssf_ -= (ppEqn.flux() * (dt/(k_ + 0.5)) );
        }
    }

    Info << "continuity error after correction" << endl;
    calcContinuityError(ssf_);

    vvf_ -= (fvc::grad(pp_, "grad(p)") * (dt/(k_+ 0.5)) );
    vvf_.correctBoundaryConditions();

    Info << endl;

    return;
}


volVectorField Foam::RegularizationModel::getConvectionTerm
(
    const surfaceScalarField& phie_,
    const volVectorField& Ue_
)
{
    volVectorField convTerm_
    (
        IOobject
        (
            "convectionTerm",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("convTerm",dimAcceleration,vector::zero)
    );

    // regular projection
    if(! regOn_)
    {
        convTerm_ = convOperator(phie_, Ue_);

    }
    else // C4 regularization
    {
        // Filter velocity and flux using polynomial Laplace filter
        tmp<volVectorField> tUef_
        (
            new volVectorField
            (
                IOobject
                (
                    "Uf",
                    Ue_.instance(), // runTime_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                filter_(Ue_)
             )
        );
        volVectorField& Uef_ = tUef_.ref();
        Uef_.correctBoundaryConditions();

        surfaceScalarField phief_("phif", surfaceFieldFilter(phie_));

        // This approach add a little more divergeence error,
        // take more iterations to be divergence free > more time required
        //
        // phief is derived from filtered-extrapolated velocity
        // instead of filtering. Hence making phief divergence-free is essential
        // surfaceScalarField phief_("phif", fvc::flux(Uef_));

        // Make filtered flux (and correcponding velocity) divergence free
        if(extpFilterFieldDivFree_)
        {
            setDivergenceFree (phief_, Uef_);
        }

        if(regOrder_ == "C4" || regOrder_ == "c4")
        {
            // Fist term: (C(us_f) uc_f)
            convTerm_ = convOperator(phief_, Uef_);

            // Second term: Filt(C(us_f) uc')
            volVectorField Cint( "Cint", convOperator(phief_, (Ue_- Uef_)));
            Cint.correctBoundaryConditions();
            convTerm_ += filter_(Cint);

            // Third term: Filt(C(us') uc_f)
            Cint = convOperator((phie_- phief_), Uef_);
            Cint.correctBoundaryConditions();
            convTerm_ += filter_(Cint);
        }
        else if (regOrder_ == "C2"|| regOrder_ == "c2")
        {
            // Filter(C(us_f) uc_f)
            convTerm_ = filter_(convOperator(phief_, Uef_));
        }
        else if(regOrder_ == "C6"|| regOrder_ == "c6")
        {
            // Fist term: (C(us_f) uc_f)
            convTerm_ = convOperator(phief_, Uef_);

            // Second term: (C(us_f) uc')
            volVectorField Uprime("Uprime", (Ue_- Uef_));
            convTerm_ += convOperator(phief_, Uprime);

            // Third term: C(us', uc_f)
            surfaceScalarField phiPrime("phiPrime", (phie_- phief_));
            convTerm_ += convOperator(phiPrime, Uef_);

            // Fourth term: Filt(C(us', uc'))
            // volVectorField Cint = convOperator(phiPrime, Uprime);
            // Cint.correctBoundaryConditions();
            // convTerm_ += filter_(Cint);

            convTerm_ += filter_(convOperator(phiPrime, Uprime));
        }
        // Explicit regularization proposed by Bose et al.
        else if (regOrder_ == "Explicit"|| regOrder_ == "explicit")
        {
            convTerm_ = filter_(convOperator(phie_, Ue_));
        }
        else
        {
            FatalErrorIn("Regularization order") << regOrder_
                << " is not recognized!\n"
                << " Avaialble orders are C2, C4 C6, or Explicit"
                << abort(FatalError);
        }



        if(runTime_.outputTime())
        {
            // convTerm_.write();

            Ue_.write();
            Uef_.write();

            volVectorField Uprime("Uprime", (Ue_- Uef_));
            Uprime.write();

            volScalarField divPhif("divPhif", fvc::div(phief_));
            divPhif.write();
        }

        tUef_.clear();
        phief_.clear();
    }

    return convTerm_;
}

// add residual function header
#include "residual.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RegularizationModel::RegularizationModel
(
    volScalarField& pp,
    const label&  pRefCell,
    const scalar& pRefValue,
    const label nNonOrthCorr
)
:
    // Set the pointer to runTime
    runTime_(pp.time()),
    mesh_(pp.mesh()),

    regOn_(false),
    extpFilterFieldDivFree_(false),
    k_(0.5),

    // get regularization sub dictionary from fvSolution
    regDict_(mesh_.solutionDict().subDict("regularization")),

    // get regularization order; if available
    regOrder_(regDict_.lookupOrDefault<word>("regOrder", "C4")),

    // get LES filter specified in the regularization dictionary
    filterPtr_(LESfilter::New(mesh_, regDict_)),
    filter_(filterPtr_()),

    residualDict_(mesh_.solutionDict().subDict("regResidual")),
    residualOrder_(residualDict_.lookupOrDefault<word>("residualOrder", "C6")),
    resPtr_(LESfilter::New(mesh_, residualDict_)),
    residual_(resPtr_()),

    pp_(pp),

    pRefCell_(pRefCell),
    pRefValue_(pRefValue),

    nNonOrthCorr_(nNonOrthCorr)
{
    Info << "Regularization dictionary: " << regDict_ << endl;
    Info << "Regularization residual dictionary: " << residualDict_ << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RegularizationModel::~RegularizationModel()
{
    Info << "RegularizationModel Destructor" << endl;
}


// ************************************************************************* //
