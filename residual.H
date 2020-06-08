surfaceScalarField Foam::RegularizationModel::surfaceFieldResidual
(const surfaceScalarField& ssf_)
{
    volVectorField vvf_(ssf_.name(), fvc::reconstruct(ssf_));
    vvf_ = residual_(vvf_);

    return(fvc::flux(vvf_));
}


volVectorField Foam::RegularizationModel::convectionResidual
(
    const surfaceScalarField& phie_,
    const volVectorField& Ue_
)
{
    tmp<volVectorField> tConvRes_
    (
        new volVectorField
        (
            IOobject
            (
                "convRes",
                Ue_.instance(), // runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("",dimAcceleration,vector::zero)
            // residual_(
            //     fvc::div(
            //         surfaceFieldResidual(phie_), residual_(Ue_),
            //         "div(phiRes,URes)"
            //     )
            // )
        )
    );
    volVectorField& convRes_ = tConvRes_.ref();

    if(regOn_)
    {
        // Resisual velocity and flux using polynomial Laplace filter
        tmp<volVectorField> tURes_
        (
            new volVectorField
            (
                IOobject
                (
                    "URes",
                    Ue_.instance(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                residual_(Ue_)
             )
        );
        volVectorField& URes_ = tURes_.ref();
        URes_.correctBoundaryConditions();

        surfaceScalarField phiRes_("phiRes", surfaceFieldResidual(phie_));

        if(residualOrder_ == "C6"|| residualOrder_ == "c6")
        {
            // Make filtered flux  divergence free;
            // velocity div. error is already high, correction not required
            if(extpFilterFieldDivFree_)
            {
                setDivergenceFree (phiRes_);
                convRes_ = residual_(fvc::div(phiRes_, URes_, "div(phi,U)"));
            }
            else
            {
                convRes_ = residual_(
                    fvc::div(phiRes_, URes_, "div(phi,U)")
                    // phiRes correction make small improvement, 
                    // but no-where close to conservation wo divFree corection
                    - URes_*fvc::div(phiRes_)
                );
            }
        }
        else
        {
            FatalErrorIn("Regularization order") << residualOrder_
                << " is not recognized!\n" << " Avaialble orders is C6"
                << abort(FatalError);
        }

        tURes_.clear();
        phiRes_.clear();
    }

    return convRes_;
}
