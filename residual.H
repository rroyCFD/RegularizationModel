surfaceScalarField Foam::RegularizationModel::surfaceFieldResidual
(const surfaceScalarField& ssf_)
{
    volVectorField vvf_(ssf_.name()+"_Res", fvc::reconstruct(ssf_));
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
        
    if(residualOrder_ == "A6"|| residualOrder_ == "a6")
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

        // Make filtered flux  divergence free for KE conservation
        // - Recommended for inviscid and laminar flow
        // - For LES of hig Re (high total viscosity, it may be ignored!)
        // - if not made divFree; U*div(phi) correction provides tiny improvment
        // Velocity div error is not affected by the regularization filter
        if(extpFilterFieldDivFree_)
        {
            setDivergenceFree (phiRes_);
        }

        convRes_ = residual_(convOperator(phiRes_, URes_));

        tURes_.clear();
        phiRes_.clear();
    }
    else
    {
        FatalErrorIn("Regularization order") << residualOrder_
            << " is not recognized!\n" << " Avaialble orders is A6"
            << abort(FatalError);
    }

    return convRes_;
}
