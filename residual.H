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
    volVectorField convRes_
    (
        IOobject
        (
            "convRes",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("convRes",dimAcceleration,vector::zero)
        // residual_(fvc::div(phiRes_, URes_, "div(phi,U)"))
    );

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
            convRes_ = residual_(fvc::div(phiRes_, URes_, "div(phi,U)"));
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
