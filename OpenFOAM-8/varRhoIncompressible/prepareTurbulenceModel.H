// create rho field for constant-density turbulence models
volScalarField rhoOne
(
    IOobject
    (
        "rhoOne",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("one", dimless, 1)
);

Switch varRho_;

IOdictionary momentumTransport
(
    IOobject
    (
        "momentumTransport",
        U.time().constant(),
        U.db(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    )
);

if (momentumTransport.typeHeaderOk<IOdictionary>(true))
{
    varRho_ = 
        momentumTransport.lookupOrDefault<Switch>("varRho", true);
}

else
{
    IOdictionary turbulenceProperties
    (
        IOobject
        (
            "turbulenceProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    varRho_ =
        turbulenceProperties.lookupOrDefault<Switch>("varRho", true);
}

Info << "varRho is set to " << varRho_ << endl;

// Determine which fields should be used in the turbulence model
const volScalarField& rhoTurb = varRho_ ? rho : rhoOne;

const surfaceScalarField& rhoPhiTurb = varRho_ ? rhoPhi : phi;

