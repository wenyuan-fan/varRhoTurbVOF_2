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

IOdictionary turbulencePropertiesDict
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

Switch varRho_
(
    turbulencePropertiesDict.lookupOrDefault<Switch>("varRho", true)
);

Info << "varRho is set to " << varRho_ << endl;

// Determine which fields should be used in the turbulence model
const volScalarField& rhoTurb = varRho_ ? rho : rhoOne;

const surfaceScalarField& rhoPhiTurb = varRho_ ? rhoPhi : phi;

