/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "turbulenceDamping.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(turbulenceDamping, 0);

    addToRunTimeSelectionTable
    (
        option,
        turbulenceDamping,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

volScalarField::Internal Foam::fv::turbulenceDamping::calculateSource
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    const volScalarField& Alpha =
        mesh().lookupObject<volScalarField>(primaryPhaseName_);

    const volVectorField grad_Alpha = fvc::grad(Alpha);
    const volScalarField grad_Alpha_mag = mag(grad_Alpha);

    // calculate interfacial area density
    volScalarField::Internal A1 = 2.0*Alpha*grad_Alpha_mag;
    volScalarField::Internal A2 = 2.0*(1.0-Alpha)*grad_Alpha_mag;


    // calculate the inverse of the length scale
    const volScalarField::Internal& V = mesh_.V();

    volScalarField oneByDn
    (
        IOobject
        (
            "oneByDn",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("oneByDn", dimless/dimLength, 0.0)
    );

    if (lengthScale_ == "FA")
    {
        const labelUList& owner = mesh_.owner();
        const labelUList& neighbour = mesh_.neighbour();

        const surfaceVectorField& sf = mesh_.Sf();

        forAll(owner, facei)
        {
            oneByDn[owner[facei]] +=
                mag(sf[facei] & grad_Alpha[owner[facei]]);

            oneByDn[neighbour[facei]] +=
                mag(sf[facei] & grad_Alpha[neighbour[facei]]);
        }

        forAll(mesh_.boundary(), patchi)
        {
            const labelUList& pFaceCells =
                mesh_.boundary()[patchi].faceCells();

            const fvsPatchField<vector>& psf = sf.boundaryField()[patchi];

            forAll(mesh_.boundary()[patchi], facei)
            {
                oneByDn[pFaceCells[facei]] +=
                    mag(psf[facei] & grad_Alpha[pFaceCells[facei]]);
            }
        }

       forAll(oneByDn, celli)
       {
           if (grad_Alpha_mag[celli] > SMALL)
           {
               oneByDn[celli] *= 0.5/V[celli]/grad_Alpha_mag[celli];
           }
           else
           {
               oneByDn[celli] = 0;
           }
       }

    }

    else if (lengthScale_ == "cubeRoot")
    {
        oneByDn.ref() = pow(V,-1.0/3.0);
    }


    // calculate separate damping terms
    volScalarField::Internal coeffs = 36.0*sqr(B_)/beta_*pow(oneByDn, 3.0);
    volScalarField::Internal source1 = coeffs*A1*rho1_*sqr(nu1_);
    volScalarField::Internal source2 = coeffs*A2*rho2_*sqr(nu2_);

    // calculate the total damping term
    dimensionedScalar heavy("heavy", dimless, 0.0);

    volScalarField::Internal source = 0.0 * source1;

    if (dampingTreatment_ == "heavyNegative")
    {
        if (rho1_ > rho2_)
        {
            heavy = - rho2_/rho1_*sqr(nu2_)/sqr(nu1_);

            source = source1*heavy + source2;
        }

        else
        {
            heavy = - rho1_/rho2_*sqr(nu1_)/sqr(nu2_);

            source = source1 + source2*heavy;
        }   
    }

    else if (dampingTreatment_ == "heavyZero")
    {
        if (rho1_ > rho2_)
        {
            source = source2;
        }

        else
        {
            source = source1;
        }   
    }

    else if (dampingTreatment_ == "symmetric")
    {
        source = sign(B_)*(source1 + source2);   
    }

    // return source term for omega equation
    if (fieldNames_[0] == "omega")
    {
        return source;
    }
    // return source term for epsilon equation
    else
    {
        const volScalarField& k = mesh().lookupObject<volScalarField>("k");

        return C2_*sqr(Cmu_)/beta_*k.internalField()*source;
    }
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::turbulenceDamping::turbulenceDamping
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(sourceName, modelType, dict, mesh),
    B_(readScalar(coeffs_.lookup("B"))),
    C2_(coeffs_.lookupOrDefault<scalar>("C2", 1.92)),
    beta_(coeffs_.lookupOrDefault<scalar>("beta", 0.075)),
    Cmu_(coeffs_.lookupOrDefault<scalar>("Cmu", 0.09)),
    lengthScale_(coeffs_.lookupOrDefault<word>("lengthScale", "FA")),
    dampingTreatment_
    (
        coeffs_.lookupOrDefault<word>("dampingTreatment", "heavyNegative")
    ),
    explicitSourceTreatment_
    (
        coeffs_.lookupOrDefault<Switch>("explicitSourceTreatment", true)
    ),
    transportProperties
    (
         mesh_.lookupObject<IOdictionary>
        (
            "transportProperties"
        )
    ),
    phase1Name_(wordList(transportProperties.lookup("phases"))[0]),
    phase2Name_(wordList(transportProperties.lookup("phases"))[1]),
    primaryPhaseName_("alpha." + phase1Name_),
    rho1_("rho", dimDensity, transportProperties.subDict(phase1Name_)),
    nu1_("nu", dimViscosity, transportProperties.subDict(phase1Name_)),
    rho2_("rho", dimDensity, transportProperties.subDict(phase2Name_)),
    nu2_("nu", dimViscosity, transportProperties.subDict(phase2Name_))
{
    coeffs_.lookup("fields") >> fieldNames_;

    if (fieldNames_.size() != 1)
    {
        FatalErrorInFunction
            << "settings are:" << fieldNames_ << exit(FatalError);
    }

    // only omega or epsilon is allowed
    if (fieldNames_[0] != "omega" && fieldNames_[0] != "epsilon")
    {
        FatalErrorInFunction
            << "The field is set to" << fieldNames_ 
            << ", it should be epsilon or omega!" << exit(FatalError);
    }

    // make sure the field name is consistent with the turbulence model
    // the following line should fail if inconsistent
    const volScalarField& epsilonOrOmega =
        mesh_.lookupObject<volScalarField>(fieldNames_[0]);

    Info << "Turbulence damping works in " << epsilonOrOmega.name() 
         << " mode"<< endl;

    applied_.setSize(fieldNames_.size(), false);

    Info << "B is set to " << B_.value() << endl;

    Info << "C2 is set to " << C2_.value() << endl;

    Info << "beta is set to " << beta_.value() << endl;

    Info << "Cmu is set to " << Cmu_.value() << endl;

    Info << "lengthScale is set to " << lengthScale_ << endl;

    Info << "dampingTreatment is set to " << dampingTreatment_ << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::turbulenceDamping::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    const volScalarField& Rho = mesh().lookupObject<volScalarField>("rho");

    const volScalarField& epsilonOrOmega =
        mesh().lookupObject<volScalarField>(fieldNames_[fieldi]);

    volScalarField::Internal source = calculateSource(eqn, fieldi)/epsilonOrOmega/Rho;

    eqn += fvm::Sp(source, epsilonOrOmega);
}


void Foam::fv::turbulenceDamping::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    const dimensionSet& dimEqn = eqn.dimensions();

    const volScalarField& epsilonOrOmega =
        mesh().lookupObject<volScalarField>(fieldNames_[fieldi]);

    const dimensionSet& dimEorW = epsilonOrOmega.dimensions();

    volScalarField::Internal source = calculateSource(eqn, fieldi);

    // divide density for strict incompressible turbulence models
    if (dimEqn == dimEorW/dimTime*dimVolume)
    {
        const volScalarField& Rho = mesh().lookupObject<volScalarField>("rho");

        source /= Rho;
    }

    if (explicitSourceTreatment_)
    {
        eqn += source;
    }
    else
    {
        eqn += fvm::Sp(source/epsilonOrOmega, epsilonOrOmega);
    }
}

bool Foam::fv::turbulenceDamping::read(const dictionary& dict)
{
    NotImplemented;

    return false;
}

// ************************************************************************* //
