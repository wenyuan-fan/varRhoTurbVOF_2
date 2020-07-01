/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2015 OpenFOAM Foundation
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

#include "VOFKEpsilonEx.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
VOFKEpsilonEx<BasicTurbulenceModel>::VOFKEpsilonEx
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    kEpsilon<BasicTurbulenceModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName,
        type
    )

{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



template<class BasicTurbulenceModel>
tmp<fvScalarMatrix>
VOFKEpsilonEx<BasicTurbulenceModel>::kSource() const
{
    const volScalarField& Rho =
        this->mesh_.objectRegistry::template
        lookupObject<volScalarField>("rho");


     volScalarField product = fvc::grad(Rho)&fvc::grad(this->k_);


    return fvm::SuSp(this->DkEff()/Rho/this->k_*product, this->k_);



}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix>
VOFKEpsilonEx<BasicTurbulenceModel>::epsilonSource() const
{
    const volScalarField& Rho =
        this->mesh_.objectRegistry::template
        lookupObject<volScalarField>("rho");

    volScalarField product = fvc::grad(Rho)&fvc::grad(this->epsilon_);


    return fvm::Sp(this->DepsilonEff()/Rho/this->epsilon_*product, this->epsilon_);



}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
