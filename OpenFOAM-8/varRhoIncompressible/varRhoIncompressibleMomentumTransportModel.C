/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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

#include "varRhoIncompressibleMomentumTransportModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(varRhoIncompressibleMomentumTransportModel, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::varRhoIncompressibleMomentumTransportModel::
varRhoIncompressibleMomentumTransportModel
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi
)
:
    momentumTransportModel
    (
        U,
        alphaRhoPhi,
        phi
    ),
    rho_(rho)
{}


Foam::tmp<Foam::volScalarField>
Foam::varRhoIncompressibleMomentumTransportModel::mu() const
{
    return rho_*nu();
}


Foam::tmp<Foam::scalarField>
Foam::varRhoIncompressibleMomentumTransportModel::mu(const label patchi) const
{
    return rho_.boundaryField()[patchi]*nu(patchi);
}


Foam::tmp<Foam::volScalarField>
Foam::varRhoIncompressibleMomentumTransportModel::mut() const
{
    return rho_*nut();
}


Foam::tmp<Foam::scalarField>
Foam::varRhoIncompressibleMomentumTransportModel::mut(const label patchi) const
{
    return rho_.boundaryField()[patchi]*nut(patchi);
}


Foam::tmp<Foam::volScalarField>
Foam::varRhoIncompressibleMomentumTransportModel::muEff() const
{
    return rho_*nuEff();
}


Foam::tmp<Foam::scalarField>
Foam::varRhoIncompressibleMomentumTransportModel::muEff(const label patchi) const
{
    return rho_.boundaryField()[patchi]*nuEff(patchi);
}


// ************************************************************************* //
