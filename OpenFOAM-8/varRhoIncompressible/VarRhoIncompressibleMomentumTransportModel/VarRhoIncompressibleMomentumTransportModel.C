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

#include "VarRhoIncompressibleMomentumTransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class TransportModel>
Foam::VarRhoIncompressibleMomentumTransportModel<TransportModel>::
VarRhoIncompressibleMomentumTransportModel
(
    const word& type,
    const geometricOneField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const TransportModel& transport
)
:
    MomentumTransportModel
    <
        geometricOneField,
        volScalarField,
        varRhoIncompressibleMomentumTransportModel,
        TransportModel
    >
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport
    )
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class TransportModel>
Foam::autoPtr<Foam::VarRhoIncompressibleMomentumTransportModel<TransportModel>>
Foam::VarRhoIncompressibleMomentumTransportModel<TransportModel>::New
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const TransportModel& transport
)
{
    return autoPtr<VarRhoIncompressibleMomentumTransportModel>
    (
        static_cast<VarRhoIncompressibleMomentumTransportModel*>(
        MomentumTransportModel
        <
            geometricOneField,
            volScalarField,
            varRhoIncompressibleMomentumTransportModel,
            TransportModel
        >::New
        (
            geometricOneField(),
            rho,
            U,
            alphaRhoPhi,
            phi,
            transport
        ).ptr())
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TransportModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::VarRhoIncompressibleMomentumTransportModel<TransportModel>::devSigma() const
{
    return devTau();
}


template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::VarRhoIncompressibleMomentumTransportModel<TransportModel>::divDevSigma
(
    volVectorField& U
) const
{
    return divDevTau(U);
}


template<class TransportModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::VarRhoIncompressibleMomentumTransportModel<TransportModel>::
devTau() const
{
    NotImplemented;

    return devSigma();
}


template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::VarRhoIncompressibleMomentumTransportModel<TransportModel>::
divDevTau
(
    volVectorField& U
) const
{
    NotImplemented;

    return divDevSigma(U);
}


template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::VarRhoIncompressibleMomentumTransportModel<TransportModel>::
divDevTau
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    NotImplemented;

    return divDevSigma(U);
}


// ************************************************************************* //
