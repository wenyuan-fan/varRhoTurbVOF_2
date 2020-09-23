/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2020 OpenFOAM Foundation
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

Typedef
    Foam::incompressible::momentumTransportModel

Typedef
    Foam::incompressible::RASModel

Typedef
    Foam::incompressible::LESModel

Description
    Typedefs for turbulence, RAS and LES models for incompressible flow
    based on the standard laminar transport package. Density is explicity
    referenced to consider the variable-density effect, as described in:
    \verbatim
        Fan, W. & Anglart, H. (2020).
        varRhoTurbVOF: A new set of volume of fluid solvers for turbulent 
        isothermal multiphase flows in OpenFOAM.
        Computer Physics Communications, 247, 106876
    \endverbatim 

SourceFiles
    varRhoKinematicMomentumTransportModel.C
    varRhoKinematicMomentumTransportModels.C

\*---------------------------------------------------------------------------*/

#ifndef varRhoKinematicMomentumTransportModel_H
#define varRhoKinematicMomentumTransportModel_H

#include "VarRhoIncompressibleMomentumTransportModel.H"
#include "laminarModel.H"
#include "RASModel.H"
#include "LESModel.H"
#include "transportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace incompressible
    {
        typedef VarRhoIncompressibleMomentumTransportModel<transportModel>
            momentumTransportModel;

        typedef laminarModel<momentumTransportModel> laminarModel;
        typedef RASModel<momentumTransportModel> RASModel;
        typedef LESModel<momentumTransportModel> LESModel;

        template<class BasicCompressibleMomentumTransportModel>
        autoPtr<BasicCompressibleMomentumTransportModel> New
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const typename BasicCompressibleMomentumTransportModel::
                transportModel& transport
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "varRhoKinematicMomentumTransportModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
