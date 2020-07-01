/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2019 OpenFOAM Foundation
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

Class
    Foam::fv::turbulenceDamping

Description
    Calculates and applies the interficial turbulence damping term to the
    omega or epsilon equation. Therefore, this fvOptions can only handle 
    omega- or epsilon-based RANS models. In addition, this implementation 
    only applies to multiphase flows where the flow is treated as a (VOF)
    mixture in turbulence modeling. Also, each phase should have constant
    properties. For backward compatibiliy, the current implementaion 
    supports turbulence models with and without the variable-density effect
    being considered.

    Reference:
    \verbatim
        Original model for omega-based equations:
            Egorov, Y. (2004).
            Validation of CFD codes with PTS-relevant test cases.
            5th Euratom Framework Programme ECORA project, pp. 91–116.

        Extension to epsilon-based equations:
            Frederix, E.M.A., Mathur, A., Dovizio, D., Geurts, B.J.,
            Komen, E.M.J. (2018).
            Reynolds-averaged modeling of turbulence damping near a
            large-scale interface in two-phase flow.
            Nuclear Engineering and Design, 333, pp. 122–130.

        The current implementation is based on:
            Fan, W. & Anglart, H. (2019).
            Progress in Phenomenological Modeling of Turbulence Damping
            around a Two-Phase Interface.
            Fluids, 4(3), 136.
    \endverbatim


Usage
    Example usage:
    \verbatim
    fields             (omega);   // Name of the field: omega or epsilon

    B                  0;         // Damping coefficient, 0 means no damping

    lengthScale        FA;        // Optional parameter to specify the method
                                  // to calculate the length scale for damping.
                                  // The alternative is "cubeRoot".

    dampingTreatment   heavyNegative; // Optional parameter to specify the 
                                      // treatment for the heavier phase.
                                      // Alternatives are "heavyZero" and
                                      // "symmetric".

    explicitSourceTreatment true; // Optional parameter to specify whether the
                                  // source term should be treated explicitly.
    \endverbatim

Author
    Wenyuan Fan

SourceFiles
    turbulenceDamping.C

\*---------------------------------------------------------------------------*/

#ifndef turbulenceDamping_H
#define turbulenceDamping_H
#include "fvCFD.H"
#include "fvOption.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
               Class turbulenceDamping Declaration
\*---------------------------------------------------------------------------*/

class turbulenceDamping
:
    public option
{
    // Private Data

        dimensionedScalar B_;

        dimensionedScalar C2_;

        dimensionedScalar beta_;

        dimensionedScalar Cmu_;

        word lengthScale_;

        word dampingTreatment_;

        Switch explicitSourceTreatment_;

        const dictionary& transportProperties;

        word phase1Name_;

        word phase2Name_;

        const word primaryPhaseName_;

        dimensionedScalar rho1_;

        dimensionedScalar nu1_;

        dimensionedScalar rho2_;

        dimensionedScalar nu2_;

        //- Source term to omega or epsilon equation
        volScalarField::Internal calculateSource
        (
            fvMatrix<scalar>& eqn,
            const label fieldi
        );

public:

    //- Runtime type information
    TypeName("turbulenceDamping");


    // Constructors

        //- Construct from explicit source name and mesh
        turbulenceDamping
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        turbulenceDamping(const turbulenceDamping&) = delete;


    // Member Functions

        // Evaluate

            //- Add turbulence damping to strict incompressible equation
            virtual void addSup
            (
                fvMatrix<scalar>& eqn,
                const label fieldi
            );

            //- Add turbulence damping to full-form equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldi
            );


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const turbulenceDamping&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
