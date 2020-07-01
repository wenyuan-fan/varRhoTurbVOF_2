// Copyright held by original authors
#include "addToRunTimeSelectionTable.H"
#include "makeTurbulenceModel.H"
#include "IncompressibleTurbulenceModel.H"
#include "incompressible/transportModel/transportModel.H"
#include "fvOptions.H"
#include "RASModel.H"
#include "fvc.H"
#include "fvm.H"

namespace Foam
{
    typedef IncompressibleTurbulenceModel<transportModel>
        transportModelIncompressibleTurbulenceModel;
    typedef RASModel<transportModelIncompressibleTurbulenceModel>
        RAStransportModelIncompressibleTurbulenceModel;
}


#include "VOFKEpsilonIm.H"

makeTemplatedTurbulenceModel
(
    transportModelIncompressibleTurbulenceModel,
    RAS,
    VOFKEpsilonIm
);

