/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "fvCFD.H"
#include "upwind.H"
#include "OFstream.H"
#include "syncTools.H"
#include "vectorList.H"
#include "hexMatcher.H"
#include "prismMatcher.H"
#include "simpleControl.H"
#include "scalarMatrices.H"
#include "slicedSurfaceFields.H"
#include "volPointInterpolation.H"

#include "fvm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);
    #include "createFields.H"
    #include "createTimeControls.H"
    #include "initContinuityErrs.H"
    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Foam::Info << "\nStarting time loop\n" << endl;
    while (runTime.loop())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        Foam::Info<< "Time = " << runTime.timeName() << nl << endl;

        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
            + fvm::div(phi, U)
            ==
            fvm::laplacian(nu, U)
        );
        UEqn.relax();
        UEqn.solve();
        U.correctBoundaryConditions();

        phi = linearInterpolate(U) & mesh.Sf();
        adjustPhi(phi, U, p);
        surfaceScalarField phiU(phi);

        dimensionedScalar rAUf("rAUf", runTime.deltaT());

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvm::laplacian(rAUf,p) == fvc::div(phi)
            );
            
            pEqn.setReference(pRefCell, pRefValue);

            if (simple.finalNonOrthogonalIter())
            {
                pEqn.solve(mesh.solver("pFinal"));
            }
            else
            {
                pEqn.solve();
            }

            // if (simple.finalNonOrthogonalIter())
            // {
            //     phi -= pEqn.flux();
            // }
        }

        #include "continuityErrs.H"

        //U += fvc::reconstruct(phi - phiU);
        U += -rAUf * fvc::grad(p);
        U.correctBoundaryConditions();
        U.replace(vector::Z,0.0);

        runTime.write();

        Foam::Info << " Min(U) = " << min(U).value()
             << " Max(U) = " << max(U).value()
             << endl;

        // div = fvc::div(phi);
        // div.correctBoundaryConditions();
        // Info<< " Min(div) = " << min(div).value()
        //     << " Max(div) = " << max(div).value()
        //     << endl;

        Foam::Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Foam::Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
