    Info<< "Reading transportProperties\n" << endl;
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar nu
    (
        transportProperties.lookup("nu")
    );

    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createPhi.H" 
    volScalarField div
    (
        IOobject
        (
            "div",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::div(phi),
        zeroGradientFvPatchField<scalar>::typeName
    );
    
    // U.correctBoundaryConditions();
    // p.correctBoundaryConditions();

    // p.write();
    // U.write();
    
    // volVectorField fl_vof_0 = fvc::div(phi, U);
    // volVectorField fl_lap_0 = fvc::surfaceIntegrate(nu * (mesh.Sf() & linearInterpolate(fvc::grad(U))));

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, simple.dict(), pRefCell, pRefValue);
    mesh.setFluxRequired(p.name());
