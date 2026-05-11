IOdictionary chemistryProperties		// senbin, 14-11-2018
(
    IOobject
    (
        "chemistryProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

IOdictionary combustionProperties		// senbin, 14-11-2018
(
    IOobject
    (
        "combustionProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

Switch IsPremixed(chemistryProperties.lookup("premixed"));
if (!IsPremixed)
{
   tmp<fvScalarMatrix> tZEqn
   (
        (
	  fvm::ddt(rho, Zmean)
	+ fvm::div(phi, Zmean)
        - fvm::laplacian( (turbulence->muEff()/reaction->Sct()), Zmean) //changed by senbin
        )
    );

    fvScalarMatrix& ZEqn = tZEqn.ref();
    ZEqn.relax();

    ZEqn.solve("Zmean");

    turbulence->correctChiZ();
    turbulence->correctVarZ();
    
}    
    Zmean = max(   min( Zmean, 1.0 ), 0.0   );
    Info<< "Z min/max = " << min(Zmean).value() << "/" << max(Zmean).value() << endl;
    

    

