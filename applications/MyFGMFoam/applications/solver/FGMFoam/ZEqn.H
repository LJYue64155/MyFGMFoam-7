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

tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Zmean)")
    )
);

Switch IsPremixed(chemistryProperties.lookup("premixed"));
if (!IsPremixed)
{
   tmp<fvScalarMatrix> tZEqn
   (
        (
	  fvm::ddt(rho, Zmean)
	+ mvConvection->fvmDiv(phi, Zmean)
        - fvm::laplacian( (turbulence->muEff()/reaction->Sct()), Zmean) //changed by senbin
        )
    );

    fvScalarMatrix& ZEqn = tZEqn.ref();
    ZEqn.relax();

    ZEqn.solve("Zmean");

    turbulence->correctChiZ();
    turbulence->correctVarZ();
    // reaction->correct();
    
}    
    Zmean = max(   min( Zmean, 1.0 ), 0.0   );
    Info<< "Z min/max = " << min(Zmean).value() << "/" << max(Zmean).value() << endl;
    

    

