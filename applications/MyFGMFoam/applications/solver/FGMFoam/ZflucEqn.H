
if(!laminarBool && RASBool){
    XiZ = turbulence->epsilon() / (turbulence->k() + smallk) * Zfluc; 
}
if(!laminarBool && LESBool){
    scalarField& cellVolumeCells = cellVolumes.primitiveFieldRef();
    cellVolumeCells = mesh.V();
    XiZ = pow( CEpsilon * turbulence->mut() * pow( std::sqrt(2.0) * mag(symm(fvc::grad(U))) , 2) / ( rho * pow(cellVolumes,2.0/3.0) ), 1.0/3.0 ) * Zfluc;
}

tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Zfluc)")
    )
);

Switch useProgressVariableVariance(combustionProperties.optionalSubDict("FGMModelCoeffs").lookupOrDefault("useProgressVariableVariance",true));
Switch IsColdState(combustionProperties.optionalSubDict("FGMModelCoeffs").lookup("IsColdState"));
if (!IsPremixed and !IsColdState)
{
    // Solve the progress variable variance equation
    tmp<fvScalarMatrix> tZflucEqn
    (
        (
        fvm::ddt(rho, Zfluc)
        + mvConvection->fvmDiv(phi, Zfluc)
        - fvm::laplacian( (turbulence->muEff()/reaction->Sct()), Zfluc)
        ==
        + 2*turbulence->mut()/reaction->Sct()*magSqr(fvc::grad(Zmean))
        - 2*rho*XiZ
        )
    );
    fvScalarMatrix& ZflucEqn = tZflucEqn.ref();

    ZflucEqn.relax();
    fvOptions.constrain(ZflucEqn);

    ZflucEqn.solve("Zfluc");
    fvOptions.correct(Zfluc);

    // reaction->correct();
}

Zfluc = max(   min( Zfluc, 1.0 ), 0.0   );
Info<< "Zfluc min/max = " << min(Zfluc).value() << "/" << max(Zfluc).value() << endl;

