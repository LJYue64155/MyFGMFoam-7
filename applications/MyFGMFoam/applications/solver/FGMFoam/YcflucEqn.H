
if(!laminarBool && RASBool){
    XiPV = turbulence->epsilon() / (turbulence->k() + smallk) * Ycfluc; 
}
if(!laminarBool && LESBool){
    scalarField& cellVolumeCells = cellVolumes.primitiveFieldRef();
    cellVolumeCells = mesh.V();
    XiPV = pow( CEpsilon * turbulence->mut() * pow( std::sqrt(2.0) * mag(symm(fvc::grad(U))) , 2) / ( rho * pow(cellVolumes,2.0/3.0) ), 1.0/3.0 ) * Ycfluc;
}


if (!IsColdState)
{
    // Solve the progress variable variance equation
    tmp<fvScalarMatrix> tYcflucEqn
    (
        (
        fvm::ddt(rho, Ycfluc)
    + fvm::div(phi, Ycfluc)
    - fvm::laplacian( (turbulence->muEff()/reaction->Sct()), Ycfluc)
    ==
        2*turbulence->mut()/reaction->Sct()*magSqr(fvc::grad(Ycmean))
    - 2*rho*XiPV
    + ( turbulence->YWI() - Ycmean*reaction->SourceYc())
    + fvOptions(rho, Ycfluc)
        )
    );
    fvScalarMatrix& YcflucEqn = tYcflucEqn.ref();

    YcflucEqn.relax();
    fvOptions.constrain(YcflucEqn);

    YcflucEqn.solve("Ycfluc");
    fvOptions.correct(Ycfluc);

}

Ycfluc = max(   min( Ycfluc, 1.0 ), 0.0   );
Info<< "Ycfluc min/max = " << min(Ycfluc).value() << "/" << max(Ycfluc).value() << endl;

Info<< "T gas min/max = " << min(thermo.T()).value() << "/" << max(thermo.T()).value() << endl;