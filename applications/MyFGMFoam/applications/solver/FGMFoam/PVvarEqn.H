
if(!laminarBool && kEpsilonBool){
    XiPV = turbulence->epsilon() / (turbulence->k() + smallk) * varPV; 
}
if(!laminarBool && LESBool){
    scalarField& cellVolumeCells = cellVolumes.primitiveFieldRef();
    cellVolumeCells = mesh.V();
    XiPV = pow( CEpsilon * turbulence->mut() * pow( std::sqrt(2.0) * mag(symm(fvc::grad(U))) , 2) / ( rho * pow(cellVolumes,2.0/3.0) ), 1.0/3.0 ) * varPV;
}

// Solve the progress variable variance equation
tmp<fvScalarMatrix> tPVvarEqn
(
    (
    fvm::ddt(rho, varPV)
  + fvm::div(phi, varPV)
  - fvm::laplacian( (turbulence->muEff()/Sct), varPV)
 ==
    2*turbulence->mut()/Sct*magSqr(fvc::grad(PV))
  - 2*rho*XiPV
  + ( turbulence->YWI() - PV*reaction->SourcePV())
  + fvOptions(rho, varPV)
    )
);
fvScalarMatrix& PVvarEqn = tPVvarEqn.ref();

PVvarEqn.relax();
fvOptions.constrain(PVvarEqn);

PVvarEqn.solve("varPV");
fvOptions.correct(varPV);
varPV = max(   min( varPV, 1.0 ), 0.0   );
Info<< "PVvar min/max = " << min(varPV).value() << "/" << max(varPV).value() << endl;

reaction->correct();