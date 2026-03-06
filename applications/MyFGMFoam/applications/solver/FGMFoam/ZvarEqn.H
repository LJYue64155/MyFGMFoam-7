
if(!laminarBool && kEpsilonBool){
    XiZ = turbulence->epsilon() / (turbulence->k() + smallk) * varZ; 
}
if(!laminarBool && LESBool){
    scalarField& cellVolumeCells = cellVolumes.primitiveFieldRef();
    cellVolumeCells = mesh.V();
    XiZ = pow( CEpsilon * turbulence->mut() * pow( std::sqrt(2.0) * mag(symm(fvc::grad(U))) , 2) / ( rho * pow(cellVolumes,2.0/3.0) ), 1.0/3.0 ) * varZ;
}

// Solve the progress variable variance equation
tmp<fvScalarMatrix> tZvarEqn
(
    (
    fvm::ddt(rho, varZ)
  + fvm::div(phi, varZ)
  - fvm::laplacian( (turbulence->muEff()/Sct), varZ)
 ==
    2*turbulence->mut()/Sct*magSqr(fvc::grad(Z))
  - 2*rho*XiZ
    )
);
fvScalarMatrix& ZvarEqn = tZvarEqn.ref();

ZvarEqn.relax();
fvOptions.constrain(ZvarEqn);

ZvarEqn.solve("varZ");
fvOptions.correct(varZ);
varZ = max(   min( varZ, 1.0 ), 0.0   );
Info<< "Zvar min/max = " << min(varZ).value() << "/" << max(varZ).value() << endl;

reaction->correct();
