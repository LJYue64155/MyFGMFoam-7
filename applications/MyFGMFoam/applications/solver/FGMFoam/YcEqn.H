
tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Ycmean)")
    )
);

Switch IsColdState(combustionProperties.optionalSubDict("FGMModelCoeffs").lookup("IsColdState"));
if (!IsColdState){
    tmp<fvScalarMatrix> tYcEqn
    (
      (
        fvm::ddt(rho, Ycmean)
      + mvConvection->fvmDiv(phi, Ycmean)
            - fvm::laplacian( (turbulence->muEff()/reaction->Sct()), Ycmean)   //-- changed by senbin, 2019-06-24 muEff=muSgs()/Sct_ + mu()/Sc_, containing rho already
            ==  
            reaction->SourceYc()   // 假设隐含rho
            //--- added by senbin, 2019-04-26
      + fvOptions(rho, Ycmean)
            )
        );
        
    fvScalarMatrix& YcEqn = tYcEqn.ref();
    YcEqn.relax();
    fvOptions.constrain(YcEqn);
        
    YcEqn.solve("Ycmean");
    fvOptions.correct(Ycmean);
    Ycmean = max(   min( Ycmean, 1.0 ), 0.0   );
    Info<< "Yc min/max = " << min(Ycmean).value() << "/" << max(Ycmean).value() << endl;
        
    turbulence->correctChiYc();
    turbulence->correctVarYc();
    // reaction->correct();   //correct SourcePV, T, sensor

}

Ycmean = max(   min( Ycmean, 1.0 ), 0.0   );
Info<< "Yc min/max = " << min(Ycmean).value() << "/" << max(Ycmean).value() << endl;
