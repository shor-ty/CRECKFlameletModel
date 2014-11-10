
	//  Solving equation of variance of mixture fraction
	fvScalarMatrix csiEqn
	(
		(
			fvm::div(phi, csi)
	      	      - fvm::Sp(fvc::div(phi), csi)		
              	      - fvm::laplacian(turbulence->muEff()/sigmat, csi)
		)
    	);
	csiEqn.relax();
	csiEqn.solve(mesh.solver("csi"));
	   

	if (varianceEquation == true)
	{
		// Solving equation of variance of mixture fraction
		    fvScalarMatrix csiv2Eqn
		(
			(
				fvm::div(phi, csiv2)
			      - fvm::Sp(fvc::div(phi), csiv2)		
			      - fvm::laplacian(turbulence->mut()/sigmat, csiv2)
			      -	Cg*turbulence->mut()*magSqr(fvc::grad(csi))
			      + Cd*rho*turbulence->epsilon()/(turbulence->k()+k_small)*csiv2
			) 
		);
		csiv2Eqn.relax();
		csiv2Eqn.solve(mesh.solver("csiv2"));
	}
	else
	{
		Info << "Algebraic equation for csiv2... " << endl;
    		scalarField& csiv2Cells = csiv2.internalField();

		csiv2Cells = Cg*turbulence->mut()*magSqr(fvc::grad(csi))*turbulence->k() / 
	  	      			(Cd*rho*(turbulence->epsilon()+epsilon_small));
	}
	
	// Updating stoichiometric scalar dissipation rate
	chi_st=Cx*turbulence->epsilon()/(turbulence->k()+k_small)*csiv2;


