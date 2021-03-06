void eigensolver(PetscErrorCode ierr, params *params, Mat &H, int argc, char **argv)
{
	

	EPS		eps;             /* eigenproblem solver context */
  	EPSType		type;
  	ST             st;
	KSP            ksp;
  	PC             pc; 
	PetscReal	tol,error;
	PetscReal	lower,upper;
  	//PetscInt       nev=dim,maxit,its;
  	PetscInt       	nev,maxit,its,nconv;
  	Vec            	xr,xi;
  	PetscScalar   	kr,ki;
	PetscReal 	re,im;
  	PetscViewer    	viewer;
	PetscInt rank;
	PetscInt size;
	PetscInt slice;
	std::string eig_file_n;
	std::ofstream eig_file;	
	char ofile[100];

        MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
        MPI_Comm_size(PETSC_COMM_WORLD,&size);

  	ierr = PetscPrintf(PETSC_COMM_WORLD,"---Beginning Eigenvalue Solver---\n");CHKERRV(ierr);
 	ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRV(ierr);

        eig_file_n.append(params->ofile_n.c_str());
        eig_file_n.append("_eval");
        eig_file.open(eig_file_n.c_str(),std::ofstream::trunc);

	//Set operators. In this case, it is a standard eigenvalue problem
	ierr = EPSSetOperators(eps,H,NULL);CHKERRV(ierr);
	ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRV(ierr);  

	ierr = EPSSetType(eps,EPSKRYLOVSCHUR);CHKERRV(ierr);

  	ierr = EPSGetST(eps,&st);CHKERRV(ierr);
  	ierr = STSetType(st,STSINVERT);CHKERRV(ierr);
  
  	ierr = STGetKSP(st,&ksp);CHKERRV(ierr);
  	ierr = KSPSetType(ksp,KSPPREONLY);CHKERRV(ierr);
  	ierr = KSPGetPC(ksp,&pc);CHKERRV(ierr);
  	ierr = PCSetType(pc,PCCHOLESKY);CHKERRV(ierr);
 	slice = size/64;	
	if(params->flag_mumps!=0)
	{
		ierr = EPSKrylovSchurSetPartitions(eps,size);CHKERRV(ierr);
	}else{
        	ierr = EPSKrylovSchurSetPartitions(eps,slice);CHKERRV(ierr);
        	ierr = EPSKrylovSchurSetDetectZeros(eps,PETSC_TRUE);CHKERRV(ierr);  /* enforce zero detection */
        	ierr = PCFactorSetMatSolverPackage(pc,MATSOLVERMUMPS);CHKERRV(ierr);
        	ierr = PetscOptionsInsertString("-mat_mumps_icntl_13 1 -mat_mumps_icntl_24 1 -mat_mumps_cntl_3 1e-12");CHKERRV(ierr);
	}

        strcpy(ofile,params->ofile_n.c_str());
        strcat(ofile,"_evecr");
        ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,ofile,&viewer);CHKERRV(ierr);

	for(PetscInt i=0;i<params->nf;i++){
	lower=std::pow(2.0*params->m[i]-params->m[i]*params->alpha*params->alpha,2.0);
//	upper=4.0*params->m[i]*params->m[i];
	upper=std::pow(2.0*params->m[i]-params->m[i]*params->alpha*params->alpha/2.0/36.0,2.0);
	
	ierr = EPSSetInterval(eps,lower,upper);
	ierr = EPSSetWhichEigenpairs(eps,EPS_ALL);
    	//Set solver parameters at runtime
  	ierr = EPSSetFromOptions(eps);CHKERRV(ierr);
// 	ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL);

   	ierr = EPSSolve(eps);CHKERRV(ierr);

	ierr = MatCreateVecs(H,NULL,&xr);CHKERRV(ierr);
	ierr = MatCreateVecs(H,NULL,&xi);CHKERRV(ierr);



   	ierr = EPSGetIterationNumber(eps,&its);CHKERRV(ierr);
  	ierr = PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %D\n",its);CHKERRV(ierr);

 
  	ierr = EPSGetType(eps,&type);CHKERRV(ierr);
  	ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRV(ierr);
  	ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRV(ierr);
  	ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev);CHKERRV(ierr);
  	ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRV(ierr);
  	ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%D\n",tol,maxit);CHKERRV(ierr);

	ierr = EPSGetConverged(eps,&nconv);CHKERRV(ierr);
  	ierr = PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %D\n",nconv);CHKERRV(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD," Range of eigenpairs: %f to %f \n\n",lower,upper);CHKERRV(ierr);

	eig_file.precision(12);

	if (nconv>0) 
	{
    		ierr = PetscPrintf(PETSC_COMM_WORLD,
         		"           k          ||Ax-kx||/||kx||\n"
         		"   ----------------- ------------------\n");CHKERRV(ierr);

		for (PetscInt i=0;i<nconv;i++)
		{
			//Get converged eigenpairs: i-th eigenvalue is stored in kr (real part) and ki (imaginary part)
      			ierr = EPSGetEigenpair(eps,i,&kr,&ki,xr,xi);CHKERRV(ierr);
         		//Compute the relative error associated to each eigenpair
     			ierr = EPSComputeError(eps,i,EPS_ERROR_RELATIVE,&error);CHKERRV(ierr);
                
			#if defined(PETSC_USE_COMPLEX)
      				re = PetscRealPart(kr);
      				im = PetscImaginaryPart(kr);
			#else
		      		re = kr;
  		    		im = ki;
			#endif

      			if (im!=0.0)
			{
        		
				ierr = PetscPrintf(PETSC_COMM_WORLD," %12f%+12f j %12g\n",re,im,error);CHKERRV(ierr);
				if(rank==0) eig_file << re << " " << im << " " << error << std::endl;
			} else 
			{
        			ierr = PetscPrintf(PETSC_COMM_WORLD,"   %12f       %12g\n",re,error);CHKERRV(ierr);
				if(rank==0) eig_file << re << " " << 0 << " " << error << std::endl;
     			}

                        ierr = VecView(xr,viewer);CHKERRV(ierr);

		}
   	 	ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRV(ierr);
  	}
	}
	eig_file.close();
	ierr = EPSDestroy(&eps);CHKERRV(ierr);
	ierr = PetscViewerDestroy(&viewer);CHKERRV(ierr);
	ierr = VecDestroy(&xr);CHKERRV(ierr);
	ierr = VecDestroy(&xi);CHKERRV(ierr);
	
  	ierr = PetscPrintf(PETSC_COMM_WORLD,"---Finishing Eigenvalue Solver---\n");CHKERRV(ierr);
}

