PetscReal ct_function(PetscReal theta_local, void * int_params)
{
	struct int_params * iparams = (struct int_params *)int_params;
	h_params hparams;	
	PetscReal result;
	PetscReal mu1,mu2,th1,x1,x2,k1,k2,aa;
	PetscReal mf2;
	PetscReal denom;
	PetscReal p_b3;

	mu1 = iparams->mu1;
	th1 = theta_local;
	mu2 = iparams->mu2;
	mf2 = iparams->mf*iparams->mf;
	p_b3=iparams->p_bohr[0]*iparams->p_bohr[0]*iparams->p_bohr[0]*8.0;//Should this be /8 or times 8?

	result = 0.0;
	
	x1 = (1.0 + mu1*th1/sqrt(mf2 + mu1*mu1))/2.0;	
	k1 = mu1*sqrt(1.0-th1*th1);
	x2 = iparams->x2;	
	k2 = iparams->k2;
	
	hparams.flag_asy=iparams->flag_asy;
	hparams.flag_perm=0;
	hparams.x1=x1;
	hparams.k1=k1;
	hparams.x2=x2;
	hparams.k2=k2;
	hparams.Jz=iparams->Jzc;
	hparams.m1=iparams->mf;
	hparams.m2=iparams->mf;
	if (( std::fabs(x1-x2) >= 1e-8 || std::fabs(k1-k2) >= 1e-8 ) )
	{
		aa = (x1-x2)*(x1-x2)*mf2/2.0*(1.0/(1.0-x2)/(1.0-x1)+1.0/x1/x2) + k1*k1 + k2*k2+ (x1-x2)/2.0*(k1*k1*(1.0/(1.0-x1)-1.0/x1)-k2*k2*(1.0/(1.0-x2)-1.0/x2));

		denom = aa*aa-4.0*k1*k1*k2*k2;
	   	if( denom > 0.0 )
	    	{
                	hparams.A  = 1.0/sqrt(denom);
                	hparams.B  = (1.0-aa*hparams.A)/2.0;

			if (iparams->index_s==0 || iparams->index_s==3)   result = G1(&hparams);
			else		       result = G2(&hparams);
			result *= 2.0*mu1*mu1*x1*(1.0-x1)/sqrt(mf2 + mu1*mu1)/sqrt(mf2 + mu1*mu1)/sqrt(mf2 + mu1*mu1)*(1.0+mu2*mu2*mu2/p_b3)/(1.0+mu1*mu1*mu1/p_b3)/M_PI;
		}
	}
	return result;
}

PetscReal ct_integrand(PetscReal mu_local, void * int_params)
{
	struct int_params * iparams = (struct int_params *)int_params;

        size_t  neval;
        PetscReal a,b,abserr,result_coul;

        gsl_function F;

        iparams->mu1 = mu_local;
        a =-1.0;
        b = 1.0;
        F.function = &ct_function;
        F.params = iparams;

        gsl_integration_cquad_workspace *work = gsl_integration_cquad_workspace_alloc(iparams->N_work);
        gsl_integration_cquad(&F,a,b,iparams->epsabs,iparams->epsrel,work,&result_coul,&abserr,&neval);
        gsl_integration_cquad_workspace_free(work);
        return result_coul;
}	

void coulomb_trick(PetscErrorCode ierr, params *params)    
{       
	size_t neval;
        PetscReal a,b,abserr,coul_result;
	PetscReal mf2,mu2,th2;
	PetscInt start,end;
	PetscInt rank;
	PetscReal *_mu;
	PetscReal *_th;
  	Vec MU_SEQ;
  	Vec TH_SEQ;
	VecScatter     ctm;
	VecScatter     ctt;

	int_params int_params;

	int_params.epsabs=1e-4;
	int_params.epsrel=1e-8;
	int_params.N_work=4000;
	int_params.flag_asy=params->flag_asy;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);	
	for(PetscInt i=0;i<params->nf;i++)
	{
		int_params.p_bohr.push_back(params->p_bohr[i]);
	}

        gsl_function F;
        F.function = &ct_integrand;
        F.params=&int_params;

	ierr = PetscPrintf(PETSC_COMM_WORLD,"---Beginning Coulomb Trick---\n");CHKERRV(ierr);
       
	ierr = VecScatterCreateToAll(params->mu,&ctm,&MU_SEQ);CHKERRV(ierr);
  	ierr = VecScatterBegin(ctm,params->mu,MU_SEQ,INSERT_VALUES,SCATTER_FORWARD);CHKERRV(ierr);
  	ierr = VecScatterEnd(ctm,params->mu,MU_SEQ,INSERT_VALUES,SCATTER_FORWARD);CHKERRV(ierr);
  	ierr = VecGetArray(MU_SEQ,&_mu);CHKERRV(ierr);

        ierr = VecScatterCreateToAll(params->theta,&ctt,&TH_SEQ);CHKERRV(ierr);
        ierr = VecScatterBegin(ctt,params->theta,TH_SEQ,INSERT_VALUES,SCATTER_FORWARD);CHKERRV(ierr);
        ierr = VecScatterEnd(ctt,params->theta,TH_SEQ,INSERT_VALUES,SCATTER_FORWARD);CHKERRV(ierr);
        ierr = VecGetArray(TH_SEQ,&_th);CHKERRV(ierr);
  
	ierr = VecGetOwnershipRange(params->CT,&start,&end);CHKERRV(ierr);
        for (PetscInt i=start; i<end; i++)
	{
		if((i-start)%((end-start)/10)==0){ierr = PetscPrintf(PETSC_COMM_SELF,"----Cc: r%d -> %d\\%d\n",rank,(i-start),(end-start));CHKERRV(ierr);}

		get_index(i,&int_params,params);
//		ierr = PetscPrintf(PETSC_COMM_SELF,"CT_Index: %d %d %d %d %d\n",int_params.index_all,int_params.index_s,int_params.index_f,int_params.index_m,int_params.index_t);CHKERRV(ierr);

		int_params.mu2 = _mu[params->nm_tot[int_params.index_f]+int_params.index_m];
		int_params.th2 = _th[params->nt_tot[int_params.index_f]+int_params.index_t];
//	        ierr = PetscPrintf(PETSC_COMM_SELF,"Pts: %d %d %f %f\n",int_params.index_m,int_params.index_t,int_params.mu2,int_params.th2);CHKERRV(ierr);

               	if(int_params.index_s<2){int_params.Jzc = params->Jz;}
		else{int_params.Jzc= -params->Jz;}

		int_params.mf = params->m[int_params.index_f];

	
		mf2=int_params.mf*int_params.mf;
		mu2=int_params.mu2;
		th2=int_params.th2;		
			
	      	int_params.x2 = (1.0+mu2*th2/sqrt(mf2+mu2*mu2))/2.0;
        	int_params.k2 = mu2*sqrt(1.0-th2*th2);

		a = 0.0;
                b = params->lambda[int_params.index_f]/2.0;

                if (int_params.Jzc !=0 || int_params.index_s%2==0)
                {
	        	gsl_integration_cquad_workspace *work = gsl_integration_cquad_workspace_alloc(int_params.N_work);
	                gsl_integration_cquad(&F,a,b,int_params.epsabs,int_params.epsrel,work,&coul_result,&abserr,&neval);
                   	gsl_integration_cquad_workspace_free(work);

			ierr  = VecSetValues(params->CT,1,&i,&coul_result,ADD_VALUES);CHKERRV(ierr);
	                if (params->Jz==0)
			{ 
				PetscInt ishift;
				ishift=i+(3-2*int_params.index_s);
	      			ierr  = VecSetValues(params->CT,1,&ishift,&coul_result,ADD_VALUES);CHKERRV(ierr);
			}
		}
	}
	
	ierr = VecRestoreArray(MU_SEQ,&_mu);CHKERRV(ierr);
	ierr = VecScatterDestroy(&ctm);CHKERRV(ierr);
	ierr = VecDestroy(&MU_SEQ);CHKERRV(ierr);

        ierr = VecRestoreArray(TH_SEQ,&_th);CHKERRV(ierr);
	ierr = VecScatterDestroy(&ctt);CHKERRV(ierr);
        ierr = VecDestroy(&TH_SEQ);CHKERRV(ierr);

	ierr = VecAssemblyBegin(params->CT);CHKERRV(ierr);
	ierr = VecAssemblyEnd(params->CT);CHKERRV(ierr);
    	
	ierr = PetscPrintf(PETSC_COMM_WORLD,"---Finishing Coulomb Trick---\n");CHKERRV(ierr);
}
