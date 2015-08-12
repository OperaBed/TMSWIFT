void write_output(PetscErrorCode ierr, params *params)
{
        PetscInt start,end;
	PetscReal xtemp,ktemp,asytemp;
	PetscReal m,mu,th,wmu,wth;
	char ofile[100];
	Vec x;
	Vec k;
	Vec asy;
	
        PetscViewer     viewer;
	
	int_params iparams;

        PetscReal *_mu;
        PetscReal *_th;
        PetscReal *_wmu;
        PetscReal *_wth;
        Vec MU_SEQ;
        Vec TH_SEQ;
        Vec WMU_SEQ;
        Vec WTH_SEQ;
        VecScatter     ctm;
        VecScatter     ctwm;
        VecScatter     ctt;
        VecScatter     ctwt;

	ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRV(ierr);
        ierr = VecSetSizes(x,PETSC_DECIDE,4*params->N_tot);CHKERRV(ierr);
        ierr = VecSetFromOptions(x);CHKERRV(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD,&k);CHKERRV(ierr);
        ierr = VecSetSizes(k,PETSC_DECIDE,4*params->N_tot);CHKERRV(ierr);
        ierr = VecSetFromOptions(k);CHKERRV(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD,&asy);CHKERRV(ierr);
        ierr = VecSetSizes(asy,PETSC_DECIDE,4*params->N_tot);CHKERRV(ierr);
        ierr = VecSetFromOptions(asy);CHKERRV(ierr);

        ierr = VecScatterCreateToAll(params->mu,&ctm,&MU_SEQ);CHKERRV(ierr);
        ierr = VecScatterBegin(ctm,params->mu,MU_SEQ,INSERT_VALUES,SCATTER_FORWARD);CHKERRV(ierr);
        ierr = VecScatterEnd(ctm,params->mu,MU_SEQ,INSERT_VALUES,SCATTER_FORWARD);CHKERRV(ierr);
        ierr = VecGetArray(MU_SEQ,&_mu);CHKERRV(ierr);

        ierr = VecScatterCreateToAll(params->theta,&ctt,&TH_SEQ);CHKERRV(ierr);
        ierr = VecScatterBegin(ctt,params->theta,TH_SEQ,INSERT_VALUES,SCATTER_FORWARD);CHKERRV(ierr);
        ierr = VecScatterEnd(ctt,params->theta,TH_SEQ,INSERT_VALUES,SCATTER_FORWARD);CHKERRV(ierr);
        ierr = VecGetArray(TH_SEQ,&_th);CHKERRV(ierr);

        ierr = VecScatterCreateToAll(params->wmu,&ctwm,&WMU_SEQ);CHKERRV(ierr);
        ierr = VecScatterBegin(ctwm,params->wmu,WMU_SEQ,INSERT_VALUES,SCATTER_FORWARD);CHKERRV(ierr);
        ierr = VecScatterEnd(ctwm,params->wmu,WMU_SEQ,INSERT_VALUES,SCATTER_FORWARD);CHKERRV(ierr);
        ierr = VecGetArray(WMU_SEQ,&_wmu);CHKERRV(ierr);

        ierr = VecScatterCreateToAll(params->wtheta,&ctwt,&WTH_SEQ);CHKERRV(ierr);
        ierr = VecScatterBegin(ctwt,params->wtheta,WTH_SEQ,INSERT_VALUES,SCATTER_FORWARD);CHKERRV(ierr);
        ierr = VecScatterEnd(ctwt,params->wtheta,WTH_SEQ,INSERT_VALUES,SCATTER_FORWARD);CHKERRV(ierr);
        ierr = VecGetArray(WTH_SEQ,&_wth);CHKERRV(ierr);

        ierr = VecGetOwnershipRange(x,&start,&end);CHKERRV(ierr);

	for (PetscInt i=start; i<end; i++)
        {
		get_index(i,&iparams,params);

	 	m=params->m[iparams.index_f];
	
		mu=_mu[params->nm_tot[iparams.index_f]+iparams.index_m];
                th=_th[params->nt_tot[iparams.index_f]+iparams.index_t];

		wmu=_wmu[params->nm_tot[iparams.index_f]+iparams.index_m];
		wth=_wth[params->nt_tot[iparams.index_f]+iparams.index_t];

        	asytemp = 1.0/sqrt(mu*mu*(m*m+mu*mu*(1.0-th*th))/2.0/sqrt(m*m+mu*mu)/sqrt(m*m+mu*mu)/sqrt(m*m+mu*mu)*wth*wmu);
		xtemp=(1.0 + mu*th/sqrt(m*m + mu*mu))/2.0;
                ktemp=mu*sqrt(1.0-th*th);
	
		ierr  = VecSetValues(x,1,&i,&xtemp,ADD_VALUES);CHKERRV(ierr);	
		ierr  = VecSetValues(k,1,&i,&ktemp,ADD_VALUES);CHKERRV(ierr);	
		ierr  = VecSetValues(asy,1,&i,&asytemp,ADD_VALUES);CHKERRV(ierr);	
	}

        ierr = VecAssemblyBegin(x);CHKERRV(ierr);
        ierr = VecAssemblyEnd(x);CHKERRV(ierr);
 
        ierr = VecAssemblyBegin(k);CHKERRV(ierr);
        ierr = VecAssemblyEnd(k);CHKERRV(ierr);
  
        ierr = VecAssemblyBegin(asy);CHKERRV(ierr);
        ierr = VecAssemblyEnd(asy);CHKERRV(ierr);

        ierr = VecRestoreArray(MU_SEQ,&_mu);CHKERRV(ierr);
        ierr = VecScatterDestroy(&ctm);CHKERRV(ierr);
        ierr = VecDestroy(&MU_SEQ);CHKERRV(ierr);

        ierr = VecRestoreArray(TH_SEQ,&_th);CHKERRV(ierr);
        ierr = VecScatterDestroy(&ctt);CHKERRV(ierr);
        ierr = VecDestroy(&TH_SEQ);CHKERRV(ierr);

        ierr = VecRestoreArray(WMU_SEQ,&_wmu);CHKERRV(ierr);
        ierr = VecScatterDestroy(&ctwm);CHKERRV(ierr);
        ierr = VecDestroy(&WMU_SEQ);CHKERRV(ierr);

        ierr = VecRestoreArray(WTH_SEQ,&_wth);CHKERRV(ierr);
        ierr = VecScatterDestroy(&ctwt);CHKERRV(ierr);
        ierr = VecDestroy(&WTH_SEQ);CHKERRV(ierr);

	strcpy(ofile,params->ofile_n);
	strcat(ofile,"_x");
	
	ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,ofile,&viewer);CHKERRV(ierr);
        ierr = VecView(x,viewer);CHKERRV(ierr);
      	ierr = PetscViewerDestroy(&viewer);CHKERRV(ierr);

	strcpy(ofile,params->ofile_n);
        strcat(ofile,"_k");
        
        ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,ofile,&viewer);CHKERRV(ierr);
        ierr = VecView(k,viewer);CHKERRV(ierr);
        ierr = PetscViewerDestroy(&viewer);CHKERRV(ierr);

	strcpy(ofile,params->ofile_n);
        strcat(ofile,"_asy");
        
        ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,ofile,&viewer);CHKERRV(ierr);
        ierr = VecView(asy,viewer);CHKERRV(ierr);
        ierr = PetscViewerDestroy(&viewer);CHKERRV(ierr);

	ierr = VecDestroy(&x);CHKERRV(ierr);
        ierr = VecDestroy(&k);CHKERRV(ierr);
        ierr = VecDestroy(&asy);CHKERRV(ierr);

}


