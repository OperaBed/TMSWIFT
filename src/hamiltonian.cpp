void hamiltonian(PetscErrorCode ierr, params *params, Mat &H)
{
	PetscInt start,end;
	int_params iparams;
	int_params jparams;
	h_params hparams;
	PetscReal m1s,m2s,x1,x2,k1,k2,mu1,mu2,th1,th2;
	PetscReal h_elem;
        PetscReal *_mu;
        PetscReal *_th;
        Vec MU_SEQ;
        Vec TH_SEQ;
        PetscReal *_wmu;
        PetscReal *_wth;
        Vec WMU_SEQ;
        Vec WTH_SEQ;
        VecScatter ctm;
        VecScatter ctwm;
        VecScatter ctt;
        VecScatter ctwt;
	PetscViewer view_H;

	ierr = PetscPrintf(PETSC_COMM_WORLD,"---Beginning Hamiltonian Construction---\n");CHKERRV(ierr);
	print_progress_header(ierr);	

	hparams.flag_anni=params->flag_anni;
	hparams.flag_mix=params->flag_mix;
	hparams.flag_asy=params->flag_asy;
	hparams.alpha=params->alpha;

	ierr = MatCreate(PETSC_COMM_WORLD,&H);CHKERRV(ierr);
  	ierr = MatSetSizes(H,PETSC_DECIDE,PETSC_DECIDE,4*params->N_tot,4*params->N_tot);CHKERRV(ierr);
  	ierr = MatSetFromOptions(H);CHKERRV(ierr);
  	ierr = MatSetUp(H);CHKERRV(ierr);

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

	ierr = MatGetOwnershipRange(H,&start,&end);CHKERRV(ierr);

  	for (PetscInt i=start;i<end;i++)
	{
	
		get_index(i,&iparams,params);
		m1s=params->m[iparams.index_f]*params->m[iparams.index_f];
		mu1=_mu[params->nm_tot[iparams.index_f]+iparams.index_m];
		th1=_th[params->nt_tot[iparams.index_f]+iparams.index_t];
		x1=(1.0+mu1*th1/sqrt(m1s+mu1*mu1))/2.0;
		k1=mu1*sqrt(1.0-th1*th1);

/*		hparams.m1=params->m[iparams.index_f];
		hparams.mu1=mu1;
		hparams.th1=th1;
		hparams.x1=x1;
		hparams.k1=k1;
		hparams.wmu1=_wmu[params->nm_tot[iparams.index_f]+iparams.index_m];
		hparams.wth1=_wth[params->nt_tot[iparams.index_f]+iparams.index_t];

		hparams.index_all1=i;
		hparams.index_f1=iparams.index_f;
		hparams.index_m1=iparams.index_m;
		hparams.index_t1=iparams.index_t;
		hparams.index_s1=iparams.index_s;
		ierr = PetscPrintf(PETSC_COMM_SELF,"h_elem: %f %f %d %d %d\n",hparams.x1,x1,hparams.index_t1,hparams.index_s1,hparams.index_f1);CHKERRV(ierr);
*/		for (PetscInt j=0;j<4*params->N_tot;j++)
		{
			hparams.flag_perm=0;
			hparams.m1=params->m[iparams.index_f];
                	hparams.mu1=mu1;
                	hparams.th1=th1;
                	hparams.x1=x1;
                	hparams.k1=k1;
                	hparams.wmu1=_wmu[params->nm_tot[iparams.index_f]+iparams.index_m];
                	hparams.wth1=_wth[params->nt_tot[iparams.index_f]+iparams.index_t];
			hparams.Jz=params->Jz;
			
                	hparams.index_all1=i;
                	hparams.index_f1=iparams.index_f;
                	hparams.index_m1=iparams.index_m;
                	hparams.index_t1=iparams.index_t;
                	hparams.index_s1=iparams.index_s;
	
	                get_index(j,&jparams,params);
			m2s=params->m[jparams.index_f]*params->m[jparams.index_f];
			mu2=_mu[params->nm_tot[jparams.index_f]+jparams.index_m];
			th2=_th[params->nt_tot[jparams.index_f]+jparams.index_t];
			x2=(1.0+mu2*th2/sqrt(m2s+mu2*mu2))/2.0;
			k2=mu2*sqrt(1.0-th2*th2);
			
			hparams.m2=params->m[jparams.index_f];	
	                hparams.mu2=mu2;
        	        hparams.th2=th2;
    	           	hparams.x2=x2;
       		        hparams.k2=k2;
			hparams.wmu2=_wmu[params->nm_tot[jparams.index_f]+jparams.index_m];
			hparams.wth2=_wth[params->nt_tot[jparams.index_f]+jparams.index_t];

	                hparams.index_all2=j;
        	        hparams.index_f2=jparams.index_f;
   		        hparams.index_m2=jparams.index_m;
    	   	        hparams.index_t2=jparams.index_t;
                	hparams.index_s2=jparams.index_s;

			h_elem=params->alpha*interaction(ierr,params,&hparams);	
//			ierr = PetscPrintf(PETSC_COMM_SELF,"%f\n",h_elem);CHKERRV(ierr);
			ierr = MatSetValue(H,i,j,h_elem,ADD_VALUES);CHKERRV(ierr);
    			if(i==j)
			{	
				h_elem=4.0*(m1s+mu1*mu1);
//				ierr = PetscPrintf(PETSC_COMM_SELF,"%f\n",h_elem);CHKERRV(ierr);
				ierr = MatSetValue(H,i,j,h_elem,ADD_VALUES);CHKERRV(ierr);
			}
		}
                if((i-start)%((end-start)/4)==0){ierr = PetscPrintf(PETSC_COMM_SELF,"|");CHKERRV(ierr);}
  	}

	ierr = MatAssemblyBegin(H,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
  	ierr = MatAssemblyEnd(H,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);

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

//	ierr = PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_DENSE);//DENSE<-->COMMON
//	MatView(H,PETSC_VIEWER_STDOUT_WORLD);

	PetscViewerBinaryOpen(PETSC_COMM_WORLD,params->hfile.c_str(),FILE_MODE_WRITE,&view_H);
	MatView(H,view_H);
	PetscViewerDestroy(&view_H);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"\n---Finishing Hamiltonian Construction---\n");CHKERRV(ierr);
}
