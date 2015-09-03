void ct_discrete(PetscErrorCode ierr,params *params)
{
        PetscReal aa,result,spinfunction;
	PetscReal m,m2,x1,k1,x2,k2,mu1,mu2;
	PetscInt n;
	PetscInt start,end;
        PetscReal *_mu;
        PetscReal *_th;
	PetscReal *_wmu;
	PetscReal *_wth;
	PetscReal p_b3;
	h_params hp;	
	int_params iparams;
	int_params jparams;
        Vec MU_SEQ;
        Vec TH_SEQ;
        Vec WMU_SEQ;
        Vec WTH_SEQ;
        VecScatter     ctm;
        VecScatter     ctwm;
        VecScatter     ctt;
        VecScatter     ctwt;

        hp.flag_asy=params->flag_asy;
        hp.Jz=params->Jz;

        ierr = PetscPrintf(PETSC_COMM_WORLD,"---Beginning Coulomb Discrete---\n");CHKERRV(ierr);
	print_progress_header(ierr);

	p_b3=params->p_bohr[0]*params->p_bohr[0]*params->p_bohr[0]*8.0;//Should this be /8 or times 8?

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

        ierr = VecGetOwnershipRange(params->CT,&start,&end);CHKERRV(ierr);
	for (PetscInt i=start; i<end; i++)
        {
	        result = 0.0;
	        spinfunction=0.0;

		get_index(i,&iparams,params);

		hp.m2=params->m[iparams.index_f];
		hp.m1=params->m[iparams.index_f];
		hp.mu2=_mu[params->nm_tot[iparams.index_f]+iparams.index_m];
		hp.th2=_th[params->nt_tot[iparams.index_f]+iparams.index_t];
		hp.x2=(1.0 + hp.mu2*hp.th2/sqrt(hp.m2*hp.m2 + hp.mu2*hp.mu2))/2.0;
		hp.k2=hp.mu2*sqrt(1.0-hp.th2*hp.th2);
		mu2=hp.mu2;	
	        for (PetscInt j = 0; j < params->nm[iparams.index_f]; ++j )
       		{
			hp.mu1=_mu[params->nm_tot[iparams.index_f]+j];
			hp.wmu1=_wmu[params->nm_tot[iparams.index_f]+j];
	
            		for (PetscInt k = 0; k < params->nt[iparams.index_f]; ++k )
            		{	
				hp.th1=_th[params->nt_tot[iparams.index_f]+k];
				hp.wth1=_wth[params->nt_tot[iparams.index_f]+k];
                		hp.x1 = (1.0 + hp.mu1*hp.th1/sqrt(hp.m1*hp.m1 + hp.mu1*hp.mu1))/2.0;	
                		hp.k1 = hp.mu1*sqrt(1.0-hp.th1*hp.th1);
				hp.flag_perm=0;
		                set_hel_params(m,x1,x2,k1,k2,n,&hp);
				mu1=hp.mu1;
				m2=m*m;
				if ( std::fabs(x1-x2) >= 1e-8 || std::fabs(k1-k2) >= 1e-8 ) 
                		{
                    			aa = (x1-x2)*(x1-x2)*m2/2.0*(1.0/(1.0-x2)/(1.0-x1)+1.0/x1/x2) + k1*k1 + k2*k2 + (x1-x2)/2.0*(k1*k1*(1.0/(1.0-x1)-1.0/x1) - k2*k2*(1.0/(1.0-x2)-1.0/x2));
                    			hp.A  = 1.0/sqrt( aa*aa - 4.0*k1*k1*k2*k2 );
                			hp.B  = (1.0-aa*hp.A)/2.0;
 					switch(iparams.index_s)
					{
						case 0:
							spinfunction = -G1(&hp);
						break;

						case 1:
							spinfunction = -G2(&hp);
						break;

						case 2:
		                              		hp.Jz=-hp.Jz;
	        	                       		spinfunction = -G2(&hp);
						break;

						case 3:
		              	        		hp.Jz=-hp.Jz;
		                                	spinfunction = -G1(&hp);		
						break;
					}
					result += spinfunction*1.0/M_PI*hp.wth1*hp.wmu1*2.0*mu1*mu1*x1*(1.0-x1)/sqrt(m2 + mu1*mu1)/sqrt(m2 + mu1*mu1)/sqrt(m2 + mu1*mu1)*(1.0+mu2*mu2*mu2/p_b3)/(1.0+mu1*mu1*mu1/p_b3);
		                }	
			}
        	}
		
	        ierr  = VecSetValues(params->CT,1,&i,&result,ADD_VALUES);CHKERRV(ierr);
                if((i-start)%((end-start)/4)==0){ierr = PetscPrintf(PETSC_COMM_SELF,"|");CHKERRV(ierr);}
	}

        ierr = VecAssemblyBegin(params->CT);CHKERRV(ierr);
        ierr = VecAssemblyEnd(params->CT);CHKERRV(ierr);

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

        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n---Finishing Coulomb Discrete---\n");CHKERRV(ierr);
}

