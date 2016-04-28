void print_input(PetscErrorCode ierr, params *params)
{

        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n---TMSWIFT Parameters---");CHKERRV(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n-Output File: %s",params->ofile_n.c_str());CHKERRV(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n-Number of Flavors: %d",params->nf);CHKERRV(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n-J_z: %d",params->Jz);CHKERRV(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n-Alpha: %g",params->alpha);CHKERRV(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n-N_tot: %d",params->N_tot);CHKERRV(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n-N_mu Discretization Flag: %s",params->sflag_nm_dis.c_str());CHKERRV(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n-N_theta Discretization Flag: %s",params->sflag_nt_dis.c_str());CHKERRV(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n-Annihilation Flag: %d",params->flag_anni);CHKERRV(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n-Flavor Mixing Flag: %d",params->flag_mix);CHKERRV(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n-Asymptotic G2 Flag: %d",params->flag_asy);CHKERRV(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n-MUMPS Flag: %d",params->flag_mumps);CHKERRV(ierr);

        for( int i=0;i<params->nf;i++)
        {
                ierr = PetscPrintf(PETSC_COMM_WORLD,"\n------");CHKERRV(ierr);
                ierr = PetscPrintf(PETSC_COMM_WORLD,"\n-Mass_%d: %g",i,params->m[i]);CHKERRV(ierr);
                ierr = PetscPrintf(PETSC_COMM_WORLD,"\n-P_Bohr_%d: %g",i,params->p_bohr[i]);CHKERRV(ierr);
                ierr = PetscPrintf(PETSC_COMM_WORLD,"\n-N_mu_%d: %d",i,params->nm[i]);CHKERRV(ierr);
                ierr = PetscPrintf(PETSC_COMM_WORLD,"\n-N_theta_%d: %d",i,params->nt[i]);CHKERRV(ierr);
                ierr = PetscPrintf(PETSC_COMM_WORLD,"\n-Lambda_%d: %g",i,params->lambda[i]);CHKERRV(ierr);
        }
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n------------\n");CHKERRV(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"%d %d %g %d %d %d %d %d",params->nf,params->Jz,params->alpha,params->flag_nm_dis,params->flag_nt_dis,params->flag_anni,params->flag_mix,params->flag_asy,params->flag_mumps);CHKERRV(ierr);
	for( int i=0;i<params->nf;i++){ierr = PetscPrintf(PETSC_COMM_WORLD," %g %d %d %g",params->m[i],params->nm[i],params->nt[i],params->lambda[i]);CHKERRV(ierr);}
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n------------\n");CHKERRV(ierr);


}


