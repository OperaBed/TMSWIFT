void read_input(PetscErrorCode ierr, char **argv, params *params)
//void read_input(char **argv,char *ofile_n, PetscInt *nf, PetscInt *Jz, PetscReal *alpha,std::vector<PetscReal> &m,std::vector<PetscInt> &nm,std::vector<PetscInt> &nt,std::vector<PetscReal> &lambda)
{

        std::ifstream ifile;
        std::string ifile_n;
        std::string null;

	PetscReal datar;
	PetscReal datat;
	PetscInt  datai;
	PetscInt  N_mu_tot;
	PetscInt  N_theta_tot;

        ifile_n = argv[1];
        ifile.open(ifile_n.c_str());
        if(!ifile.is_open())
        {
	        ierr = PetscPrintf(PETSC_COMM_WORLD,"---Error:: Initialization File Doesn't Exist---\n");CHKERRV(ierr);
        }

        ifile >> params->ofile_n >> null;
        ifile >> params->nf >> null;
	ifile >> params->Jz >> null;
        ifile >> params->alpha >> null;
        ifile >> params->flag_nm_dis >> null;
        ifile >> params->flag_nt_dis >> null;
        ifile >> params->flag_anni >> null;
        ifile >> params->flag_mix >> null;
        ifile >> params->flag_asy >> null;
        ifile >> params->flag_mumps >> null;
        for(int i=0; i<params->nf; i++)
        {
		ifile >> null;
		
		ifile >> datar >> null;
                params->m.push_back(datar);
                params->p_bohr.push_back(datar*params->alpha/2.0);
		ifile >> datai >> null;
                params->nm.push_back(datai);
		ifile >> datai >> null;
                params->nt.push_back(datai);
		ifile >> datar >> null;
                params->lambda.push_back(datar);
        }
	set_sflags(params);
	ifile.close();

	N_mu_tot=0;
	N_theta_tot=0;
	params->N_tot=0;
	
	datat=0.0;
	datar=0.0;
	for( PetscInt i=0;i<params->nf;i++)
	{
		params->nm_tot.push_back(datar);
		params->nt_tot.push_back(datat);
		N_mu_tot+=params->nm[i];
		N_theta_tot+=params->nt[i];
		params->N_tot+=params->nm[i]*params->nt[i];
		params->N_tot_f.push_back(params->N_tot);
		datar+=params->nm[i];
		datat+=params->nt[i];
	}
	get_outfile_name(params);

	params->hfile.append(params->ofile_n);
        params->hfile.append("_H");

	
	ierr = VecCreate(PETSC_COMM_WORLD,&params->mu);CHKERRV(ierr);
        ierr = VecSetSizes(params->mu,PETSC_DECIDE,N_mu_tot);CHKERRV(ierr);
        ierr = VecSetFromOptions(params->mu);CHKERRV(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD,&params->wmu);CHKERRV(ierr);
        ierr = VecSetSizes(params->wmu,PETSC_DECIDE,N_mu_tot);CHKERRV(ierr);
        ierr = VecSetFromOptions(params->wmu);CHKERRV(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD,&params->theta);CHKERRV(ierr);
        ierr = VecSetSizes(params->theta,PETSC_DECIDE,N_theta_tot);CHKERRV(ierr);
        ierr = VecSetFromOptions(params->theta);CHKERRV(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD,&params->wtheta);CHKERRV(ierr);
        ierr = VecSetSizes(params->wtheta,PETSC_DECIDE,N_theta_tot);CHKERRV(ierr);
        ierr = VecSetFromOptions(params->wtheta);CHKERRV(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD,&params->CT);CHKERRV(ierr);
        ierr = VecSetSizes(params->CT,PETSC_DECIDE,4*params->N_tot);CHKERRV(ierr);
        ierr = VecSetFromOptions(params->CT);CHKERRV(ierr);
}

