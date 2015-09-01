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

	params->hfile.append(params->ofile_n);
        for(int i = 0; i<params->nf; i++)
        { 
                params->hfile.append(std::to_string(static_cast<long long>(params->nm[i])));
                params->hfile.append(std::to_string(static_cast<long long>(params->nt[i])));
        }
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

void print_input(PetscErrorCode ierr, params *params)
{

        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n---TMSWIFT Parameters---");CHKERRV(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n-Output File: %s",params->ofile_n);CHKERRV(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n-Number of Flavors: %d",params->nf);CHKERRV(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n-J_z: %d",params->Jz);CHKERRV(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n-Alpha: %g",params->alpha);CHKERRV(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n-N_tot: %d",params->N_tot);CHKERRV(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n-N_mu Discretization Flag: %s",params->sflag_nm_dis.c_str());CHKERRV(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n-N_theta Discretization Flag: %s",params->sflag_nt_dis.c_str());CHKERRV(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n-Annihilation Flag: %d",params->flag_anni);CHKERRV(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n-Flavor Mixing Flag: %d",params->flag_mix);CHKERRV(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n-Asymptotic G2 Flag: %d",params->flag_asy);CHKERRV(ierr);

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

}

void set_sflags(params *params)
{

	switch (params->flag_nm_dis)
	{
		case 1:
		params->sflag_nm_dis="legendre";
		break;

		case 2:
		params->sflag_nm_dis="clenshaw_curtis";
		break;

		case 4:
		params->sflag_nm_dis="laguerre";
		break;

		case 3:
		params->sflag_nm_dis="chebyshev1";
		break;

		default:
		params->sflag_nm_dis="legendre";
		break;

	}

        switch (params->flag_nt_dis)
        {
                case 1:
                params->sflag_nt_dis="legendre";
                break;

                case 2:
                params->sflag_nt_dis="clenshaw_curtis";
                break;

                case 3:
                params->sflag_nt_dis="chebyshev1";
                break;

                default:
                params->sflag_nt_dis="legendre";
                break;

        }


}


