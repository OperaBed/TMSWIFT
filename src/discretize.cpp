void discretize(PetscErrorCode ierr, params *params)
{
	ierr = PetscPrintf(PETSC_COMM_WORLD,"---Beginning Discretization---\n");CHKERRV(ierr);

	std::vector<PetscReal> weight;
	std::vector<PetscReal> abscis;

	PetscReal low_limit;
	PetscReal v_w;
	PetscReal v_n;
	
	PetscInt index_nf;
	index_nf=0;

	PetscInt index;

	for(PetscInt i=0; i<params->nf;i++)
	{
		low_limit=1.0/(1.0+params->lambda[i]/2.0/params->p_bohr[0]);		

		weight.resize(params->nm[i]);
		abscis.resize(params->nm[i]);

		switch(params->flag_nm_dis)
		{

			case 1:
			webbur::legendre_compute(params->nm[i],&abscis[0],&weight[0]);
			for (PetscInt j=0; j<params->nm[i]; j++)
                	{	
				index=j+index_nf;
				abscis[j]=(1.0-low_limit)/2.0*abscis[j]+(1.0+low_limit)/2.0;
		         	v_w = weight[j]*(1.0-low_limit)/2.0/abscis[j]/abscis[j]*params->p_bohr[0];
                		v_n = (1/abscis[j]-1)*params->p_bohr[0];
				ierr  = VecSetValues(params->wmu,1,&index,&v_w,INSERT_VALUES);CHKERRV(ierr);
				ierr  = VecSetValues(params->mu,1,&index,&v_n,INSERT_VALUES);CHKERRV(ierr);
                	}
			break;

			case 2:
	       		webbur::clenshaw_curtis_compute(params->nm[i],&abscis[0],&weight[0]);
			for (PetscInt j=0; j<params->nm[i]; j++)
                        {
                                index=j+index_nf;
                                v_w = weight[j];
                                abscis[j]=(1.0-low_limit)/2.0*abscis[j]+(1.0+low_limit)/2.0;
                                v_n = (1/abscis[j]-1)*params->p_bohr[0];
                                ierr  = VecSetValues(params->wmu,1,&index,&v_w,INSERT_VALUES);CHKERRV(ierr);
                                ierr  = VecSetValues(params->mu,1,&index,&v_n,INSERT_VALUES);CHKERRV(ierr);
                        }
			break;
			
			case 3:
			webbur::chebyshev1_compute(params->nm[i],&abscis[0],&weight[0]);
			for (PetscInt j=0; j<params->nm[i]; j++)
                        {
                                index=j+index_nf;
                                v_w = weight[j]*std::sqrt(1.0-abscis[j]*abscis[j]);
                                abscis[j]=(1.0-low_limit)/2.0*abscis[j]+(1.0+low_limit)/2.0;
                                v_n = (1/abscis[j]-1)*params->p_bohr[0];
                                ierr  = VecSetValues(params->wmu,1,&index,&v_w,INSERT_VALUES);CHKERRV(ierr);
                                ierr  = VecSetValues(params->mu,1,&index,&v_n,INSERT_VALUES);CHKERRV(ierr);
                        }
			break;

                        case 4:
                        webbur::gen_laguerre_compute(params->nm[i],0,&abscis[0],&weight[0]);
                        for (PetscInt j=0; j<params->nm[i]; j++)
                        {
                                index=j+index_nf;
                                v_w = weight[j]*std::exp(abscis[j]);
                                abscis[j]=(1.0-low_limit)/2.0*abscis[j]+(1.0+low_limit)/2.0;
                                v_n = (1/abscis[j]-1)*params->p_bohr[0];
                                ierr  = VecSetValues(params->wmu,1,&index,&v_w,INSERT_VALUES);CHKERRV(ierr);
                                ierr  = VecSetValues(params->mu,1,&index,&v_n,INSERT_VALUES);CHKERRV(ierr);
                        }

                        break;

			case 5:
			for (PetscInt j=0; j<params->nm[i]; j++)
                        {
                                index=j+index_nf;
                                v_w = 1.0;
                                v_n = 1.0*j*params->lambda[i]/(1.0*params->nm[i]);
                                ierr  = VecSetValues(params->wmu,1,&index,&v_w,INSERT_VALUES);CHKERRV(ierr);
                                ierr  = VecSetValues(params->mu,1,&index,&v_n,INSERT_VALUES);CHKERRV(ierr);
                        }

                        break;


			

			default:
			webbur::legendre_compute(params->nm[i],&abscis[0],&weight[0]);
			for (PetscInt j=0; j<params->nt[i]; j++)
                        {
                        	index=j+index_nf;
                                abscis[j]=(1.0-low_limit)/2.0*abscis[j]+(1.0+low_limit)/2.0;
                                v_w = weight[j]*(1.0-low_limit)/2.0/abscis[j]/abscis[j]*params->p_bohr[0];
                                v_n = (1/abscis[j]-1)*params->p_bohr[0];
                                ierr  = VecSetValues(params->wmu,1,&index,&v_w,INSERT_VALUES);CHKERRV(ierr);
                                ierr  = VecSetValues(params->mu,1,&index,&v_n,INSERT_VALUES);CHKERRV(ierr);
			}
			break;


		
		}
		

		
		index_nf+=params->nm[i];
	}

	index_nf=0;

        ierr = VecAssemblyBegin(params->wmu);CHKERRV(ierr);
        ierr = VecAssemblyEnd(params->wmu);CHKERRV(ierr);
        ierr = VecAssemblyBegin(params->mu);CHKERRV(ierr);
        ierr = VecAssemblyEnd(params->mu);CHKERRV(ierr);

	for(PetscInt i=0; i<params->nf;i++)
	{
		weight.resize(params->nt[i]);
		abscis.resize(params->nt[i]);

		switch(params->flag_nt_dis)
		{

			case 1:
			webbur::legendre_compute(params->nt[i],&abscis[0],&weight[0]);
			for (PetscInt j=0; j<params->nt[i]; j++)
                	{	
				index=j+index_nf;
		         	v_w = weight[j];
                		v_n = abscis[j];
				ierr  = VecSetValues(params->wtheta,1,&index,&v_w,INSERT_VALUES);CHKERRV(ierr);
				ierr  = VecSetValues(params->theta,1,&index,&v_n,INSERT_VALUES);CHKERRV(ierr);
                	}
			break;

			case 2:
	       		webbur::clenshaw_curtis_compute(params->nt[i],&abscis[0],&weight[0]);
			 for (PetscInt j=0; j<params->nt[i]; j++)
                        {
                                index=j+index_nf;
                                v_w = weight[j];
                                v_n = abscis[j];
                                ierr  = VecSetValues(params->wtheta,1,&index,&v_w,INSERT_VALUES);CHKERRV(ierr);
                                ierr  = VecSetValues(params->theta,1,&index,&v_n,INSERT_VALUES);CHKERRV(ierr);
                        }
			break;
			
			case 3:
			webbur::chebyshev1_compute(params->nt[i],&abscis[0],&weight[0]);
                        for (PetscInt j=0; j<params->nt[i]; j++)
                        {
                                index=j+index_nf;
                                v_w = weight[j]*std::sqrt(1.0-abscis[j]*abscis[j]);
                                v_n = abscis[j];
                                ierr  = VecSetValues(params->wtheta,1,&index,&v_w,INSERT_VALUES);CHKERRV(ierr);
                                ierr  = VecSetValues(params->theta,1,&index,&v_n,INSERT_VALUES);CHKERRV(ierr);
                        }

			break;			

			default:
			webbur::legendre_compute(params->nt[i],&abscis[0],&weight[0]);
		 	for (PetscInt j=0; j<params->nt[i]; j++)
                        {
                                index=j+index_nf;
                                v_w = weight[j];
                                v_n = abscis[j];
                                ierr  = VecSetValues(params->wtheta,1,&index,&v_w,INSERT_VALUES);CHKERRV(ierr);
                                ierr  = VecSetValues(params->theta,1,&index,&v_n,INSERT_VALUES);CHKERRV(ierr);
                        }
			break;

		
		}
		
		index_nf+=params->nt[i];
	}

	ierr = VecAssemblyBegin(params->wtheta);CHKERRV(ierr);
        ierr = VecAssemblyEnd(params->wtheta);CHKERRV(ierr);
	ierr = VecAssemblyBegin(params->theta);CHKERRV(ierr);
        ierr = VecAssemblyEnd(params->theta);CHKERRV(ierr);


	ierr = PetscPrintf(PETSC_COMM_WORLD,"---Finishing Discretization---\n");CHKERRV(ierr);

}
