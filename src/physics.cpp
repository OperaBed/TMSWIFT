PetscReal Int(PetscInt n, h_params *hp)
{
        PetscReal vz,result;
	PetscReal A,B,k1,k2;
	A=hp->A;
	B=hp->B;
	k1=hp->k1;
	k2=hp->k2;
        vz=1.0;
        if (n%2==0) vz = -1.0;
        if( std::fabs(k1)<1e-15 || std::fabs(k2)<1e-15 )
        {
                if(n==0){return 1.0;
                }else{return 0.0;}
        }else{
		result=vz*std::pow(A,1.0-n)*std::pow(B/k1/k2,1.0*n);
                return result;
        }
}


void set_hel_params(PetscReal &m,PetscReal &x1,PetscReal &x2,PetscReal &k1,PetscReal &k2,PetscInt &n,h_params *params)
{

        m=params->m1;            
        n=params->Jz;
	if((params->flag_perm==0))
	{
        	x1=params->x1;
        	x2=params->x2;
        	k1=params->k1;
        	k2=params->k2;
	}else
	{
                x2=params->x1;
                x1=params->x2;
                k2=params->k1;
                k1=params->k2;
	}

}

PetscReal general_spin(h_params *hp)
{
	PetscReal result;

	result=0.0;

	if(hp->index_f1==hp->index_f2)
	{
		if(hp->index_s2==0)
		{
			if(hp->index_s1==0){result=G1(hp);}
                	if(hp->index_s1==1){result=G3_star(hp);}
                	if(hp->index_s1==2){result=G3(hp);}
                	if(hp->index_s1==3){result=0.0;}
		}
		if(hp->index_s2==1)
		{	
			hp->flag_perm=1;
                        if(hp->index_s1==0){result=G3_star(hp);}
                        if(hp->index_s1==1){result=G2(hp);}
                        if(hp->index_s1==2){result=G4(hp);}
			hp->Jz=-hp->Jz;
                        if(hp->index_s1==3){result=-G3(hp);}
		}
		if(hp->index_s2==2)
		{
			hp->flag_perm=1;
                        if(hp->index_s1==0){result=G3(hp);}
                        if(hp->index_s1==1){result=G4(hp);}
			hp->Jz=-hp->Jz;
                        if(hp->index_s1==2){result=G2(hp);}
                        if(hp->index_s1==3){result=-G3_star(hp);}
		}
		if(hp->index_s2==3)
		{
			hp->Jz=-hp->Jz;
                        if(hp->index_s1==0){result=0.0;}
                        if(hp->index_s1==1){result=-G3(hp);}
                        if(hp->index_s1==2){result=-G3_star(hp);}
                        if(hp->index_s1==3){result=G1(hp);}
		}
	}

	if(hp->flag_anni==0 && std::abs(hp->Jz)<2){result+=anni_spin(hp);}

	return result;
}


PetscReal interaction(PetscErrorCode ierr,params *params, h_params *hp)
{
	PetscReal h_elem;
	PetscReal ct_term;
	PetscReal jacobian;
	PetscReal aa,denom;
	PetscReal m1s,m2s;
	PetscReal mu1,mu2,wmu1,wmu2,wth1,wth2,x1,x2,k1,k2;

	m1s=hp->m1*hp->m1;
	m2s=hp->m2*hp->m2;
        mu1 = hp->mu1;
        mu2 = hp->mu2;

	wmu1=hp->wmu1;
	wmu2=hp->wmu2;
	wth1=hp->wth1;
	wth2=hp->wth2;

        x1 = hp->x1;
        k1 = hp->k1;
        x2 = hp->x2;
        k2 = hp->k2;
	
	h_elem=0.0;

	hp->omega = ((m1s+k1*k1)/x1/(1.0-x1)+(m2s+k2*k2)/x2/(1.0-x2))/2.0;
	if(hp->index_all1==hp->index_all2)
	{
		VecGetValues(params->CT,1,&hp->index_all1,&ct_term);

		h_elem+=ct_term;
	
	}else{
		if(hp->index_f1==hp->index_f2 && hp->index_t1==hp->index_t2 && hp->index_m1==hp->index_m2){ return 0.0;}

		aa = (x1-x2)*(x1-x2)*m1s/2.0*(1.0/(1.0-x2)/(1.0-x1)+1.0/x1/x2) + k1*k1 + k2*k2+ (x1-x2)/2.0*(k1*k1*(1.0/(1.0-x1)-1.0/x1)-k2*k2*(1.0/(1.0-x2)-1.0/x2));

                denom = aa*aa-4.0*k1*k1*k2*k2;
                hp->A  = 1.0/sqrt(denom);
                hp->B  = (1.0-aa*hp->A)/2.0;

		jacobian = sqrt(wmu1*wmu2*wth1*wth2)*mu2*mu1*sqrt( 4.0*x1*(1.0-x1)*x2*(1.0-x2)/sqrt(m1s + mu1*mu1)/sqrt(m2s + mu2*mu2));
		
		h_elem=jacobian/M_PI*general_spin(hp);

	}	
	return h_elem;
}
