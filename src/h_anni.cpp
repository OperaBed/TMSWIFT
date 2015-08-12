PetscReal F1(h_params *params)
{
	PetscReal h_elem,m,x1,x2,k1,k2;
        PetscInt n;

	h_elem=0.0;

        set_hel_params(m,x1,x2,k1,k2,n,params);

	if(std::abs(n)==1)
	{
		h_elem = 2.0*params->m1*params->m2/params->omega*(1.0/x1+1/(1-x1))*(1/x2+1.0/(1-x2));
	}
	
	return h_elem;
}

PetscReal F2(h_params *params)
{
        PetscReal h_elem,m,x1,x2,k1,k2;
        PetscInt n;

        h_elem=0.0;

        set_hel_params(m,x1,x2,k1,k2,n,params);
	if(n==0){h_elem=4.0;}
        if(std::abs(n)==1){h_elem = 2.0/params->omega*k1*k2/x1/x2;}

        return h_elem;
}

PetscReal F2_star(h_params *params)
{
        PetscReal h_elem,m,x1,x2,k1,k2;
        PetscInt n;

        h_elem=0.0;

        set_hel_params(m,x1,x2,k1,k2,n,params);
        if(n==0){h_elem=4.0;}
        if(std::abs(n)==1){h_elem = 2.0/params->omega*k1*k2/(1.0-x1)/(1.0-x2);}

        return h_elem;
}

PetscReal F3(h_params *params)
{
        PetscReal h_elem,m,x1,x2,k1,k2;
        PetscInt n;

        h_elem=0.0;

        set_hel_params(m,x1,x2,k1,k2,n,params);
        if(std::abs(n)==1){h_elem = 2.0*params->m2/params->omega*k1/(1-x1)*(1/x2+1.0/(1-x2));}

        return h_elem;
}


PetscReal F3_star(h_params *params)
{
        PetscReal h_elem,m,x1,x2,k1,k2;
        PetscInt n;

        h_elem=0.0;

        set_hel_params(m,x1,x2,k1,k2,n,params);
        if(std::abs(n)==1){h_elem = -2.0*params->m2/params->omega*k1/x1*(1/x2+1.0/(1-x2));}

        return h_elem;
}

PetscReal F4(h_params *params)
{
        PetscReal h_elem,m,x1,x2,k1,k2;
        PetscInt n;

        h_elem=0.0;

        set_hel_params(m,x1,x2,k1,k2,n,params);
        if(n==0){h_elem=4.0;}
        if(std::abs(n)==1){h_elem = -2.0/params->omega*k1*k2/x1/(1-x2);}

        return h_elem;
}

PetscReal anni_spin(h_params *hp)
{
        PetscReal result;

        result=0.0;

        if(hp->Jz==-1)
        {
                hp->index_s1=3-hp->index_s1;
                hp->index_s2=3-hp->index_s2;
        }

        if(hp->index_s1==3 || hp->index_s2==3){return 0.0;}
        else{
                if(hp->index_s2==0)
                {
                        if(hp->index_s1==0){result=F1(hp);}
                        if(hp->index_s1==1){result=F3(hp);}
                        if(hp->index_s1==2){result=F3_star(hp);}
                }

                if(hp->index_s2==1)
                {
                        if(hp->index_s1==0){hp->flag_perm=1;result=F3_star(hp);}
                        if(hp->index_s1==1){result=F2_star(hp);}
                        if(hp->index_s1==2){result=F4(hp);}
                }

                if(hp->index_s2==2)
                {
                        if(hp->index_s1==0){hp->flag_perm=1;result=F3_star(hp);}
                        if(hp->index_s1==1){hp->flag_perm=1;result=F4(hp);}
                        if(hp->index_s1==2){result=F2(hp);}
                }

        }

        return result;
}

