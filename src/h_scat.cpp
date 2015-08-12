PetscReal G1(h_params *params)
/* Scattering Channel, Diagonal Element, Parallel Spins */
{
        PetscReal h_elem;
        PetscReal m,x1,x2,k1,k2;
        PetscInt n;

        set_hel_params(m,x1,x2,k1,k2,n,params);

        h_elem = m*m*(1.0/x1/x2+1/(1-x2)/(1-x1))*Int(std::abs(1-n),params)+ k1*k2/x2/x1/(1-x2)/(1-x1)*Int(std::abs(n),params);
        return h_elem;
}

PetscReal G2(h_params *params)
/* diagonal matrix element, anti-parallel spins */
{
        PetscReal h_elem;
        PetscReal m,x1,x2,k1,k2;
        PetscInt n;

        set_hel_params(m,x1,x2,k1,k2,n,params);

        h_elem = 1/x1/x2 + 1.0/(1-x2)/(1-x1);
        h_elem = (m*m*h_elem + k2*k2/x2/(1-x2)+ k1*k1/x1/(1-x1))*Int(std::abs(n),params)+k1*k2*(Int(std::abs(1-n),params)/x1/x2+Int(std::abs(1+n),params)/(1-x1)/(1-x2));
        if(n==0 && params->flag_asy==0){h_elem+=2.0/(x1+x2-2.0*x1*x2);}

        return h_elem;
}

PetscReal G3(h_params *params)
{
        PetscReal h_elem;
        PetscReal m,x1,x2,k1,k2;
        PetscInt n;

        set_hel_params(m,x1,x2,k1,k2,n,params);

	h_elem = -m/x1/x2*(k1*Int(std::abs(1-n),params)-k2*(1-x1)/(1-x2)*Int(std::abs(n),params));

        return h_elem;
}

PetscReal G3_star(h_params *params)
{
        PetscReal h_elem;
        PetscReal m,x1,x2,k1,k2;
        PetscInt n;

        set_hel_params(m,x1,x2,k1,k2,n,params);

        h_elem = m/(1.0-x1)/(1.0-x2)*(k1*Int(std::abs(1-n),params)-k2*x1/x2*Int(std::abs(n),params));

        return h_elem;
}

PetscReal G4(h_params *params)
{
        PetscReal h_elem;
        PetscReal m,x1,x2,k1,k2;
        PetscInt n;

        set_hel_params(m,x1,x2,k1,k2,n,params);

        h_elem = -m*m*(x1-x2)*(x1-x2)/(1-x1)/(1-x2)/x1/x2*Int(std::abs(n),params);

        return h_elem;
}

