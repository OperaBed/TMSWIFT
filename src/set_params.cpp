void get_index(PetscInt index, int_params *int_params, params *params)
{
        PetscInt ct,ft;
        ft=0;
//        int_params->index_s=(index-index%params->N_tot)/params->N_tot;
        ct=index;
        while(ct>=(params->N_tot_f[ft]*4)){ft++;}
        int_params->index_f=ft;
        if(ft!=0){ct=ct-params->N_tot_f[ft-1]*4;}
        int_params->index_m=(ct-ct%(4*params->nt[ft]))/(4*params->nt[ft]);
	ct=ct-int_params->index_m*(4*params->nt[ft]);
        int_params->index_t=(ct-ct%4)/4;
	ct=ct-4*int_params->index_t;
	int_params->index_s=ct;
        int_params->index_all=index;
}
