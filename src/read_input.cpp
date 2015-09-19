void read_input(char **argv, params *params)
{

        std::string ifile_n;

        ifile_n = argv[1];
        params->ofile_n.append(ifile_n.c_str());
        get_params_from_file(ifile_n,params);
        set_sflags(params);

        params->N_tot=0;

        for(int i=0;i<params->nf;i++)
        {
                params->p_bohr.push_back(params->m[i]*params->alpha/2.0);
                params->N_tot_f.push_back(params->N_tot);
                params->N_tot+=params->nm[i]*params->nt[i];
                params->N_tot_e.push_back(params->nm[i]*params->nt[i]);
        }
                params->N_tot_f.push_back(params->N_tot);
}

