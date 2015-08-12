struct params{

	char ofile_n[100];

        PetscInt nf;
	PetscInt Jz;
	PetscInt N_tot;
        PetscReal alpha;
	PetscInt flag_nm_dis;
	PetscInt flag_nt_dis;
	PetscInt flag_anni,flag_mix,flag_asy;

	std::string sflag_nm_dis;
	std::string sflag_nt_dis;
	std::string hfile;

        std::vector<PetscReal> m;
        std::vector<PetscInt> nm;
        std::vector<PetscInt> nt;
        std::vector<PetscInt> N_tot_f;
        std::vector<PetscInt> nm_tot;
        std::vector<PetscInt> nt_tot;
        std::vector<PetscReal> lambda;
	std::vector<PetscReal> p_bohr;

	Vec mu;
	Vec wmu;

	Vec theta;
	Vec wtheta;

	Vec CT;
};

struct int_params{

	PetscInt Jzc;
	PetscInt nf;
	PetscInt N_tot;	
	PetscInt flag_asy;
	PetscReal alpha;

        std::vector<PetscReal> m;
        std::vector<PetscInt> nm;
        std::vector<PetscInt> nt;
        std::vector<PetscReal> lambda;
        std::vector<PetscReal> p_bohr;

        PetscReal mu1,mu2,th1,th2,x1,x2,k1,k2,mf,wth1,wth2,wmu1,wmu2;
        PetscInt index_all,index_m, index_t, index_f,index_s;

        PetscReal epsabs,epsrel;
        size_t N_work;

};

struct h_params{

	PetscReal alpha,omega,mu1,mu2,th1,th2,x1,x2,k1,k2,m1,m2,wth1,wth2,wmu1,wmu2;
        PetscInt index_all1,index_m1, index_t1, index_f1,index_s1;
        PetscInt index_all2,index_m2, index_t2, index_f2,index_s2;
	PetscInt Jz;
	PetscReal A,B;
	PetscInt flag_perm,flag_anni,flag_mix,flag_asy;
};
