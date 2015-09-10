//Found in input.cpp
void set_sflags(params *params);
void read_input(PetscErrorCode ierr, char **argv, params *params);
void print_input(PetscErrorCode ierr, params *params);

//Found in string.cpp
std::string toStr(unsigned long long num);

//Found in discretize.cpp
void discretize(PetscErrorCode ierr, params *params);

//Found in set_params.cpp
void get_index(PetscInt index, int_params *int_params, params *params);
void set_params(params *params, int_params *iparams);

//Found in coulomb_cont.cpp
void coulomb_trick(PetscErrorCode ierr, params *params);
PetscReal ct_function(PetscReal theta_local, void *params);
PetscReal ct_integrand(PetscReal mu_local, void *params);

//Found in coulomb_discrete.cpp
void ct_discrete(PetscErrorCode ierr, params *params);

//Found in counterterm.cpp
PetscReal self_counterterm(params *params, h_params *hp);

//Found in physics.cpp
PetscReal Int(PetscInt n, h_params *hp);
void set_hel_params(PetscReal &m,PetscReal &x1,PetscReal &x2,PetscReal &k1,PetscReal &k2,PetscInt &n,h_params *params);
PetscReal interaction(PetscErrorCode ierr, params *params, h_params *hp);
PetscReal general_spin(h_params *hp);

//Found in h_scat.cpp
PetscReal G1(h_params *params);
PetscReal G2(h_params *params);
PetscReal G3(h_params *params);
PetscReal G3_star(h_params *params);
PetscReal G4(h_params *params);

//Found in h_anni.cpp
PetscReal F1(h_params *params);
PetscReal F2(h_params *params);
PetscReal F2_star(h_params *params);
PetscReal F3(h_params *params);
PetscReal F3_star(h_params *params);
PetscReal F4(h_params *params);
PetscReal anni_spin(h_params *hp);

//Found in hamiltonian.cpp
void hamiltonian(PetscErrorCode ierr, params *params, Mat &H);

//Found in solver.cpp
void eigensolver(PetscErrorCode ierr, params *params, Mat &H, int argc, char **argv);

//Found in output.cpp
void write_output(PetscErrorCode ierr, params *params);
void print_progress_header(PetscErrorCode ierr);
//Found in cleanup.cpp
void cleanup(PetscErrorCode ierr, params *params);
