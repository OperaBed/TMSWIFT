#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <slepceps.h>

#include "params.cpp"
#include "set_flags.cpp"
#include "print_params.cpp"
#include "string.cpp"
#include "read_input.cpp"
#include "read_file.cpp"
static char help[] = "Analysis code for extracting single wave functions and computing probablilties for the problem of QFT bound states on the light front.\n\n";

int main(int argc, char *argv[])
{
	std::vector<double> x;
	std::vector<double> k;
	std::vector<double> asy;
	std::vector<double> evec_r;
	std::vector<double> eval;
	std::vector<double> nm;
	std::vector<double> nt;
	std::vector<double> n_tot;
	std::string line;
	int nev;
	int index;
	int index_p;
	PetscErrorCode ierr;

	std::ofstream f0,f1,f2,f3;
	std::vector<std::string> f_name;
	std::vector<double> prob;
	double k_exp;
	double k2_exp;
	double dc_scalar;
	double dc_pseudo;
	params p;

	SlepcInitialize(&argc,&argv,(char*)0,help);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"----------------------------------\n");CHKERRQ(ierr);

        read_input(argv,&p);
	print_input(ierr,&p);	
        ierr = PetscPrintf(PETSC_COMM_WORLD,"----------------------------------\n");CHKERRQ(ierr);
	f_name.resize(4);
	prob.resize(p.nf);

	line=p.ofile_n;
	line.append("x");
	read_ofile(line,&x);
      
  	line=p.ofile_n;
        line.append("k");
	read_ofile(line,&k);

        line=p.ofile_n;
        line.append("asy");
	read_ofile(line,&asy);

        line=p.ofile_n;
        line.append("evecr");
	read_ofile(line,&evec_r);

	line=p.ofile_n;
	line.append("eval");
	read_ofile(line,&eval);

	nev=evec_r.size()/x.size();

	for(int i=0;i<nev;++i)
	{
		dc_scalar=0;
		dc_pseudo=0;
		k_exp=0;
		k2_exp=0;

		for(int o=0;o<p.nf;++o)
		{
			prob[o]=0;
		}

		for(int l=0;l<p.nf;++l)
		{	
			for(int m=0;m<4;++m)
			{	
				std::stringstream ss;
				ss << p.ofile_n << "n" << i << "_f" << l << "_s" << m;
				ss >> f_name[m];	
			}
			f0.open(f_name[0].c_str());	
			f1.open(f_name[1].c_str());	
			f2.open(f_name[2].c_str());	
			f3.open(f_name[3].c_str());	
			for(int j=0;j<p.N_tot_e[l];++j)
			{
				index=4*j+i*x.size()+4*p.N_tot_f[l];
				index_p=4*j+4*p.N_tot_f[l];
				for(int o=0; o<4; o++){
					prob[l]+=evec_r[o+index]*evec_r[o+index];
					k_exp+=k[0+index_p]*evec_r[o+index]*evec_r[o+index];
					k2_exp+=k[0+index_p]*k[0+index_p]*evec_r[o+index]*evec_r[o+index];
					dc_scalar+=1.0/asy[0+index_p]*(evec_r[1+index]-evec_r[2+index]);
					dc_pseudo+=1.0/asy[0+index_p]*(evec_r[1+index]+evec_r[2+index]);
				}
				f0 << x[0+index_p] << " " << k[0+index_p] << " " << asy[0+index_p] << " " << evec_r[0+index] << std::endl;
				f1 << x[1+index_p] << " " << k[1+index_p] << " " << asy[1+index_p] << " " << evec_r[1+index] << std::endl;
				f2 << x[2+index_p] << " " << k[2+index_p] << " " << asy[2+index_p] << " " << evec_r[2+index] << std::endl;
				f3 << x[3+index_p] << " " << k[3+index_p] << " " << asy[3+index_p] << " " << evec_r[3+index] << std::endl;
				if((j+1)%p.nt[l]==0)
				{
					f0 << std::endl;
					f1 << std::endl;
					f2 << std::endl;
					f3 << std::endl;
				}
			}
			f0.close();
			f1.close();
			f2.close();
			f3.close();
		}
		std::cout << std::setw(5) << i << ": " << std::setw(10) << eval[i] << " -- ";
		for(int o=0;o<p.nf;++o)
		{
			std::cout << std::setw(12) << prob[o] << "   ";
		}
		if(std::fabs(dc_scalar)<1e-10){dc_scalar=0;}
		if(std::fabs(dc_pseudo)<1e-10){dc_pseudo=0;}
		std::cout << " f_S: " << std::setw(12) << std::fabs(dc_scalar) << " f_P: " << std::setw(12) << std::fabs(dc_pseudo);
		std::cout << " <k>: " << std::setw(12) << k_exp << " <k^2>: " << std::setw(12) << k2_exp;
		std::cout << std::endl;
	}

	        ierr = SlepcFinalize();

	return 0;
}
