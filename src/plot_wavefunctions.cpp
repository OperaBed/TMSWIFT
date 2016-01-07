#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <slepceps.h>

#include "params.cpp"
#include "set_flags.cpp"
#include "print_params.cpp"
#include "string.cpp"
#include "read_input.cpp"
#include "read_file.cpp"
static char help[] = "Plotting code for extracting single wave functions and computing probablilties for the problem of QFT bound states on the light front.\n\n";

double round_to_digits(double value, int digits)
{
    if (value == 0.0) // otherwise it will return 'nan' due to the log10() of zero
        return 0.0;

    double factor = pow(10.0, digits - ceil(log10(fabs(value))));
    return round(value * factor) / factor;   
}

int main(int argc, char *argv[])
{
	std::vector<double> min;
	std::vector<double> max;
	std::vector<double> increm;
	std::string f_name;
	std::string line;
	int nev;
	int s;
	PetscErrorCode ierr;

        SlepcInitialize(&argc,&argv,(char*)0,help);

	std::ifstream ifile;
	std::ofstream gfile;
	params p;

        ierr = PetscPrintf(PETSC_COMM_WORLD,"----------------------------------\n");CHKERRQ(ierr);
        read_input(argv,&p);
        print_input(ierr,&p);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"----------------------------------\n");CHKERRQ(ierr);

	nev=std::atoi(argv[2]);
	min.resize(12*p.nf);
	max.resize(12*p.nf);
	increm.resize(12*p.nf);
	for(int i=0;i<12*p.nf;i++)
	{
		min[i]=0.0;
		max[i]=0.0;
		increm[i]=0.0;
	}

	for(int l=0;l<p.nf;++l)
	{	
		for(int m=0;m<4;++m)
		{	
			std::stringstream ss;
			ss << p.ofile_n << "n" << nev << "_f" << l << "_s" << m;
			ss >> f_name;	
		
//			std::cout << f_name << std::endl;
			ifile.open(f_name.c_str());	
	        	while (std::getline(ifile, line))
        		{
                		std::istringstream iss(line);
                		double x,k,a,prob;
                		if ((iss >> x >> k >> a >> prob))
                		{
					if(x<min[12*l+3*m+0]) min[12*l+3*m+0]=x;	       
                 			if(k<min[12*l+3*m+1]) min[12*l+3*m+1]=k;	       
                 			if(prob*prob<min[12*l+3*m+2]) min[12*l+3*m+2]=prob*prob;

					if(x>max[12*l+3*m+0]) max[12*l+3*m+0]=x;
                                        if(k>max[12*l+3*m+1]) max[12*l+3*m+1]=k;
                                        if(prob*prob>max[12*l+3*m+2]) max[12*l+3*m+2]=prob*prob;	       
                		}
        		}
			ifile.close();
		}
	}

	for(int i=0;i<12*p.nf;i++)
        {
        	increm[i]=round_to_digits((max[i]-min[i])/5.0,2);
//	        std::cout << i << " " << round_to_digits(min[i],2) << " " << round_to_digits(max[i],2) << " " << round_to_digits(increm[i],2) << std::endl;
        }
	if(p.Jz==0){s=2;}else{s=4;}
	gfile.open("plot_wave",std::ofstream::trunc);

	gfile << "set term postscript eps enhanced color dashed font \"Times-Roman\" size "<< 3.75*s << "," << 10.0/3.0*p.nf  << "\nset output ";
	gfile << "'wave_"<< nev  << ".eps'\nset hidden3d\nset border 1023-128\nset ztics\nset grid z\nset ticslevel 0\nset xlabel 'x'\nset ylabel 'k_{/Symbol \\136}'\nset format z \"%.1tx10^{%T}\"\nunset key\n";
	gfile << "set xrange [0:1]\nset multiplot layout " << p.nf << "," << s << std::endl;

       for(int i=0;i<4*p.nf;i++)
        {	
		int f;
		std::stringstream ss;
		if(i<4){f=0;}else if(i<8){f=1;}else if(i<12){f=2;}
                ss << p.ofile_n << "n" << nev << "_f" << f  << "_s" << i%4;
                ss >> f_name;

		if(i%4<2)
		{
			gfile << "tit = \'";
			if(i<4){gfile << "{/Symbol mm} ";}
			else if(i<8) {gfile << "ee ";}
			else if(i<12) {gfile << "{/Symbol tt} ";}
			if(i%4==0){gfile << "{/Symbol \\255 }{/Symbol \\255 } ";}
			if(i%4==1){gfile << "{/Symbol \\255 }{/Symbol \\257 } ";}
			gfile << "component of state " << nev << " for J_{z}=" << p.Jz << "\'\n";
			gfile << "set title font \'Times-Roman,20\'\nset title tit\n";
			gfile << "set xtics 0.1,.2,.9 offset -2\n";
			gfile << "set ytics " << round_to_digits(min[1+3*i],2) << "," << round_to_digits(increm[1+3*i],2) << "," << round_to_digits(max[1+3*i],2) <<" offset 2\n";
			gfile << "set ztics " << round_to_digits(min[2+3*i],2) << "," << round_to_digits(increm[2+3*i],2) << "," << round_to_digits(max[2+3*i],2) <<" offset -1\n";
			gfile << "splot \'" << f_name << "\' u 1:2:($4*$4) w lines\n";
		}
		if(p.Jz!=0 && i%4>1)
		{
			gfile << "tit = \'";
                        if(i<4){gfile << "{/Symbol mm} ";}
                        else if(i<8) {gfile << "ee ";}
                        else if(i<12) {gfile << "{/Symbol tt} ";}
			if(i%4==2){gfile << "{/Symbol \\257 }{/Symbol \\255 } ";}
			if(i%4==3){gfile << "{/Symbol \\257 }{/Symbol \\257 } ";}		
			gfile << "component of state " << nev << " for J_{z}=" << p.Jz << "\'\n";
			gfile << "set title font \'Times-Roman,20\'\nset title tit\n";
                        gfile << "set xtics 0.1,.2,.9 offset -2\n";
                        gfile << "set ytics " << round_to_digits(min[1+3*i],2) << "," << round_to_digits(increm[1+3*i],2) << "," << round_to_digits(max[1+3*i],2) <<" offset 2\n";
                        gfile << "set ztics " << round_to_digits(min[2+3*i],2) << "," << round_to_digits(increm[2+3*i],2) << "," << round_to_digits(max[2+3*i],2) <<" offset -1\n";
			gfile << "splot \'./" << f_name << "\' u 1:2:($4*$4) w lines\n";
		}

	}	
	gfile.close();	
	s=std::system("gnuplot plot_wave");
//	s=std::system("rm plot_wave");
        ierr = SlepcFinalize();

	return 0;
}
