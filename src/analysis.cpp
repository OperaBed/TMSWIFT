#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
struct params{

        char ofile_n[100];

        int nf;
        int Jz;
        int N_tot;

        std::vector<double> m;
        std::vector<int> nm;
        std::vector<int> nt;
        std::vector<int> N_tot_f;
        std::vector<int> N_tot_e;

};

void read_input(char **argv, params *params)
{

        std::ifstream ifile;
        std::string ifile_n;
        std::string null;

        double datar;
        int  datai;

        ifile_n = argv[1];
        ifile.open(ifile_n.c_str());

        ifile >> params->ofile_n >> null;
        ifile >> params->nf >> null;
        ifile >> params->Jz >> null;
        ifile >> null >> null;
        ifile >> null >> null;
        ifile >> null >> null;
        ifile >> null >> null;
        ifile >> null >> null;
        ifile >> null >> null;
        for(int i=0; i<params->nf; i++)
        {
                ifile >> null;

                ifile >> datar >> null;
                params->m.push_back(datar);
                ifile >> datai >> null;
                params->nm.push_back(datai);
                ifile >> datai >> null;
                params->nt.push_back(datai);
                ifile >> null >> null;
        }
        ifile.close();

        params->N_tot=0;

        for(int i=0;i<params->nf;i++)
        {
                params->N_tot_f.push_back(params->N_tot);
                params->N_tot+=params->nm[i]*params->nt[i];
                params->N_tot_e.push_back(params->nm[i]*params->nt[i]);
        }
                params->N_tot_f.push_back(params->N_tot);
}

void read_ofile(std::string file_name, std::vector<double> *vec)
{
	std::string line;
	std::fstream infile;
        infile.open(file_name.c_str());
        while (std::getline(infile, line))
        {
                std::istringstream iss(line);
                double a;
                if ((iss >> a)) 
                {
			vec->push_back(a);
		}
        }
        infile.close();
}



int main(int argc, char *argv[])
{
	std::vector<double> x;
	std::vector<double> k;
	std::vector<double> asy;
	std::vector<double> evec_r;
	std::vector<double> nm;
	std::vector<double> nt;
	std::vector<double> n_tot;
	std::string line;
	int nev;
	int index;
	int index_p;
	

	std::ofstream f0,f1,f2,f3;
	std::vector<std::string> f_name;
	std::vector<double> prob;
	params p;

        read_input(argv,&p);
	
	f_name.resize(4);
	prob.resize(p.nf);

	line=p.ofile_n;
	line.append("_x");
	read_ofile(line,&x);
      
  	line=p.ofile_n;
        line.append("_k");
	read_ofile(line,&k);

        line=p.ofile_n;
        line.append("_asy");
	read_ofile(line,&asy);

        line=p.ofile_n;
        line.append("_evecr");
	read_ofile(line,&evec_r);

	nev=evec_r.size()/x.size();
	std::cout << x.size() << " " << k.size() << " " << asy.size() << " " << evec_r.size() << std::endl;
	std::cout << p.N_tot_f[0] << " " << p.N_tot_f[1] << " " << p.N_tot_f[2] << std::endl;
	std::cout << p.N_tot_e[0] << " " << p.N_tot_e[1]  << std::endl;
	for(int i=0;i<nev;++i)
	{
		for(int o=0;o<p.nf;++o)
		{
			prob[o]=0;
		}

		for(int l=0;l<p.nf;++l)
		{	
			for(int m=0;m<4;++m)
			{	
				std::stringstream ss;
				ss << p.ofile_n << "_n" << i << "_f" << l << "_s" << m;
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
				for(int o=0; o<4; o++){prob[l]+=evec_r[o+index]*evec_r[o+index];}
				f0 << x[0+index_p] << " " << k[0+index_p] << " " << asy[0+index_p] << " " << evec_r[0+index] << std::endl;
				f1 << x[1+index_p] << " " << k[1+index_p] << " " << asy[1+index_p] << " " << evec_r[1+index] << std::endl;
				f2 << x[2+index_p] << " " << k[2+index_p] << " " << asy[2+index_p] << " " << evec_r[2+index] << std::endl;
				f3 << x[3+index_p] << " " << k[3+index_p] << " " << asy[3+index_p] << " " << evec_r[3+index] << std::endl;
				if((j+1)%p.nm[l]==0)
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
		std::cout << i << ": ";
		for(int o=0;o<p.nf;++o)
		{
			std::cout << prob[o] << "   ";
		}
		std::cout << std::endl;
	}
	return 0;
}
