#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
/*
0 = number of flavours 
1 = J_z
2 = alpha
3 = N_mu_dis
4 = N_th_dis
5 = Flag_Anni
6 = Flag_Mix_Flavours
7 = Flag_Asy_G2
8 = Flag_Mumps
9,13,17,21... = mass of flavours
10,14,18,22...= N_mu
11,15,19,23...= N_th
12,16,20,24...= Lambda
*/

void write_file(std::string filename,std::string outfilename, std::vector<double> var)
{
	std::ofstream ofile;
	int counter;
	counter = 0;
        ofile.open(filename.c_str());

	ofile << outfilename << "	#outputfile" << std::endl;
	ofile << (int)var[0] <<	"	#number_of_flavors" << std::endl;
	ofile << (int)var[1] << "	#J_z" << std::endl;
	ofile << var[2] << "	#alpha" << std::endl;
	ofile << (int)var[3] << "	#N_mu_dis" << std::endl;
	ofile << (int)var[4] << "	#N_theta_dis" << std::endl;
	ofile << (int)var[5] << "	#Flag_Annihilation" << std::endl;
	ofile << (int)var[6] << "	#Flag_Mix_Flavors" << std::endl;
	ofile << (int)var[7] << "	#Flag_Asy_G2" << std::endl;
	ofile << (int)var[8] << "	#Flag_MUMPS" << std::endl;
	ofile << "---------------------" << std::endl;
	for(int i=0;i<var[0];i++)
	{
		ofile << var[4*counter+1+8] << "	#mass_" << i << std::endl;
		ofile << var[4*counter+2+8] << "	N_mu_" << i <<std::endl;
		ofile << var[4*counter+3+8] << "	N_theta_" << i <<std::endl;
		ofile << var[4*counter+4+8] << "	#Lambda_" << i << std::endl;
		ofile << "---------------------" << std::endl;
		counter++;
	}

	ofile.close();
}

int main(int argc, char *argv[])
{

	std::string outfilename;
	std::string runfilename;
	std::string filename;
	std::ofstream run;
	std::vector<double> var;
	int counter = 0;
	filename="init_";
	runfilename="run_file";
	outfilename="OUT/truem";
	
	run.open(runfilename.c_str());
	var.resize(21);

	var[0]=1;//flavors
	var[1]=0;//J_z
	var[2]=0.3;//alpha
	var[3]=1;//N_mu_dis
	var[4]=1;//N_th_dis
	var[5]=0;//flag anni
	var[6]=0;//Flag mix flavors
	var[7]=1;//Flag Asy g2
	var[8]=1;//Flag mumps

	auto var_alpha = std::vector<double> {0.01,0.07,0.1};
	auto var_dis = std::vector<double> {1};
	auto var_asy = std::vector<double> {0};

	auto var_m = std::vector<std::vector<double>>
	{
		{1.0},
		{0.5},
		{2.0}
	};
	
	auto var_nm = std::vector<std::vector<double>>
	{
		{13},
		{0},
		{0}
	};
	
	auto var_lam = std::vector<std::vector<double>>
	{
		{5,10,20},
		{0},
		{0}
	};

	
	for (auto i_alpha=var_alpha.begin();i_alpha!=var_alpha.end();i_alpha++){
		var[2]=*i_alpha;
	for (auto i_dis=var_dis.begin();i_dis!=var_dis.end();i_dis++){
		var[3]=*i_dis;
		var[4]=*i_dis;
	for (auto i_asy=var_asy.begin();i_asy!=var_asy.end();i_asy++){
		var[7]=*i_asy;
	for (auto i_m=var_m[0].begin();i_m!=var_m[0].end();i_m++){
		var[9]=*i_m;
	for (auto j_m=var_m[1].begin();j_m!=var_m[1].end();j_m++){
		var[13]=*j_m;
	for (auto k_m=var_m[2].begin();k_m!=var_m[2].end();k_m++){
		var[17]=*k_m;
	for (auto i_nm=var_nm[0].begin();i_nm!=var_nm[0].end();i_nm++){
		var[10]=*i_nm;
		var[11]=*i_nm;
	for (auto j_nm=var_nm[1].begin();j_nm!=var_nm[1].end();j_nm++){
		var[14]=*j_nm;
		var[15]=*j_nm;
	for (auto k_nm=var_nm[2].begin();k_nm!=var_nm[2].end();k_nm++){
		var[18]=*k_nm;
		var[19]=*k_nm;
	for (auto i_lam=var_lam[0].begin();i_lam!=var_lam[0].end();i_lam++){
		var[12]=*i_lam;
	for (auto j_lam=var_lam[1].begin();j_lam!=var_lam[1].end();j_lam++){
		var[16]=*j_lam;
	for (auto k_lam=var_lam[2].begin();k_lam!=var_lam[2].end();k_lam++){
		var[20]=*k_lam;
		counter++;
		filename="init_";
		filename.append(std::to_string(counter));		
		//run << "~/Software/petsc-3.6.0/arch-linux2-c-opt/bin/mpiexec -np 4 ./tmswift " << filename << std::endl;
		run << "mpiexec -np 2 ./tmswift " << filename << std::endl;
		write_file(filename,outfilename,var);
	}}}}}}}}}}}}
	std::cout << "------- " << counter << " input files produced --------" << std::endl;
	return 0;
}
