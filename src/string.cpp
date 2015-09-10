std::string toStrBase82(unsigned long long num) {
   std::string charset = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz+-.,?<>;:=!@#$^&*()~`";
   int base = charset.length();
   std::string str = num ? "" : "0";

   while (num) {
      str = charset.substr(num % base, 1) + str;
      num /= base;
   }

   return str;
}

void get_outfile_name(params *p)
{
	unsigned long long num;

	std::string string_temp;

	num=p->nf*10000000000+p->Jz*1000000000+std::round(p->alpha*100)*1000000+p->flag_nm_dis*100000+p->flag_nt_dis*10000+p->flag_anni*1000+p->flag_mix*100+p->flag_asy*10+p->flag_mumps;
                string_temp.append(toStrBase82(num));

	for(PetscInt i=0;i<p->nf;i++)
	{
		num=std::round(p->m[i]*10)*10000000000+p->nm[i]*10000000+p->nt[i]*10000+std::round(p->lambda[i]*10);
		string_temp.append("_");   
		string_temp.append(toStrBase82(num));	
	}
	p->ofile_n.append(string_temp);
}
