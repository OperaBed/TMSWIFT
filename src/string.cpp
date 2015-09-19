unsigned long long uintpow(unsigned value, unsigned exp) {
        unsigned long long result = 1;
        for(unsigned int i=0; i<exp; ++i) {
                result *= value;
        }
        return result;
}


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

unsigned long long frStrBase(std::string str)
{
        unsigned long long val;
        unsigned int i,j;
        std::string charset = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz+-.,?<>;:=!@#$^&*()~`";
        j=0;
        val=0;
        while(j<str.length())
        {
        i=0;
        while (str.substr(j,1)!=charset.substr(i,1))
        {
                i++;
        }
        val+=uintpow(charset.length(),str.length()-(j+1))*i;
        j++;
        }
        return val;
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
	p->ofile_n.append("_");
	p->ofile_n.append(string_temp);
}

void get_params_gen(unsigned long long num, params *p)
{
        unsigned long long b10;
        b10=1e10;
        for(int i=0;i<2;i++)
        {
                if(i==0){p->nf=(num-num%b10)/b10;}
                if(i==1){p->Jz=(num-num%b10)/b10;}
                num=num%b10;
                b10=b10/10;
        }
        p->alpha=1.0*(num-num%1000000)/100000000;
        num=num%1000000;
        b10=1e5;
        for(int i=0;i<6;i++)
        {
                if(i==0){p->flag_nm_dis=(num-num%b10)/b10;}
                else if(i==1){p->flag_nt_dis=(num-num%b10)/b10;}
                else if(i==2){p->flag_anni=(num-num%b10)/b10;}
                else if(i==3){p->flag_mix=(num-num%b10)/b10;}
                else if(i==4){p->flag_asy=(num-num%b10)/b10;}
                else if(i==5){p->flag_mumps=(num-num%b10)/b10;}
                num=num%b10;
                b10=b10/10;
        }
}


void get_params_flavor(unsigned long long num, params *p)
{
        unsigned long long b10;
        b10=1e10;
        for(int i=0;i<4;i++)
        {
                if(i==0){p->m.push_back(1.0/10.0*(num-num%b10)/b10);}
		else if(i==1){p->nm.push_back((num-num%b10)/b10);}
		else if(i==2){p->nt.push_back((num-num%b10)/b10);}
                else if(i==3){p->lambda.push_back(num/10.0);}
                num=num%b10;
                b10=b10/1e3;
        }

}

void get_params_from_file(std::string str, params *p)
{
        unsigned int i=0;
        int start=0;
        int end=0;
        int flag=0;
        while(i<str.length())
        {
                if(str.substr(i,1)=="_")
                {
                        start=end;
                        end=i;
                        if(start!=0)
                        {
                                if(flag==0){get_params_gen(frStrBase(str.substr(start+1,end-start-1)),p); flag++;}
                                else{get_params_flavor(frStrBase(str.substr(start+1,end-start-1)),p);}
                        }
                }
                i++;
        }
}


