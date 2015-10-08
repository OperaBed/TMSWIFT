void set_sflags(params *params)
{

	switch (params->flag_nm_dis)
	{
		case 1:
		params->sflag_nm_dis="legendre";
		break;

		case 2:
		params->sflag_nm_dis="clenshaw_curtis";
		break;

		case 4:
		params->sflag_nm_dis="laguerre";
		break;

		case 3:
		params->sflag_nm_dis="chebyshev1";
		break;

		case 5:
		params->sflag_nm_dis="uniform";
		break; 
		
		default:
		params->sflag_nm_dis="legendre";
		break;

	}

        switch (params->flag_nt_dis)
        {
                case 1:
                params->sflag_nt_dis="legendre";
                break;

                case 2:
                params->sflag_nt_dis="clenshaw_curtis";
                break;

                case 3:
                params->sflag_nt_dis="chebyshev1";
                break;

                default:
                params->sflag_nt_dis="legendre";
                break;

        }


}


