/*--------------------------------------------------------------------------------------*/
/*		 _______ __  __  _______          _______ ______ _______ 		*/
/*		|__   __|  \/  |/ ____\ \        / /_   _|  ____|__   __|		*/
/*		   | |  | \  / | (___  \ \  /\  / /  | | | |__     | |   		*/
/*		   | |  | |\/| |\___ \  \ \/  \/ /   | | |  __|    | |   		*/
/*		   | |  | |  | |____) |  \  /\  /   _| |_| |       | |   		*/
/*		   |_|  |_|  |_|_____/    \/  \/   |_____|_|       |_|   		*/
/*		 True Muonium Solver with Iterative Front-Form Techniques		*/
/*											*/
/*--------------------------------------------------------------------------------------*/

#include <slepceps.h>
#include <gsl/gsl_integration.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>

#include "params.cpp"
#include "set_params.cpp"
#include "check_file.cpp"
//#include "sandia_rules_pos.h"
//#include "sandia_rules_pos.c"
#include "sandia_rules.hpp"
#include "sandia_rules.cpp"
#include "prototypes.hpp"
#include "input.cpp"
#include "set_flags.cpp"
#include "print_params.cpp"
#include "string.cpp"
#include "discretize.cpp"
#include "coulomb_cont.cpp"
#include "coulomb_discrete.cpp"
#include "physics.cpp"
#include "h_scat.cpp"
#include "h_anni.cpp"
#include "hamiltonian.cpp"
#include "solver.cpp"
#include "output.cpp"
#include "cleanup.cpp"
static char help[] = "Eigenvalue solver for the problem of QFT bound states on the light front.\n\n"
  "The command line options are not sorted out yet:\n";



int main(int argc, char *argv[])
{
 	PetscErrorCode ierr;
	params params;

	Mat H; 

	SlepcInitialize(&argc,&argv,(char*)0,help);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------------------\n");CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"               _______ __  __  _______          _______ ______ _______                \n");CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"              |__   __|  \\/  |/ ____\\ \\        / /_   _|  ____|__   __|               \n");CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"                 | |  | \\  / | (___  \\ \\  /\\  / /  | | | |__     | |                  \n");CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"                 | |  | |\\/| |\\___ \\  \\ \\/  \\/ /   | | |  __|    | |                  \n");CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"                 | |  | |  | |____) |  \\  /\\  /   _| |_| |       | |                  \n");CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"                 |_|  |_|  |_|_____/    \\/  \\/   |_____|_|       |_|                  \n");CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"               True Muonium Solver with Iterative Front-Form Techniques               \n");CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"                                                                                      \n");CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------------------------------\n");CHKERRQ(ierr);


	read_input(ierr,argv,&params);
	print_input(ierr,&params);
	
	if(check_file(params.hfile))
	{
		PetscViewer    viewer_H;
		PetscViewerBinaryOpen(PETSC_COMM_WORLD,params.hfile.c_str(),FILE_MODE_READ,&viewer_H);
		MatCreate(PETSC_COMM_WORLD,&H);
		MatSetFromOptions(H);
		MatLoad(H,viewer_H);
		PetscViewerDestroy(&viewer_H);		

	}else
	{
        	discretize(ierr,&params);
  	    	ierr = VecView(params.mu,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
//        	ierr = VecView(params.theta,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

		write_output(ierr,&params);	

        	coulomb_trick(ierr,&params);	
//		ierr = VecView(params.CT,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

		ct_discrete(ierr,&params);
//		ierr = VecView(params.CT,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

		hamiltonian(ierr,&params,H);
// 		ierr = PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_DENSE);//DENSE<-->COMMON
//      	MatView(H,PETSC_VIEWER_STDOUT_WORLD);
	}

	cleanup(ierr,&params);	
	eigensolver(ierr,&params,H,argc,argv);

	ierr = MatDestroy(&H);CHKERRQ(ierr);

	ierr = SlepcFinalize();
	return 0;
}

