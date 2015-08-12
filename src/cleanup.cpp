void cleanup(PetscErrorCode ierr, params *params)
{

        ierr = VecDestroy(&params->mu);CHKERRV(ierr);
        ierr = VecDestroy(&params->wmu);CHKERRV(ierr);
        ierr = VecDestroy(&params->theta);CHKERRV(ierr);
        ierr = VecDestroy(&params->wtheta);CHKERRV(ierr);
        ierr = VecDestroy(&params->CT);CHKERRV(ierr);

}
