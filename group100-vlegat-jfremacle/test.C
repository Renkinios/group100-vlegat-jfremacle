// femProblem *femElasticityCreate(const char *filename,femGeo* theGeometry) 
// {
//     femProblem *theProblem = malloc(sizeof(femProblem));
//     femElasticCase iCase; 
//     FILE* file = fopen(filename,"r");
//     if (file == NULL) Error("No mesh file !");
//     ErrorScan(fscanf(file,"Type of problem : %c \n", &iCase));
//     ErrorScan(fscanf(file,"Young modulus : %d \n", &theProblem->E));
//     ErrorScan(fscanf(file,"Poisson ratio : %d \n", &theProblem->nu));
//     ErrorScan(fscanf(file,"Mass density : %d \n", &theProblem->rho));
//     ErrorScan(fscanf(file,"Gravity : %d \n", &theProblem->g));
//     double E = theProblem->E;
//     double nu = theProblem->nu;
//     double rho = theProblem->rho;
//     double g = theProblem->g;