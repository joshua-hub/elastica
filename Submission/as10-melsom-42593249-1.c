// phys3071 as10 melsom 42593249

/*_____________________________________________________________________________
Description: This program solves the shape of a rod bending under compression.
It utilises both the bisection method and the midpoint integration method to 
solve the differential equation

Inputs: The user is requested to enter the applied pressure, the step sizes for
the integration steps. The user must also enter a file name for an output. 

Calculations: This program uses the bisection method to solve for which initial
dtheta/dS solves for a symmetric bend in the rod about s = 0.5. The value theta
must be zero at the point and the values are calculated using the midpoint 
integration method. Wrong initial dtheta/dS values will result in an asymmetric
bend about s = 0.5, which is not a value solution. 

outputs:The created file will have column headers {s (double), theta (double), 
x displacement (double), y displacement (double)}. The name of the file is 
entered by the user. The program also outputs to terminal the maximum value of 
y along the integration.

Compiled as gcc as10-melsom-42593249-1.c -o as10 -Wall -lm
_____________________________________________________________________________*/

// Included Libraries ---------------------------------------------------------
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

// Function Prototypes --------------------------------------------------------
double integrator (double DeltaS, double P, double init_ThetaDash);
double bisection (double DeltaS, double P);
void output (double Pressure, double dTheta, FILE *dat);

// Begin main function --------------------------------------------------------
int main () {
  double dtheta_root; // The value of dTheta such that theta(0.5) = 0
  double Pressure; // User provided pressure value
  double DeltaS; // User provided step sizes
  char fname[100]; // User provided output file name
  FILE *dat; // File pointer
  
  // User input section -------------------------------------------------------
  // User input of Pressure.
  printf("Please enter the load to be placed through the beam: ");
  do { // Making sure only correct values can be provided
    scanf("%lf", &Pressure);
    if (Pressure< 0.0) { // Condition, positive pressure.
      printf("Error, please input a positive number: ");
    }
  } while (Pressure< 0.0);
  
  // User input of step size
  printf("Please enter the step size Delta s: ");
  do { // Making sure only correct values can be provided
    scanf("%lf", &DeltaS);
    if (DeltaS< 0.0) { // Condition, positive step sizes.
      printf("Error, please input a positive number: ");
    }
  } while (DeltaS< 0.0);
  
  // User input of output file name
  printf("Please type the output file name: ");
  scanf("%s", fname);

  // Call functions -----------------------------------------------------------
  dtheta_root = bisection(DeltaS, Pressure);
  
  // Open the file before calling output function. Is done this way to not have
  // to pass extra items to output function.
  dat= fopen(fname, "w"); 
  output(Pressure, dtheta_root, dat);
  fclose(dat);
  
  return (EXIT_SUCCESS);
}

// FUNCTIONS ------------------------------------------------------------------

// Bisection method function --------------------------------------------------
// This function performs the bisection method on Theta(s= 0.5) as a function 
// of initial condition dTheta. The upper and lower values of dTheta converge
// to the root, which is returned.
double bisection(double DeltaS, double Pressure) {
  double theta_root;     // The root of Theta as a function of dTheta/dS
  double thresh = 1e-8;  // Not given, but went to sufficiently small value
  double dThetaA= 0.5;   // As required in assignment sheet
  double dThetaB= 30.0;  // As required in assignment sheet
  double dThetaBisecting= (dThetaA+ dThetaB)/ 2.0; // Bisecting value
  double lower;          // Theta evaluated with initial condition dThetaA
  
  lower= integrator(DeltaS, Pressure, dThetaA); 
  // Loop for bisection method, calling for integrating function to return the 
  // values of Theta(S= 0.5).
  do {
    // If the bisection is in the lower half
    if (lower* (integrator(DeltaS, Pressure, dThetaBisecting)) < 0.0) {
      (dThetaB= dThetaBisecting); // move upper bound
    }
    
    // If the bisection is in the upper half
    else {
      (dThetaA= dThetaBisecting); // move lower bound
      lower= integrator(DeltaS, Pressure, dThetaBisecting); //recalculate lower
    }
  
  dThetaBisecting= (dThetaA+ dThetaB)/ 2.0; // Recalculate bisecting value
  } while (fabs(dThetaA- dThetaB) > thresh); // Bisecting until convergence
  
  theta_root=(dThetaA+ dThetaB)/ 2.0;
  return (theta_root);
}

// Integrator function --------------------------------------------------------
// This function ingrates the differential equation from s= 0to s= 0.5 given
// the applied pressure, the step sizes and the initial condition of dTheta/dS.
double integrator(double DeltaS, double Pressure, double dTheta) {
  int steps= 1;
  double t_old= dTheta, t_half, t_new;
  double theta_old= 0.0, theta_half, theta_new; // initial angle is zero
  double s_position= 0.0; // starting position
  
  do {
    // mid point method of integration to solve the coupled equations.
    theta_half= theta_old + (DeltaS/ 2.0)* t_old;
    t_half= t_old- (DeltaS/ 2.0)* Pressure* sin(theta_old);    
    theta_new= theta_old+ DeltaS* t_half;
    t_new= t_old- DeltaS* Pressure* sin(theta_half);
    
    // Reinitialise the 'old values' for the next step
    theta_old = theta_new;
    t_old = t_new;
    
    steps++; // Integer step counter to avoid round off from summing doubles
    s_position= (1.0* steps)* DeltaS;
  } while (s_position<= 0.5);
  
  return (theta_new);
}

// outputting function --------------------------------------------------------
// This function is used only to output x and y position values. The root of 
// the equation has already been solved and is used as the initial value in 
// this loop. This output uses fixed step sizes of 1e-3 and outputs to file 
// every 10 steps to produce 100 evenly spaced values of S and the 
// corresponding Theta, x and y.
void output (double Pressure, double dTheta, FILE *dat) {
  double DeltaS= 1e-3, s_position= 0.0, y_max= 0.0;
  double x_old= 0.0, x_new; // x positions
  double y_old= 0.0, y_new; // y positions
  int steps= 1;
  double t_old= dTheta, t_half, t_new;
  double theta_old= 0.0, theta_half, theta_new;
  
  do {
    // mid point method of integration to solve the coupled equations.
    theta_half= theta_old + (DeltaS/ 2.0)* t_old;
    t_half= t_old- (DeltaS/ 2.0)* Pressure* sin(theta_old);    
    theta_new= theta_old+ DeltaS* t_half;
    t_new= t_old- DeltaS* Pressure* sin(theta_half);
    
    // Also solve for x and y displacement. No need for x_half or y_half, 
    // because dx/dS or dy/dS and not functions of x or y respectively.
    x_new= x_old+ DeltaS* cos(theta_half); 
    y_new= y_old+ DeltaS* sin(theta_half);
    
    if (((steps- 1)% 10)== 0) // Used to print 10 values to file 
     (fprintf(dat, "%lf\t%lf\t%lf\t%lf\n",s_position, theta_old, x_new, y_new));
    if (y_new> y_max) (y_max= y_new); // Finds the maximum y deflection
    
    // Reinitialise the 'old values' for the next step
    x_old= x_new;
    y_old= y_new;
    theta_old= theta_new;
    t_old= t_new;
    
    steps++; // Integer step counter to avoid round off from summing doubles
    s_position= (1.0* steps)* DeltaS;
  } while (s_position<= 1.0); 
  printf("The height of maximum deflection (y_max) is %lf\n\n", y_max);
}