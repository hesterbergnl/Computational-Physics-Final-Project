#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "comphys.c"
#include "comphys.h"

double vxdot(double x,double y,double vy);
double vydot(double x,double y,double vx);
void initcond();
double x1,vx1,vy1,xsat,xmimas,omega;
double yy1;
double G = 6.67384e-11;
double msat = 568.3e24;
double mmimas = 3.7493e19;
double a = 189176000.0; //separation of mimas + saturn in meters
double pi = 3.141592653;


int main(){
  double k1,l1,m1,n1,k2,l2,m2,n2;
  double k3,l3,m3,n3,k4,l4,m4,n4;
  double t1,t2,tf,x2,y2,vx2,vy2;
  int i,j,k;
  FILE *out,*rsp;
  double dt = 10.0;
  
  
  /*STARTS HERE*/
  double rn,theta,radius,omega2,v,v2;
  long idum = -1;
  time_t now;
  
  now = time(NULL);
  idum = -1*now;
  rn = ran1(&idum);  //initialize
  /*ENDS HERE*/
  
  tf = 1000.0;    //This value is the time that the function will solve the DE for
  
  xsat = (mmimas/(msat+mmimas))*a;
  xmimas = (msat/(msat+mmimas))*a;
  omega = pow(((mmimas+msat)*G)/(a*a*a),0.5);
  
  if((out = fopen("hesterbergnl_project_data.txt","w")) == NULL){
    printf("\nCannot open output file... \n");
	exit(1);
  }
  
  if((rsp = fopen("gnuplot.rsp","w")) == NULL){
    printf("\nCannot open output file for gnuplot... \n");
	exit(1);
  }
  
  for(i=0;i<1000;i++){
  /*start*/

  rn = ran1(&idum);
  theta = 2*pi*rn;
  rn = ran1(&idum);
  radius = (7.45e7 + 6.15e7*rn);
  x1 = radius* cos(theta) + xsat;
  yy1 = radius* sin(theta);
  
  printf("%lf %lf\n",x1,yy1);
  
    //v = 29400.0/(pow(radius,0.5));
    vx1 = ((omega*radius))*cos(theta);
    vy1 = ((omega*radius))*sin(theta);
  /*end*/
  
    //initcond();
	t1 = 0.0;
	while(t1<tf){
	  fprintf(out,"%lf %lf\n",x1,yy1);
	  if(t1 == tf) break;
	  
	  k1 = vxdot(x1,yy1,vy1)*dt;
	  l1 = vydot(x1,yy1,vx1)*dt;
	  m1 = vx1*dt;
	  n1 = vy1*dt;
	  
	  k2 = vxdot((x1+.5*m1),(yy1+.5*n1),(vy1 + .5*l1))*dt;
	  l2 = vydot((x1+.5*m1),(yy1+.5*n1),(vx1 + .5*k1))*dt;
	  m2 = (vx1 + .5*k1)*dt;
	  n2 = (vy1 + .5*l1)*dt;
	  
	  k3 = vxdot((x1+.5*m2),(yy1+.5*n2),(vy1 + .5*l2))*dt;
	  l3 = vydot((x1+.5*m2),(yy1+.5*n2),(vx1 + .5*k2))*dt;
	  m3 = (vx1 + .5*k2)*dt;
	  n3 = (vy1 + .5*l2)*dt;
	  
	  k4 = vxdot((x1+m3),(yy1+n3),(vy1+l3))*dt;
	  l4 = vydot((x1+m3),(yy1+n3),(vx1+k3))*dt;
	  m4 = (vx1 + k3)*dt;
	  n4 = (vy1 + l3)*dt;
	  
	  x2 = x1 + (m1+2.0*m2+2.0*m3+m4)/6.0;
	  y2 = yy1 + (n1+2.0*n2+2.0*n3+n4)/6.0;
	  vx2 = vx1 +(k1+2.0*k2+2.0*k3+k4)/6.0;
	  vy2 = vy1 +(l1+2.0*l2+2.0*l3+l4)/6.0;
	  t2 = t1 + dt;
	  
	  t1 = t2;
	  x1 = x2;
	  yy1 = y2;
	  vx1 = vx2;
	  vy1 = vy2;
	}
  }
  
  fclose(out);
  
  fprintf(rsp,"plot 'hesterbergnl_project_data.txt' using 1:2 with lines\n");
  fprintf(rsp,"pause mouse\n");
  fprintf(rsp,"replot\n");
  fclose(rsp);
  if(system("gnuplot gnuplot.rsp") == -1){
    printf("\nCommand could not be executed\n");
	exit(1);
  }
  return(0);
} 

double vxdot(double x,double y,double vy){
  double temp1,temp2,temp3,temp4;
  
  temp1 = (msat*G*(x-xsat))/pow((((x-xsat)*(x-xsat))+(y*y)),(3.0/2.0));
  temp2 = (mmimas*G*(x-xmimas))/pow((((x-xmimas)*(x-xmimas))+(y*y)),(3.0/2.0));
  temp3 = ((mmimas+msat)*G*x)/(a*a*a);
  temp4 = (2*omega*vy);
  return(-temp1-temp2+temp3+temp4);
}

double vydot(double x,double y,double vx){
  double temp1,temp2,temp3,temp4;
  temp1 = (msat*G*y)/pow((((x-xsat)*(x-xsat))+(y*y)),(3.0/2.0));
  temp2 = (mmimas*G*y)/pow((((x-xmimas)*(x-xmimas))+(y*y)),(3.0/2.0));
  temp3 = ((mmimas+msat)*G*y)/(a*a*a);
  temp4 = (2*omega*vx);
  return(-temp1-temp2+temp3-temp4);
}

/*
void initcond(){
  double rn,theta,radius,omega2,v;
  long i;
  long idum = -1;
  time_t now;
  
  now = time(NULL);
  idum = -1*now;
  rn = ran1(&idum);  //initialize
  
  rn = ran1(&idum);
  theta = 2*pi*rn;
  rn = ran1(&idum);
  radius = (7.45e7 + 6.15e7*rn);
  x1 = radius*cos(theta) + xsat;
  yy1 = radius*sin(theta);
  
  printf("%lf %lf\n",x1,yy1);
  
  v = 29400.0/(pow(radius,0.5));
  //omega2 = v/radius-omega;
  //v2 = omega2*radius;
  vx1 = v*cos(theta);
  vy1 = v*sin(theta);
}
*/
  
  
  
	
	  
	  
