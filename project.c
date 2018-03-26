#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "comphys.c"
#include "comphys.h"

double vxdot(double x,double y,double vy);
double vydot(double x,double y,double vx);
void initcond();
double x1,y1,vx1,vy1,xsat,xmimas,omega;
double dt = 1.0;
double G = 6.67384e-11;
double msat = 568.3e24;
double mmimas = 3.7493e19;
double a = 189176000; 

xsat = mmimas/(msat+mmimas);
xmimas = msat/(msat+mmimas);
omega = pow(((mmimas+msat)*G)/(a*a*a),0.5);

int main(){
  double k1,l1,m1,n1,k2,l2,m2,n2;
  double k3,l3,m3,n3,k3,l4,m4,n4;
  double t1,t2,tf,x2,y2,vx2,vy2;
  int i,j,k;
  FILE *out,*rsp;
  
  tf = 10000.0;    //This value is the time that the function will solve the DE for
  
  if((out = fopen("hesterbergnl_project_data.dat","w")) == NULL){
    printf("\nCannot open output file... \n");
	exit(1);
  }
  
  if((rsp = fopen("gnuplot.rsp","w")) == NULL){
    printf("\nCannot open output file for gnuplot... \n");
	exit(1);
  }
  
  for(i=0;i<tf;i++){
    initcond();
	while(t<tf){
	  fprintf(out,"%lf %lf\n",t1,x1,y1);
	  if(t1 == tf) break;
	  
	  k1 = vxdot(x1,y1,vy1);
	  l1 = vydot(x1,y1,vx1);
	  m1 = vx;
	  n1 = vy;
	  
	  k2 = vxdot((x1+.5*m1),(y1+.5*n1),(vy1 + .5*l1));
	  l2 = vydot((x1+.5*m1),(y1+.5*n1),(vx1 + .5*k1));
	  m2 = vx + .5*k1;
	  n2 = vy + .5*l1;
	  
	  k3 = vxdot((x1+.5*m2),(y1+.5*n2),(vy1 + .5*l2));
	  l3 = vydot((x1+.5*m2),(y1+.5*n2),(vx1 + .5*k2));
	  m3 = vx + .5*k2;
	  n3 = vy + .5*l2;
	  
	  k4 = vxdot((x1+m3),(y1+n3),(vy1+l3));
	  l4 = vydot((x1+m3),(y1+n3),(vx1+k3));
	  m4 = vx + k3;
	  n4 = vy + l3;
	  
	  x2 = x1 + (m1+2.0*m2+2.0*m3+m4);
	  y2 = y1 + (n1+2.0*n2+2.0*n3+n4);
	  vx2 = vx +(k1+2.0*k2+2.0*k3+k4);
	  vy2 = vy +(l1+2.0*l2+2.0*l3+l4);
	  t2 = t1 + dt;
	  
	  t1 = t2;
	  x1 = x2;
	  y1 = y2;
	  vx1 = vx2;
	  vy1 = vy2;
	}
  }
  
  fclose(out);
  
  fprintf(rsp,"plot 'hesterbergnl_project_data.dat' using 1:2 with lines\n");
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
  return(temp1-temp2+temp3+temp4);
}

double vydot(double x,double y,double vx){
  double temp1,temp2,temp3,temp4;
  temp1 = (msat*G*y)/pow((((x-xsat)*(x-xsat))+(y*y)),(3.0/2.0));
  temp2 = (mmimas*G*y)/pow((((x-xmimas)*(x-xmimas))+(y*y)),(3.0/2.0));
  temp3 = ((mmimas+msat)*G*y)/(a*a*a);
  temp4 = (2*omega*vx);
  return(temp1-temp2+temp3-temp4);
}

void initcond(){
  double rn;
  long i;
  long idum = -1;
  time_t now;
  
  now = time(NULL);
  idum = -1*now;
  rn = ran1(&idum);  //initialize
  rn = ran1(&idum);
  
  
  x1 = 0.0;
  y1 = 60000000.0 + rn*125000000.0;   //picking random starting value between surface of saturn and moon mimas
  vy1 = 0.0;    //starts at rest in y direction so it is moving tangential to saturn;
  rn = ran1(&idum);
  vx1 = -30000.0 + rn*60000.0;
}

  
  
  
	
	  
	  