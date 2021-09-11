#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define nx 101
#define ny 101
#define L 2.0
#define W 1.0
#define pi M_PI
#define alpha 1.0
#define Ti 1
double Told[nx+2][ny+2],Tnew[nx+2][ny+2],K1[nx+2][ny+2],K2[nx+2][ny+2],K3[nx+2][ny+2];
double dt,dx,dy,Tana[nx][ny],t,temp;
double solver();
double initializer();
double output_files();
double analytic();
double output_files0_5();
double output_files0_1();
double output_files0_2();
double output_files1_0();
int main()
{
	initializer();
	solver();
	output_files();
	return 0;
}
double initializer()
{
	int i,j;
	for(i=0;i<=nx+1;i++)
	{
		for(j=0;j<=ny+1;j++)
		{
			Told[i][j]=0.0;
			Tnew[i][j]=0.0;
			K1[i][j]=0.0;
			K3[i][j]=0.0;
			K2[i][j]=0.0;
		}
	}
	for(i=1;i<=nx;i++)
	{
		for(j=1;j<=ny;j++)
		{
			Told[i][j]=1.0;
			Tnew[i][j]=1.0;
			if(i==nx)
			{
				Told[i][j]=0.0;
				Tnew[i][j]=0.0;
			}
			else if(j==ny)
			{
				Told[i][j]=0.0;
				Tnew[i][j]=0.0;	
			}
			
		}
	}	
}
double solver()
{
	int i,j,k,m,n,iteration;iteration=0;
	FILE *steady;
	steady=fopen("Timeintegration_steady_variation.dat","w");
	fprintf(steady,"VARIABLES=\"Time\",\"coded_Temp\",\"analytical_Temp\"\n");
	double error;
	error=1.0;double x,y,errorn;dt=0.00001;
	dx=L/(double)(nx-1);
	dy=W/(double)(ny-1);
	t=0.0;
	while(error>1.0e-6)
	{
		iteration++;
		for(j=1;j<=ny;j++)
		{
			for(i=1;i<=nx;i++)
			{
				K1[i][j]=Told[i][j]+
				        ((dt*alpha/2)*((Told[i-1][j]-2*Told[i][j]+Told[i+1][j])/(dx*dx)+
				        (Told[i][j-1]-2*Told[i][j]+Told[i][j+1])/(dy*dy)));
			}
		}
		for(j=1;j<=ny;j++)
		{
			for(i=1;i<=nx;i++)
			{
				if(i==1)
				{
					K1[i][j]=K1[i+1][j];
				}
				else if(j==1)
				{
					K1[i][j]=K1[i][j+1];
				}
				else if(i==nx)
				{
					K1[i][j]=0.0;
				}
				else if(j==ny)
				{
					K1[i][j]=0.0;
				}
			}
		}				
		for(j=1;j<=ny;j++)
		{
			for(i=1;i<=nx;i++)
			{
				K2[i][j]=Told[i][j]+
				        ((dt*alpha/2)*((K1[i-1][j]-2*K1[i][j]+K1[i+1][j])/(dx*dx)+
				        (K1[i][j-1]-2*K1[i][j]+K1[i][j+1])/(dy*dy)));
			}
		}
		for(j=1;j<=ny;j++)
		{
			for(i=1;i<=nx;i++)
			{
				if(i==1)
				{
					K2[i][j]=K2[i+1][j];
				}
				else if(j==1)
				{
					K2[i][j]=K2[i][j+1];
				}
				else if(i==nx)
				{
					K2[i][j]=0.0;
				}
				else if(j==ny)
				{
					K2[i][j]=0.0;
				}
			}
		}
		for(j=1;j<=ny;j++)
		{
			for(i=1;i<=nx;i++)
			{
				K3[i][j]=Told[i][j]+
				        ((dt*alpha)*((K2[i-1][j]-2*K2[i][j]+K2[i+1][j])/(dx*dx)+
				        (K2[i][j-1]-2*K2[i][j]+K2[i][j+1])/(dy*dy)));
			}
		}
		for(j=1;j<=ny;j++)
		{
			for(i=1;i<=nx;i++)
			{
				if(i==1)
				{
					K3[i][j]=K3[i+1][j];
				}
				else if(j==1)
				{
					K3[i][j]=K3[i][j+1];
				}
				else if(i==nx)
				{
					K3[i][j]=0.0;
				}
				else if(j==ny)
				{
					K3[i][j]=0.0;
				}
			}
		}
		for(j=1;j<=ny;j++)
		{
			for(i=1;i<=nx;i++)
			{
				Tnew[i][j]=Told[i][j]+
				        ((dt*alpha/6)*((K3[i-1][j]-2*K3[i][j]+K3[i+1][j])/(dx*dx)+
				        (K3[i][j-1]-2*K3[i][j]+K3[i][j+1])/(dy*dy)+
						((Told[i-1][j]-2*Told[i][j]+Told[i+1][j])/(dx*dx)+
				        (Told[i][j-1]-2*Told[i][j]+Told[i][j+1])/(dy*dy))+
						2*((K2[i-1][j]-2*K2[i][j]+K2[i+1][j])/(dx*dx)+
				        (K2[i][j-1]-2*K2[i][j]+K2[i][j+1])/(dy*dy))+
						2*((K1[i-1][j]-2*K1[i][j]+K1[i+1][j])/(dx*dx)+
				        (K1[i][j-1]-2*K1[i][j]+K1[i][j+1])/(dy*dy))));
			}
		}
		for(j=1;j<=ny;j++)
		{
			for(i=1;i<=nx;i++)
			{
				if(i==1)
				{
					Tnew[i][j]=Tnew[i+1][j];
				}
				else if(j==1)
				{
					Tnew[i][j]=Tnew[i][j+1];
				}
				else if(i==nx)
				{
					Tnew[i][j]=0.0;
				}
				else if(j==ny)
				{
					Tnew[i][j]=0.0;
				}
			}
		}
		errorn=0.0;
		for(j=1;j<=ny;j++)
		{
			for(i=1;i<=nx;i++)
			{
				errorn+=fabs(Tnew[i][j]-Told[i][j]);
			}
		}error=errorn;
		for(j=1;j<=ny;j++)
		{
			for(i=1;i<=nx;i++)
			{
				Told[i][j]=Tnew[i][j];
			}
		}
		t=dt*iteration;
		temp=Tnew[(nx-1)/4][(ny-1)/4];	
		if(iteration%1000==1)
		{
			analytic();
			printf("time=%lf\titeration=%d\terror=%lf\t\n",t,iteration,error,temp);
			fprintf(steady,"%5.8f\t%5.8f\t%5.8f\n",t,temp,Tana[(nx-1)/4][(ny-1)/4]);
		}
		if(t==0.1)
		{
			analytic();
			output_files0_1();
		}
		if(t==0.2)
		{
			analytic();
			output_files0_2();
		}
		if(t==0.5)
		{
			analytic();
			output_files0_5();
		}
		if(t==1.0)
		{
			analytic();
			output_files1_0();
		}
	}
	printf("t=%lf",t);
	analytic();
	fclose(steady);
}
double output_files0_1()
{
	int i,j;double x,y;
	dx=L/(double)(nx-1);
	dy=W/(double)(ny-1);
	FILE *output,*output1;
	output=fopen("Timeintegration@0_1.dat","w");
	output1=fopen("Timeintegration_x_0.5@0_1.dat","w");
	fprintf(output,"VARIABLES=\"X\",\"Y\",\"coded_Temp\",\"analytical_Temp\"\n");
	fprintf(output,"ZONE F=POINT\n");
	fprintf( output, "I=%d, J=%d\t\n", nx, ny );
	fprintf(output1,"VARIABLES=\"Y\",\"coded_Temp\",\"analytical_Temp\"\n");
	for(j=1;j<=ny;j++)
	{
		for(i=1;i<=nx;i++)
		{
			x=(i-1)*dx;
			y=(j-1)*dy;
			fprintf(output,"%5.8f\t%5.8f\t%5.8f\t%5.8f\n",x,y,Tnew[i][j],Tana[i][j]);	
		}
	}
	for(j=1;j<=ny;j++)
	{
		y=(j-1)*dy;
		fprintf(output1,"%5.8f\t%5.8f\t%5.8f\n",y,Tnew[(nx-1)/2][j],Tana[(nx-1)/2][j]);
	}
	fclose(output);
	fclose(output1);
}
double output_files1_0()
{
	int i,j;double x,y;
	dx=L/(double)(nx-1);
	dy=W/(double)(ny-1);
	FILE *output;
	output=fopen("Timeintegration@1.dat","w");
	fprintf(output,"VARIABLES=\"X\",\"Y\",\"coded_Temp\",\"analytical_Temp\"\n");
	fprintf(output,"ZONE F=POINT\n");
	fprintf( output, "I=%d, J=%d\t\n", nx, ny );
	
	for(j=1;j<=ny;j++)
	{
		for(i=1;i<=nx;i++)
		{
			x=(i-1)*dx;
			y=(j-1)*dy;
			fprintf(output,"%5.8f\t%5.8f\t%5.8f\t%5.8f\n",x,y,Tnew[i][j],Tana[i][j]);	
		}
	}
	fclose(output);
}
double output_files0_2()
{
	int i,j;double x,y;
	dx=L/(double)(nx-1);
	dy=W/(double)(ny-1);
	FILE *output;
	output=fopen("Timeintegration@0_2.dat","w");
	fprintf(output,"VARIABLES=\"X\",\"Y\",\"coded_Temp\",\"analytical_Temp\"\n");
	fprintf(output,"ZONE F=POINT\n");
	fprintf( output, "I=%d, J=%d\t\n", nx, ny );
	for(j=1;j<=ny;j++)
	{
		for(i=1;i<=nx;i++)
		{
			x=(i-1)*dx;
			y=(j-1)*dy;
			fprintf(output,"%5.8f\t%5.8f\t%5.8f\t%5.8f\n",x,y,Tnew[i][j],Tana[i][j]);
		}
	}
	fclose(output);
}
double output_files0_5()
{
	int i,j;double x,y;
	dx=L/(double)(nx-1);
	dy=W/(double)(ny-1);
	FILE *output,*output1,*output2,*output3;
	output=fopen("Timeintegration@0_5.dat","w");
	fprintf(output,"VARIABLES=\"X\",\"coded_Temp\",\"analytical_Temp\"\n");
		for(i=1;i<=nx;i++)
	{
		x=(i-1)*dx;
		fprintf(output,"%5.8f\t%5.8f\t%5.8f\n",x,Tnew[i][(ny-1)/2],Tana[i][(ny-1)/2]);		
	}
	fclose(output);
}
double analytic()
{
	int i,j,m,n; 
	double lamda_m,lamda_n,sum_m,sum_n,x,y;
	dx=L/(double)(nx-1);
	dy=W/(double)(ny-1);
	for(j=1;j<=ny;j++)
	{
		for(i=1;i<=nx;i++)
		{
			x=(i-1)*dx;
			y=(j-1)*dy;
			sum_m=0.0;sum_n=0.0;
			for(m=0;m<=200;m++)
			{
				lamda_m=(2*m+1)*pi/(2*W);
				sum_m+=(pow(-1,m)/(lamda_m*W))*(exp(-alpha*t*lamda_m*lamda_m))*cos(lamda_m*y);
			}
			for(n=0;n<=200;n++)
			{
				lamda_n=(2*n+1)*pi/(2*L);
				sum_n+=((pow(-1,n)/(lamda_n*L)))*(exp(-alpha*t*lamda_n*lamda_n))*cos(lamda_n*x);
			}
			Tana[i][j]=sum_m*sum_n*4;
		}
	}
}
double output_files()
{
	int i,j;double x,y;
	dx=L/(double)(nx-1);
	dy=W/(double)(ny-1);
	FILE *output;
	output=fopen("Timeintegration_steady.dat","w");
	fprintf(output,"VARIABLES=\"X\",\"Y\",\"coded_Temp\",\"analytical_Temp\"\n");
	fprintf(output,"ZONE F=POINT\n");
	fprintf( output, "I=%d, J=%d\t\n", nx, ny );
	for(j=1;j<=ny;j++)
	{
		for(i=1;i<=nx;i++)
		{
			x=(i-1)*dx;
			y=(j-1)*dy;
			fprintf(output,"%5.8f\t%5.8f\t%5.8f\t%5.8f\n",x,y,Tnew[i][j],Tana[i][j]);
		}
	}
	fclose(output);
}
