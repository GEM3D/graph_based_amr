#include "typedefs.h"

void Struct_3D(Cube& cube,int id,int N, int M,int L,double **X, double **Y,double **Z);
void VTK_solution(double *X, double*Y, double*Z,int N, int M,int L,int ID,double *q);
void save(Cube& cube,Cube& cube_temp, int taged_cube);
void add_kids(Cube& cube,int id);
void AMR_Refine_Cube(Cube& cube,int id, int L, int M, int N);
void Amr_four_to_one_enforce(Cube& cube,Vector_Int& refine_list,int L,int M, int N);
void AMR_Derefine_Cube(Cube& cube,int parent_id);
void write_hdf5_serial(Cube& cube, int L, int M, int N);
void sphere(Cube& cube, Vector_Int& refine_list,double C);
void elem_conn_check(Cube& cube);
void circle(double **x, double **y,int n);
void identify(Cube& cube, Vector_Int& refine_list,double *x,double *y,int n);
void check_coordinates(Cube&cube);

int main(int argcs, char* pArgs[])
{

Cube cube;
cube.push_back(cube_data());
// reserve your initial guess of the vector
cube.reserve(50000);

int nface=6;
int nchild=8;
unsigned int i,j,k,l;


//cube[0].ID=0;
cube[0].xyz[0]=-1.0;
cube[0].xyz[1]=1.0;
cube[0].xyz[2]=-1.0;
cube[0].xyz[3]=1.0;
cube[0].xyz[4]=-1.0;
cube[0].xyz[5]=1.0;


for(i=0;i<6;i++)
{
	cube[0].nbr[i].push_back(-i-1);
	printf("%d\n",cube[0].nbr[i].at(0));
}


double *X=NULL;
double *Y=NULL;
double *Z=NULL;

int L=10;
int M=L;
int N=L;

// allocate for the original 

double *q=NULL;
q=(double*)calloc(N*M*L,sizeof(double));

int q_size=L*M*N;

//cube[0].q=(double*)malloc(q_size*sizeof(double));


//VTK_solution(X,Y,Z,N,M,L,0,q);

// divide each elem by 8
/*
Cube cube_temp;
cube_temp.push_back(cube_data());
*/
//save(cube,cube_temp,0);

AMR_Refine_Cube(cube,0,L,M,N);

double *x=NULL;
double *y=NULL;
int n=200;

circle(&x,&y,n);
/*
for(i=0;i<n;i++)
{
	printf("%lf %lf\n",x[i],y[i]);
}
*/

//AMR_Refine_Cube(cube,0,L,M,N);
//AMR_Refine_Cube(cube,8,L,M,N);
Vector_Int refine_list;
/*
refine_list.push_back(8);
Amr_four_to_one_enforce(cube,refine_list,L,M,N);
refine_list.push_back(8);
Amr_four_to_one_enforce(cube,refine_list,L,M,N);
*/
double c=1.;
int l1,l2;

for(j=0;j<0;j++)
{
//sphere(cube,refine_list,c);
identify(cube,refine_list,x,y,n);
/*
for(i=0;i<refine_list.size();i++)
	{
		printf("refine_list=%d\n",refine_list.at(i));
	 }
*/
Amr_four_to_one_enforce(cube,refine_list,L,M,N);

for(l1=0;l1<cube.size();l1++)
{
	//printf("============================ %d \n ",l1);
	//printf("Adaptation Level%d \n ",cube[l1].level);
	
	for(l2=0;l2<nface;l2++)
	{
		if(cube[l1].nbr[l2].size()>4)
		//for(k=0;k<cube[l1].nbr[l2].size();k++)
		{
			//printf("%d\t",cube[l1].nbr[l2].at(k));
			printf("element at fault%d\n",l1);
		}
		//printf("\n");
	}
	
}

#if(0)
 check_coordinates(cube);
 elem_conn_check(cube);
#endif

	
	// c=c*.5;
	 printf("%d\n",cube.size());
}


/*

refine_list.push_back(0);
refine_list.push_back(8);

Amr_four_to_one_enforce(cube,refine_list,L,M,N);
*/

// for debugging the following routine, I mess up the connectivity to see if I can catch it
//cube[0].nbr[5].at(0)=2; 

#if(0)
{
elem_conn_check(cube);

}
#endif

//printf("right here\n");

/*
for(l1=0;l1<cube.size();l1++)
{
	printf("============================ %d \n ",l1);
	printf("Adaptation Level%d \n ",cube[l1].level);
	
	for(l2=0;l2<nface;l2++)
	{
		for(k=0;k<cube[l1].nbr[l2].size();k++)
		{
			printf("%d\t",cube[l1].nbr[l2].at(k));
		}
		printf("\n");
	}
	
}
*/
//AMR_Derefine_Cube(cube,1);

//remove_kids_compact(cube,0);



/*
save(cube,cube_temp,1);

add_kids(cube,cube_temp,1);


save(cube,cube_temp,2);

add_kids(cube,cube_temp,2);
*/

for(i=0;i<0;i++)
{
//save(cube,cube_temp,i);

//add_kids(cube,cube_temp,i);
	
AMR_Refine_Cube(cube,i,L,M,N);	
}

/*
for(l=0;l<cube.size();l++)
{
	printf("============================ %d \n ",l);
	for(j=0;j<nface;j++)
	{
		for(k=0;k<cube[l].nbr[j].size();k++)
		{
			printf("%d\t",cube[l].nbr[j].at(k));
		}
		printf("\n");
	}
	
}
*/
//

int id;
/*
for(id=0;id<cube.size();id++)
{
Struct_3D(cube,id,N,M,L,&X,&Y,&Z); 

for(i=0;i<N*M*L;i++)
{
	q[i]=(double)id;
}
VTK_solution(X,Y,Z,N,M,L,id,q);
}
*/

write_hdf5_serial(cube,L,M,N);
#if(0)
for(i=0;i<cube.size();i++)
{
	printf("============================ %d \n ",i);
	for(j=0;j<nface;j++)
	{
		for(k=0;k<cube[i].nbr[j].size();k++)
		{
			printf("%d\t",cube[i].nbr[j].at(k));
		}
		printf("\n");
	}
	
}
#endif

/*
for(i=0;i<cube.size();i++)
{
	printf("============================ %d \n ",i);
	for(j=0;j<nface;j++)
	{
		printf("%lf\t",cube[i].xyz[j]);	
	}
	printf("\n");
	
}
*/



// clean up 
/*
for(i=0;i<cube.size();i++)
{
	free(cube[i].q);
}
*/
Cube cube_dealloc;

cube.swap(cube_dealloc);


free(q);
free(X);
free(Y);
free(Z);


}
//

void Struct_3D(Cube& cube,int id,int N, int M,int L,double **X, double **Y,double **Z)
{
	int i;
	
	double Xa=cube[id].xyz[0];
	double Xb=cube[id].xyz[1];
	double Ya=cube[id].xyz[2];
	double Yb=cube[id].xyz[3];
	double Za=cube[id].xyz[4];
	double Zb=cube[id].xyz[5];
		
	double hx=L-1.0;
	double hy=M-1.0;
	double hz=N-1.0;
	double Xh=(Xb-Xa)/(hx);
	double Yh=(Yb-Ya)/(hy);
	double Zh=(Zb-Za)/(hz);
	
	(*X)=(double*)malloc(L*sizeof(double));
	(*Y)=(double*)malloc(M*sizeof(double));	
	(*Z)=(double*)malloc(N*sizeof(double));
	
	for(i=0;i<L;i++)
	{
	(*X)[i]=Xa+Xh*i;		
	//printf("%lf %lf\n",X[i],Y[i]);	
	}

	for(i=0;i<M;i++)
	{
	(*Y)[i]=Ya+Yh*i;	
	}
		
	for(i=0;i<N;i++)
	{
	(*Z)[i]=Za+Zh*i;	
	}
		
}


//======================================================================

void VTK_solution(double *X, double*Y, double*Z,int N, int M,int L,int ID,double *q)
{
int i,j,k;	
	
	FILE  *fp=NULL;	
	
	// here we get some data into variable data
    char filename[64];
    sprintf (filename, "out%d.vtk",ID);
    fp = fopen(filename, "w");
   
	
   // fp=fopen("out.vtk","w");
    fprintf(fp,"# vtk DataFile Version 2.0 \n");
    fprintf(fp,"Grid\n");
    fprintf(fp,"ASCII\n");
    fprintf(fp,"DATASET STRUCTURED_GRID\n");
    fprintf(fp,"DIMENSIONS %d %d %d\n",N,M,L); 
    fprintf(fp,"POINTS %d float\n",M*N*L);
    
    for(k=0;k<N;k++)
    {
    for(j=0;j<M;j++)
    {
    for(i=0;i<L;i++)	
	{		
	fprintf(fp,"%lf %lf %lf\n",X[i],Y[j],Z[k]);		
	}
    }
    }
   fprintf(fp,"POINT_DATA %d\n",M*N*L);
    fprintf(fp,"SCALARS q float 1\n");
    fprintf(fp,"LOOKUP_TABLE default\n");
	for(j=0;j<M*N*L;j++)
	{		
	fprintf(fp,"%lf\n",q[j]);		
	}
	fclose(fp);
	
	
}

//=======================================================================

void save(Cube& cube,Cube& cube_temp, int tagged_cube)
{
	unsigned int i,j;
	
	for(i=0;i<6;i++)
	{
	cube_temp[0].xyz[i]=cube[tagged_cube].xyz[i];
    }
	
	for(i=0;i<6;i++)
   {
	   for(j=0;j<cube[tagged_cube].nbr[i].size();j++)
	   {
	cube_temp[0].nbr[i].push_back(-1);
       }
	//printf("%d\n",cube[0].nbr[i].at(0));
   }
	
  for(i=0;i<6;i++)
	{
		for(j=0;j<cube[tagged_cube].nbr[i].size();j++)
		{
	     cube_temp[0].nbr[i].at(j)=cube[tagged_cube].nbr[i].at(j);
       }
    }
}
//======================================================================

void add_kids(Cube& cube,int id)
{
 unsigned int i,j;
	
 
 Cube cube_temp;
 
 cube_temp.push_back(cube_data());
 
 for(i=0;i<6;i++)
	{
	cube_temp[0].xyz[i]=cube[id].xyz[i];
    }
	
	for(i=0;i<6;i++)
   {
	   for(j=0;j<cube[id].nbr[i].size();j++)
	   {
	cube_temp[0].nbr[i].push_back(-1);
       }
	//printf("%d\n",cube[0].nbr[i].at(0));
   }
	
  for(i=0;i<6;i++)
	{
		for(j=0;j<cube[id].nbr[i].size();j++)
		{
	     cube_temp[0].nbr[i].at(j)=cube[id].nbr[i].at(j);
       }
    }
 
 
 int nquad=cube.size();
 int index=nquad;	
	
 for(i=0;i<7;i++)
 {
  cube.push_back(cube_data());
 }
 
 // eliminate the old neighbor data from the parent element since now parent will be ith kid
 
 for(i=0;i<6;i++)
 {
 cube[id].nbr[i].clear();
}

double xmid=(cube_temp[0].xyz[0]+cube_temp[0].xyz[1])*0.5;
double ymid=(cube_temp[0].xyz[2]+cube_temp[0].xyz[3])*0.5;
double zmid=(cube_temp[0].xyz[4]+cube_temp[0].xyz[5])*0.5;

printf("%lf %lf %lf\n",xmid,ymid,zmid);
// inside the refined element update 
//======================================================================

//                     kid 0 is located at 0,0,0 

//======================================================================

  cube[id].xyz[0]=cube_temp[0].xyz[0];  	
  cube[id].xyz[1]=xmid;
  
  cube[id].xyz[2]=cube_temp[0].xyz[2];
  cube[id].xyz[3]=ymid;
  
  cube[id].xyz[4]=cube_temp[0].xyz[4];
  cube[id].xyz[5]=zmid;
  
  // update the kid neighbor, neigbors z=0, y=0 and x=0 dont change
  // update the other 3 neighbors
  // sit in the middle of the element, 
  // 0 >> bottom 1>> top 
  // 2 >> right  3>>left
  // 4 >> back   5>> front
  // some of the neghbirs will be inherited from the parent element
  
  unsigned int elem_id;
  double x0,x1,y0,y1,z0,z1;
  double xc,yc,zc;
  double xc_id,yc_id,zc_id; 
   
  cube[id].nbr[1].push_back(nquad+3);
  cube[id].nbr[3].push_back(nquad+2);  
  cube[id].nbr[5].push_back(nquad);
  
  Vector_Dbl center;
  double r; 
  int idx;
   
 if(cube_temp[0].nbr[0].size()==1)
 {
   cube[id].nbr[0].push_back(cube_temp[0].nbr[0].at(0));  
   
    if(cube_temp[0].nbr[0].at(0)>0)
  {
   cube[cube_temp[0].nbr[0].at(0)].nbr[1].push_back(nquad);
   cube[cube_temp[0].nbr[0].at(0)].nbr[1].push_back(nquad+1);
   cube[cube_temp[0].nbr[0].at(0)].nbr[1].push_back(nquad+2);
  }
 }
 else
 {
	 // call splitting function
	     xc_id= 0.5*(cube[id].xyz[0]+cube[id].xyz[1]);		 
		 yc_id= 0.5*(cube[id].xyz[2]+cube[id].xyz[3]);  
		 
		 r=1000.0;
		
	 for(i=0;i<cube_temp[0].nbr[0].size();i++)
	 {
		 x0= cube[cube_temp[0].nbr[0].at(i)].xyz[0];
		 x1= cube[cube_temp[0].nbr[0].at(i)].xyz[1]; 
		 
		 y0= cube[cube_temp[0].nbr[0].at(i)].xyz[2];
		 y1= cube[cube_temp[0].nbr[0].at(i)].xyz[3]; 
		 
		 xc=0.5*(x0+x1);
		 yc=0.5*(y0+y1);
		 
		 if((xc-xc_id)*(xc-xc_id)+(yc-yc_id)*(yc-yc_id)<r)
		 {
			idx=i;
			r=(xc-xc_id)*(xc-xc_id)+(yc-yc_id)*(yc-yc_id); 
		 }		 		 
	 }
	 //printf("idx =%d\n",idx);
	 cube[id].nbr[0].push_back(cube_temp[0].nbr[0].at(idx));
 }
 
 // due to my numbering convention the neighboring element for id'th element is not going to change
 /*
 for(i=0;i<cube[cube_temp[0].nbr[0].at(idx)].nbr[1].size();i++)
 {
	 if(cube[cube_temp[0].nbr[0].at(idx)].nbr[1].at(i)==id)
	 {
         		 
	 }
 }
 */
 
 if(cube_temp[0].nbr[2].size()==1)
 {
 cube[id].nbr[2].push_back(cube_temp[0].nbr[2].at(0));
 
  if(cube_temp[0].nbr[2].at(0)>0)
  {
   cube[cube_temp[0].nbr[2].at(0)].nbr[3].push_back(nquad+3);
   cube[cube_temp[0].nbr[2].at(0)].nbr[3].push_back(nquad+4);
   cube[cube_temp[0].nbr[2].at(0)].nbr[3].push_back(nquad+1);
  }
 }
 else
 {
	 // call splitting function
	 
	     xc_id= cube[id].xyz[0]+cube[id].xyz[1];		 
		 zc_id= cube[id].xyz[4]+cube[id].xyz[5];  
		 
		 r=1000.0;
		
	 for(i=0;i<cube_temp[0].nbr[2].size();i++)
	 {		 
		 x0= cube[cube_temp[0].nbr[2].at(i)].xyz[0];
		 x1= cube[cube_temp[0].nbr[2].at(i)].xyz[1]; 
		 
		 z0= cube[cube_temp[0].nbr[2].at(i)].xyz[4];
		 z1= cube[cube_temp[0].nbr[2].at(i)].xyz[5]; 
		 
		 xc=0.5*(x0+x1);
		 zc=0.5*(z0+z1);
		 
		 if((xc-xc_id)*(xc-xc_id)+(zc-zc_id)*(zc-zc_id)<r)
		 {
			idx=i;
			r=(xc-xc_id)*(xc-xc_id)+(zc-zc_id)*(zc-zc_id); 
		 }		 		 
	 }
	 
	 cube[id].nbr[2].push_back(cube_temp[0].nbr[2].at(idx));
 }
 

 if(cube_temp[0].nbr[4].size()==1)
 {
  cube[id].nbr[4].push_back(cube_temp[0].nbr[4].at(0));
 
 //printf("%d\n",cube_temp[0].nbr[4].at(0));
 
  if(cube_temp[0].nbr[4].at(0)>0)
  {  
  cube[cube_temp[0].nbr[4].at(0)].nbr[5].push_back(nquad+2);
  cube[cube_temp[0].nbr[4].at(0)].nbr[5].push_back(nquad+6);
  cube[cube_temp[0].nbr[4].at(0)].nbr[5].push_back(nquad+3);
  }
 }
 else
 {
	 // call splitting function	 
	     yc_id= cube[id].xyz[2]+cube[id].xyz[3];
		 zc_id= cube[id].xyz[4]+cube[id].xyz[5];  
		 
		 r=1000.0;
		
	 for(i=0;i<cube_temp[0].nbr[4].size();i++)
	 {
		 
		 y0= cube[cube_temp[0].nbr[4].at(i)].xyz[2];
		 y1= cube[cube_temp[0].nbr[4].at(i)].xyz[3]; 
		 
		 z0= cube[cube_temp[0].nbr[4].at(i)].xyz[4];
		 z1= cube[cube_temp[0].nbr[4].at(i)].xyz[5]; 
		 
		 yc=0.5*(y0+y1);
		 zc=0.5*(z0+z1);
		 
		 if((yc-yc_id)*(yc-yc_id)+(zc-zc_id)*(zc-zc_id)<r)
		 {
			idx=i;
			r=(yc-yc_id)*(yc-yc_id)+(zc-zc_id)*(zc-zc_id); 
		 }		 		 
	 }

	 cube[id].nbr[4].push_back(cube_temp[0].nbr[4].at(idx));
 }
 
 
 
 
 //=====================================================================
   
  //           kid 1 which is nquad'th element located at 1,0,0 
 
 //=====================================================================
  int nbr_id;
 
  cube[index].xyz[0]=xmid;  	
  cube[index].xyz[1]=cube_temp[0].xyz[1];
  
  cube[index].xyz[2]=cube_temp[0].xyz[2];
  cube[index].xyz[3]=ymid;
  
  cube[index].xyz[4]=cube_temp[0].xyz[4];
  cube[index].xyz[5]=zmid;
      
  
  cube[index].nbr[1].push_back(nquad+4); 
  
  cube[index].nbr[3].push_back(nquad+1); 
  
  cube[index].nbr[4].push_back(id); 
  
  if(cube_temp[0].nbr[0].size()==1)
  { 
	cube[index].nbr[0].push_back(cube_temp[0].nbr[0].at(0));   
  }
  else
  {
	     xc_id= cube[index].xyz[0]+cube[index].xyz[1];		 
		 yc_id= cube[index].xyz[2]+cube[index].xyz[3];  
		 
		 r=1000.0;
		
	 for(i=0;i<cube_temp[0].nbr[0].size();i++)
	 {
		 x0= cube[cube_temp[0].nbr[0].at(i)].xyz[0];
		 x1= cube[cube_temp[0].nbr[0].at(i)].xyz[1]; 
		 
		 y0= cube[cube_temp[0].nbr[0].at(i)].xyz[2];
		 y1= cube[cube_temp[0].nbr[0].at(i)].xyz[3]; 
		 
		 xc=0.5*(x0+x1);
		 yc=0.5*(y0+y1);
		 
		 if((xc-xc_id)*(xc-xc_id)+(yc-yc_id)*(yc-yc_id)<r)
		 {
			idx=i;
			r=(xc-xc_id)*(xc-xc_id)+(yc-yc_id)*(yc-yc_id); 
		 }		 		 
	 }
	 
	 nbr_id=cube_temp[0].nbr[0].at(idx);
	 cube[index].nbr[0].push_back(nbr_id);
  
	  for(i=0;i<cube[nbr_id].nbr[1].size();i++)
	 {
		 if(cube[nbr_id].nbr[1].at(i)==id)
		 {
			cube[nbr_id].nbr[1].at(i)=index;
		 }
	 }  
  }
  // update the neighbors nbr vector
   
  if(cube_temp[0].nbr[2].size()==1)
  {
	cube[index].nbr[2].push_back(cube_temp[0].nbr[2].at(0)); 
  }
  else
  {	  
	     xc_id= cube[index].xyz[0]+cube[index].xyz[1];
		 
		 zc_id= cube[index].xyz[4]+cube[index].xyz[5];  
		 
		 r=1000.0;
		
	 for(i=0;i<cube_temp[0].nbr[2].size();i++)
	 {
		 
		 x0= cube[cube_temp[0].nbr[2].at(i)].xyz[0];
		 x1= cube[cube_temp[0].nbr[2].at(i)].xyz[1]; 
		 
		 z0= cube[cube_temp[0].nbr[2].at(i)].xyz[4];
		 z1= cube[cube_temp[0].nbr[2].at(i)].xyz[5]; 
		 
		 xc=0.5*(x0+x1);
		 zc=0.5*(z0+z1);
		 
		 if((xc-xc_id)*(xc-xc_id)+(zc-zc_id)*(zc-zc_id)<r)
		 {
			idx=i;
			r=(xc-xc_id)*(xc-xc_id)+(zc-zc_id)*(zc-zc_id); 
		 }		 		 
	 }
	 
	  nbr_id=cube_temp[0].nbr[2].at(idx);
	  cube[index].nbr[2].push_back(nbr_id); 
	  
	   for(i=0;i<cube[nbr_id].nbr[3].size();i++)
	 {
		 if(cube[nbr_id].nbr[3].at(i)==id)
		 {
			cube[nbr_id].nbr[3].at(i)=index;
		 }
	 }  
	  
	  
  }
  
   if(cube_temp[0].nbr[5].size()==1)
  {
	cube[index].nbr[5].push_back(cube_temp[0].nbr[5].at(0)); 
	
	if(cube_temp[0].nbr[5].at(0)>0)
	{
	cube[cube_temp[0].nbr[5].at(0)].nbr[4].pop_back();	
	
   cube[cube_temp[0].nbr[5].at(0)].nbr[4].push_back(nquad);
   cube[cube_temp[0].nbr[5].at(0)].nbr[4].push_back(nquad+1);
   cube[cube_temp[0].nbr[5].at(0)].nbr[4].push_back(nquad+5);
   cube[cube_temp[0].nbr[5].at(0)].nbr[4].push_back(nquad+4);
  
   }
  }
  else
  {
	     yc_id= cube[index].xyz[2]+cube[index].xyz[3];
		 zc_id= cube[index].xyz[4]+cube[index].xyz[5];  
		 
		 r=1000.0;
		
	 for(i=0;i<cube_temp[0].nbr[5].size();i++)
	 {
		 
		 y0= cube[cube_temp[0].nbr[5].at(i)].xyz[2];
		 y1= cube[cube_temp[0].nbr[5].at(i)].xyz[3]; 
		 
		 z0= cube[cube_temp[0].nbr[5].at(i)].xyz[4];
		 z1= cube[cube_temp[0].nbr[5].at(i)].xyz[5]; 
		 
		 yc=0.5*(y0+y1);
		 zc=0.5*(z0+z1);
		 
		 if((yc-yc_id)*(yc-yc_id)+(zc-zc_id)*(zc-zc_id)<r)
		 {
			idx=i;
			r=(yc-yc_id)*(yc-yc_id)+(zc-zc_id)*(zc-zc_id); 
		 }		 		 
	 }
	 
	 nbr_id=cube_temp[0].nbr[5].at(idx);
	 cube[index].nbr[5].push_back(nbr_id);
	 
	  for(i=0;i<cube[nbr_id].nbr[4].size();i++)
	 {
		 if(cube[nbr_id].nbr[4].at(i)==id)
		 {
			cube[nbr_id].nbr[4].at(i)=index;
		 }
	 }    
  }
  
 //=====================================================================    

  //          kid 2 is (nquad+1)'th element located at 1,1,0 
  
  //====================================================================
 
  index++;
 
  cube[index].xyz[0]=xmid;  	
  cube[index].xyz[1]=cube_temp[0].xyz[1];
  
  cube[index].xyz[2]=ymid;
  cube[index].xyz[3]=cube_temp[0].xyz[3];
  
  cube[index].xyz[4]=cube_temp[0].xyz[4];
  cube[index].xyz[5]=zmid;
  
  
  cube[index].nbr[1].push_back(nquad+5); 
  cube[index].nbr[2].push_back(nquad); 
  cube[index].nbr[4].push_back(nquad+2); 
   
  
  if(cube_temp[0].nbr[0].size()==1)
  {
	cube[index].nbr[0].push_back(cube_temp[0].nbr[0].at(0)); 
  }
  else
  {
	     xc_id= cube[index].xyz[0]+cube[index].xyz[1];
		 
		 yc_id= cube[index].xyz[2]+cube[index].xyz[3];  
		 
		 r=1000.0;
		
	 for(i=0;i<cube_temp[0].nbr[0].size();i++)
	 {
		 
		 x0= cube[cube_temp[0].nbr[0].at(i)].xyz[0];
		 x1= cube[cube_temp[0].nbr[0].at(i)].xyz[1]; 
		 
		 y0= cube[cube_temp[0].nbr[0].at(i)].xyz[2];
		 y1= cube[cube_temp[0].nbr[0].at(i)].xyz[3]; 
		 
		 xc=0.5*(x0+x1);
		 yc=0.5*(y0+y1);
		 
		 if((xc-xc_id)*(xc-xc_id)+(yc-yc_id)*(yc-yc_id)<r)
		 {
			idx=i;
			r=(xc-xc_id)*(xc-xc_id)+(yc-yc_id)*(yc-yc_id); 
		 }		 		 
	 }
	 
	 nbr_id=cube_temp[0].nbr[0].at(idx);
	 cube[index].nbr[0].push_back(nbr_id);
	 
	  for(i=0;i<cube[nbr_id].nbr[1].size();i++)
	 {
		 if(cube[nbr_id].nbr[1].at(i)==id)
		 {
			cube[nbr_id].nbr[1].at(i)=index;
		 }
	 }    
	 
	 
  }
  
  if(cube_temp[0].nbr[3].size()==1)
  {
	cube[index].nbr[3].push_back(cube_temp[0].nbr[3].at(0));
	
	if(cube_temp[0].nbr[3].at(0)>0)
	{
   cube[cube_temp[0].nbr[3].at(0)].nbr[2].pop_back();	
   cube[cube_temp[0].nbr[3].at(0)].nbr[2].push_back(nquad+2);
   cube[cube_temp[0].nbr[3].at(0)].nbr[2].push_back(nquad+1);
   cube[cube_temp[0].nbr[3].at(0)].nbr[2].push_back(nquad+5);
   cube[cube_temp[0].nbr[3].at(0)].nbr[2].push_back(nquad+6);
   }
  }
  else
  {
	    xc_id= cube[index].xyz[0]+cube[index].xyz[1];
		 
		 zc_id= cube[index].xyz[4]+cube[index].xyz[5];  
		 
		 r=1000.0;
		
	 for(i=0;i<cube_temp[0].nbr[3].size();i++)
	 {
		 
		 x0= cube[cube_temp[0].nbr[3].at(i)].xyz[0];
		 x1= cube[cube_temp[0].nbr[3].at(i)].xyz[1]; 
		 
		 z0= cube[cube_temp[0].nbr[3].at(i)].xyz[4];
		 z1= cube[cube_temp[0].nbr[3].at(i)].xyz[5]; 
		 
		 xc=0.5*(x0+x1);
		 zc=0.5*(z0+z1);
		 
		 if((xc-xc_id)*(xc-xc_id)+(zc-zc_id)*(zc-zc_id)<r)
		 {
			idx=i;
			r=(xc-xc_id)*(xc-xc_id)+(zc-zc_id)*(zc-zc_id); 
		 }		 		 
	 }
	 
	   nbr_id=cube_temp[0].nbr[3].at(idx); 
	   cube[index].nbr[0].push_back(nbr_id);
	 
	  for(i=0;i<cube[nbr_id].nbr[2].size();i++)
	 {
		 if(cube[nbr_id].nbr[2].at(i)==id)
		 {
			cube[nbr_id].nbr[2].at(i)=index;
		 }
	 }    
	  
	  
  }
  
  if(cube_temp[0].nbr[5].size()==1)
  {
	cube[index].nbr[5].push_back(cube_temp[0].nbr[5].at(0)); 
  }
  else
  {
	  
	  yc_id= cube[index].xyz[2]+cube[index].xyz[3];
		 zc_id= cube[index].xyz[4]+cube[index].xyz[5];  
		 
		 r=1000.0;
		
	 for(i=0;i<cube_temp[0].nbr[5].size();i++)
	 {
		 
		 y0= cube[cube_temp[0].nbr[5].at(i)].xyz[0];
		 y1= cube[cube_temp[0].nbr[5].at(i)].xyz[1]; 
		 
		 z0= cube[cube_temp[0].nbr[5].at(i)].xyz[4];
		 z1= cube[cube_temp[0].nbr[5].at(i)].xyz[5]; 
		 
		 yc=0.5*(y0+y1);
		 zc=0.5*(z0+z1);
		 
		 if((yc-yc_id)*(yc-yc_id)+(zc-zc_id)*(zc-zc_id)<r)
		 {
			idx=i;
			r=(yc-yc_id)*(yc-yc_id)+(zc-zc_id)*(zc-zc_id); 
		 }		 		 
	 }
	 
	 nbr_id=cube_temp[0].nbr[5].at(idx);
	 cube[index].nbr[5].push_back(nbr_id);
	 
	  for(i=0;i<cube[nbr_id].nbr[4].size();i++)
	 {
		 if(cube[nbr_id].nbr[4].at(i)==id)
		 {
			cube[nbr_id].nbr[4].at(i)=index;
		 }
	 }     	  
  }
  
 //=====================================================================
 
 //        kid 3 is (nquad+2)th element located at 1,1,0 

//======================================================================

  index++;
 
  cube[index].xyz[0]=cube_temp[0].xyz[0];  	
  cube[index].xyz[1]=xmid;
  
  cube[index].xyz[2]=ymid;
  cube[index].xyz[3]=cube_temp[0].xyz[3];
  
  cube[index].xyz[4]=cube_temp[0].xyz[4];
  cube[index].xyz[5]=zmid;
  
  cube[index].nbr[1].push_back(nquad+6); 
  cube[index].nbr[2].push_back(id); 
  cube[index].nbr[5].push_back(nquad+1); 
   
  if(cube_temp[0].nbr[0].size()==1)
  {
	cube[index].nbr[0].push_back(cube_temp[0].nbr[0].at(0)); 
  }
  else
  {
	     xc_id= cube[index].xyz[0]+cube[index].xyz[1];
		 
		 yc_id= cube[index].xyz[2]+cube[index].xyz[3];  
		 
		 r=1000.0;
		
	 for(i=0;i<cube_temp[0].nbr[0].size();i++)
	 {
		 
		 x0= cube[cube_temp[0].nbr[0].at(i)].xyz[0];
		 x1= cube[cube_temp[0].nbr[0].at(i)].xyz[1]; 
		 
		 y0= cube[cube_temp[0].nbr[0].at(i)].xyz[2];
		 y1= cube[cube_temp[0].nbr[0].at(i)].xyz[3]; 
		 
		 xc=0.5*(x0+x1);
		 yc=0.5*(y0+y1);
		 
		 if((xc-xc_id)*(xc-xc_id)+(yc-yc_id)*(yc-yc_id)<r)
		 {
			idx=i;
			r=(xc-xc_id)*(xc-xc_id)+(yc-yc_id)*(yc-yc_id); 
		 }		 		 
	 }
	 nbr_id=cube_temp[0].nbr[0].at(idx);
	 cube[index].nbr[0].push_back(nbr_id);
	 
	  for(i=0;i<cube[nbr_id].nbr[1].size();i++)
	 {
		 if(cube[nbr_id].nbr[1].at(i)==id)
		 {
			cube[nbr_id].nbr[1].at(i)=index;
		 }
	 }     	    
  }
  
  if(cube_temp[0].nbr[3].size()==1)
  {
	cube[index].nbr[3].push_back(cube_temp[0].nbr[3].at(0)); 
  }
  else
  {
	     xc_id= cube[index].xyz[0]+cube[index].xyz[1];		 
		 zc_id= cube[index].xyz[4]+cube[index].xyz[5];  
		 
		 r=1000.0;
		
	 for(i=0;i<cube_temp[0].nbr[3].size();i++)
	 {
		 
		 x0= cube[cube_temp[0].nbr[3].at(i)].xyz[0];
		 x1= cube[cube_temp[0].nbr[3].at(i)].xyz[1]; 
		 
		 z0= cube[cube_temp[0].nbr[3].at(i)].xyz[4];
		 z1= cube[cube_temp[0].nbr[3].at(i)].xyz[5]; 
		 
		 xc=0.5*(x0+x1);
		 zc=0.5*(z0+z1);
		 
		 if((xc-xc_id)*(xc-xc_id)+(zc-zc_id)*(zc-zc_id)<r)
		 {
			idx=i;
			r=(xc-xc_id)*(xc-xc_id)+(zc-zc_id)*(zc-zc_id); 
		 }		 		 
	 }
	 
	  nbr_id=cube_temp[0].nbr[3].at(idx);
	  cube[index].nbr[3].push_back(nbr_id);
	  
	   for(i=0;i<cube[nbr_id].nbr[2].size();i++)
	 {
		 if(cube[nbr_id].nbr[2].at(i)==id)
		 {
			cube[nbr_id].nbr[2].at(i)=index;
		 }
	 }     	     	  
  }
  
  if(cube_temp[0].nbr[4].size()==1)
  {
	cube[index].nbr[4].push_back(cube_temp[0].nbr[4].at(0)); 
  }
  else
  {
	   yc_id= cube[index].xyz[2]+cube[index].xyz[3];
		 zc_id= cube[index].xyz[4]+cube[index].xyz[5];  
		 
		 r=1000.0;
		
	 for(i=0;i<cube_temp[0].nbr[4].size();i++)
	 {
		 
		 y0= cube[cube_temp[0].nbr[4].at(i)].xyz[0];
		 y1= cube[cube_temp[0].nbr[4].at(i)].xyz[1]; 
		 
		 z0= cube[cube_temp[0].nbr[4].at(i)].xyz[4];
		 z1= cube[cube_temp[0].nbr[4].at(i)].xyz[5]; 
		 
		 yc=0.4*(y0+y1);
		 zc=0.4*(z0+z1);
		 
		 if((yc-yc_id)*(yc-yc_id)+(zc-zc_id)*(zc-zc_id)<r)
		 {
			idx=i;
			r=(yc-yc_id)*(yc-yc_id)+(zc-zc_id)*(zc-zc_id); 
		 }		 		 
	 }	 
	 
	  nbr_id=cube_temp[0].nbr[4].at(idx);
	  cube[index].nbr[4].push_back(nbr_id);
	  
	   for(i=0;i<cube[nbr_id].nbr[5].size();i++)
	 {
		 if(cube[nbr_id].nbr[5].at(i)==id)
		 {
			cube[nbr_id].nbr[5].at(i)=index;
		 }
	 }     	  	  
  }
   
   
  //====================================================================
   
  //          kid 4 is (nquad+3)th element  located at 1,1,0 
 
 //=====================================================================
 
  index++;
 
  cube[index].xyz[0]=cube_temp[0].xyz[0];  	
  cube[index].xyz[1]=xmid;
  
  cube[index].xyz[2]=cube_temp[0].xyz[2];
  cube[index].xyz[3]=ymid;
  
  cube[index].xyz[4]=zmid;
  cube[index].xyz[5]=cube_temp[0].xyz[5];
  
   
   
  cube[index].nbr[0].push_back(id); 
  cube[index].nbr[3].push_back(nquad+6); 
  cube[index].nbr[5].push_back(nquad+4); 
   
  if(cube_temp[0].nbr[1].size()==1)
  {
	cube[index].nbr[1].push_back(cube_temp[0].nbr[1].at(0)); 
	
	if(cube_temp[0].nbr[1].at(0)>0)
	{
   cube[cube_temp[0].nbr[1].at(0)].nbr[0].pop_back();	
   cube[cube_temp[0].nbr[1].at(0)].nbr[0].push_back(nquad+3);
   cube[cube_temp[0].nbr[1].at(0)].nbr[0].push_back(nquad+4);
   cube[cube_temp[0].nbr[1].at(0)].nbr[0].push_back(nquad+5);
   cube[cube_temp[0].nbr[1].at(0)].nbr[0].push_back(nquad+6);
}
	
  }
  else
  {
	     xc_id= cube[index].xyz[0]+cube[index].xyz[1];		 
		 yc_id= cube[index].xyz[2]+cube[index].xyz[3];  
		 
		 r=1000.0;
		
	 for(i=0;i<cube_temp[0].nbr[1].size();i++)
	 {
		 
		 x0= cube[cube_temp[0].nbr[1].at(i)].xyz[0];
		 x1= cube[cube_temp[0].nbr[1].at(i)].xyz[1]; 
		 
		 y0= cube[cube_temp[0].nbr[1].at(i)].xyz[2];
		 y1= cube[cube_temp[0].nbr[1].at(i)].xyz[3]; 
		 
		 xc=0.5*(x0+x1);
		 yc=0.5*(z0+z1);
		 
		 if((xc-xc_id)*(xc-xc_id)+(yc-yc_id)*(yc-yc_id)<r)
		 {
			idx=i;
			r=(xc-xc_id)*(xc-xc_id)+(yc-yc_id)*(yc-yc_id); 
		 }		 		 
	 }
	  
	   nbr_id=cube_temp[0].nbr[1].at(idx);
	  cube[index].nbr[1].push_back(nbr_id);
	  
	   for(i=0;i<cube[nbr_id].nbr[0].size();i++)
	 {
		 if(cube[nbr_id].nbr[0].at(i)==id)
		 {
			cube[nbr_id].nbr[0].at(i)=index;
		 }
	 }     	  	 
	  
  }
  
  if(cube_temp[0].nbr[2].size()==1)
  {
	cube[index].nbr[2].push_back(cube_temp[0].nbr[2].at(0)); 
  }
  else
  {
	  
	     xc_id= cube[index].xyz[0]+cube[index].xyz[1];
		 
		 zc_id= cube[index].xyz[4]+cube[index].xyz[5];  
		 
		 r=1000.0;
		
	 for(i=0;i<cube_temp[0].nbr[2].size();i++)
	 {
		 
		 x0= cube[cube_temp[0].nbr[2].at(i)].xyz[0];
		 x1= cube[cube_temp[0].nbr[2].at(i)].xyz[1]; 
		 
		 z0= cube[cube_temp[0].nbr[2].at(i)].xyz[4];
		 z1= cube[cube_temp[0].nbr[2].at(i)].xyz[5]; 
		 
		 xc=0.5*(x0+x1);
		 zc=0.5*(z0+z1);
		 
		 if((xc-xc_id)*(xc-xc_id)+(zc-zc_id)*(zc-zc_id)<r)
		 {
			idx=i;
			r=(xc-xc_id)*(xc-xc_id)+(zc-zc_id)*(zc-zc_id); 
		 }		 		 
	 }
	   
	  nbr_id=cube_temp[0].nbr[2].at(idx);
	  cube[index].nbr[2].push_back(nbr_id);
	  
	   for(i=0;i<cube[nbr_id].nbr[3].size();i++)
	 {
		 if(cube[nbr_id].nbr[3].at(i)==id)
		 {
			cube[nbr_id].nbr[3].at(i)=index;
		 }
	 }     	  		  
  }
  
  
  if(cube_temp[0].nbr[4].size()==1)
  {
	cube[index].nbr[4].push_back(cube_temp[0].nbr[4].at(0)); 
  }
  else
  {
	     yc_id= cube[index].xyz[2]+cube[index].xyz[3];
		 zc_id= cube[index].xyz[4]+cube[index].xyz[5];  
		 
		 r=1000.0;
		
	 for(i=0;i<cube_temp[0].nbr[4].size();i++)
	 {
		 
		 y0= cube[cube_temp[0].nbr[4].at(i)].xyz[0];
		 y1= cube[cube_temp[0].nbr[4].at(i)].xyz[1]; 
		 
		 z0= cube[cube_temp[0].nbr[4].at(i)].xyz[4];
		 z1= cube[cube_temp[0].nbr[4].at(i)].xyz[5]; 
		 
		 yc=0.4*(y0+y1);
		 zc=0.4*(z0+z1);
		 
		 if((yc-yc_id)*(yc-yc_id)+(zc-zc_id)*(zc-zc_id)<r)
		 {
			idx=i;
			r=(yc-yc_id)*(yc-yc_id)+(zc-zc_id)*(zc-zc_id); 
		 }		 		 
	 }
	  
	  nbr_id=cube_temp[0].nbr[4].at(idx);
	  cube[index].nbr[4].push_back(nbr_id);
	  
	   for(i=0;i<cube[nbr_id].nbr[5].size();i++)
	 {
		 if(cube[nbr_id].nbr[5].at(i)==id)
		 {
			cube[nbr_id].nbr[5].at(i)=index;
		 }
	 }     	  		
	 
	 
  }
   
   
//======================================================================
    
//        kid 5 is (nquad+4)th element  located at 1,1,0 
 
//======================================================================


  index++;
 
  cube[index].xyz[0]=xmid;  	
  cube[index].xyz[1]=cube_temp[0].xyz[1];
  
  cube[index].xyz[2]=cube_temp[0].xyz[2];
  cube[index].xyz[3]=ymid;
  
  cube[index].xyz[4]=zmid;
  cube[index].xyz[5]=cube_temp[0].xyz[5];
  
  
  cube[index].nbr[0].push_back(nquad); 
  cube[index].nbr[3].push_back(nquad+5); 
  cube[index].nbr[4].push_back(nquad+3); 
   
  if(cube_temp[0].nbr[1].size()==1)
  {
	cube[index].nbr[1].push_back(cube_temp[0].nbr[1].at(0)); 
  }
  else
  {
	     xc_id= cube[index].xyz[0]+cube[index].xyz[1];
		 
		 yc_id= cube[index].xyz[2]+cube[index].xyz[3];  
		 
		 r=1000.0;
	  
	   for(i=0;i<cube_temp[0].nbr[1].size();i++)
	 {
		 
		 x0= cube[cube_temp[0].nbr[1].at(i)].xyz[0];
		 x1= cube[cube_temp[0].nbr[1].at(i)].xyz[1]; 
		 
		 y0= cube[cube_temp[0].nbr[1].at(i)].xyz[2];
		 y1= cube[cube_temp[0].nbr[1].at(i)].xyz[3]; 
		 
		 xc=0.5*(x0+x1);
		 yc=0.5*(z0+z1);
		 
		 if((xc-xc_id)*(xc-xc_id)+(yc-yc_id)*(yc-yc_id)<r)
		 {
			idx=i;
			r=(xc-xc_id)*(xc-xc_id)+(yc-yc_id)*(yc-yc_id); 
		 }		 		 
	 }
	 	  
	  nbr_id=cube_temp[0].nbr[1].at(idx);
	  cube[index].nbr[1].push_back(nbr_id);
	  
	   for(i=0;i<cube[nbr_id].nbr[0].size();i++)
	 {
		 if(cube[nbr_id].nbr[0].at(i)==id)
		 {
			cube[nbr_id].nbr[0].at(i)=index;
		 }
	 }     	  	
	  
	  
  }
  
  if(cube_temp[0].nbr[2].size()==1)
  {
	cube[index].nbr[2].push_back(cube_temp[0].nbr[2].at(0)); 
  }
  else
  {
	     xc_id= cube[index].xyz[0]+cube[index].xyz[1];
		 
		 zc_id= cube[index].xyz[4]+cube[index].xyz[5];  
		 
		 r=1000.0;
		
	 for(i=0;i<cube_temp[0].nbr[2].size();i++)
	 {
		 
		 x0= cube[cube_temp[0].nbr[2].at(i)].xyz[0];
		 x1= cube[cube_temp[0].nbr[2].at(i)].xyz[1]; 
		 
		 z0= cube[cube_temp[0].nbr[2].at(i)].xyz[4];
		 z1= cube[cube_temp[0].nbr[2].at(i)].xyz[5]; 
		 
		 xc=0.5*(x0+x1);
		 zc=0.5*(z0+z1);
		 
		 if((xc-xc_id)*(xc-xc_id)+(zc-zc_id)*(zc-zc_id)<r)
		 {
			idx=i;
			r=(xc-xc_id)*(xc-xc_id)+(zc-zc_id)*(zc-zc_id); 
		 }		 		 
	 }
	   
	   nbr_id=cube_temp[0].nbr[2].at(idx);
	  cube[index].nbr[2].push_back(nbr_id);
	  
	   for(i=0;i<cube[nbr_id].nbr[3].size();i++)
	 {
		 if(cube[nbr_id].nbr[3].at(i)==id)
		 {
			cube[nbr_id].nbr[3].at(i)=index;
		 }
	 }     	  	
	  
  }
  
  if(cube_temp[0].nbr[5].size()==1)
  {
	cube[index].nbr[5].push_back(cube_temp[0].nbr[5].at(0)); 
  }
  else
  {
	    yc_id= cube[index].xyz[2]+cube[index].xyz[3];
		 zc_id= cube[index].xyz[5]+cube[index].xyz[5];  
		 
		 r=1000.0;
		
	 for(i=0;i<cube_temp[0].nbr[5].size();i++)
	 {
		 
		 y0= cube[cube_temp[0].nbr[5].at(i)].xyz[0];
		 y1= cube[cube_temp[0].nbr[5].at(i)].xyz[1]; 
		 
		 z0= cube[cube_temp[0].nbr[5].at(i)].xyz[5];
		 z1= cube[cube_temp[0].nbr[5].at(i)].xyz[5]; 
		 
		 yc=0.5*(y0+y1);
		 zc=0.5*(z0+z1);
		 
		 if((yc-yc_id)*(yc-yc_id)+(zc-zc_id)*(zc-zc_id)<r)
		 {
			idx=i;
			r=(yc-yc_id)*(yc-yc_id)+(zc-zc_id)*(zc-zc_id); 
		 }		 		 
	 }
	 
	 
	  nbr_id=cube_temp[0].nbr[5].at(idx);
	  cube[index].nbr[5].push_back(nbr_id);
	  
	   for(i=0;i<cube[nbr_id].nbr[4].size();i++)
	 {
		 if(cube[nbr_id].nbr[4].at(i)==id)
		 {
			cube[nbr_id].nbr[4].at(i)=index;
		 }
	 }     	  	
	 
	 
  }
  
 //=====================================================================
  
             //kid 6 is (nquad+5)th element located at 1,1,0 
 
 //=====================================================================
 
  index++;
 
  cube[index].xyz[0]=xmid;  	
  cube[index].xyz[1]=cube_temp[0].xyz[1];
  
  cube[index].xyz[2]=ymid;
  cube[index].xyz[3]=cube_temp[0].xyz[3];
  
  cube[index].xyz[4]=zmid;
  cube[index].xyz[5]=cube_temp[0].xyz[5];
  
  cube[index].nbr[0].push_back(nquad+1); 
  cube[index].nbr[2].push_back(nquad+4); 
  cube[index].nbr[4].push_back(nquad+6); 
   
  if(cube_temp[0].nbr[1].size()==1)
  {
	cube[index].nbr[1].push_back(cube_temp[0].nbr[1].at(0)); 
  }
  else
  {
	     xc_id= cube[index].xyz[1]+cube[index].xyz[2];
		 yc_id= cube[index].xyz[2]+cube[index].xyz[3];  
		 
		 r=1000.0;
		 
	    for(i=0;i<cube_temp[0].nbr[1].size();i++)
	 {
		 
		 x0= cube[cube_temp[0].nbr[1].at(i)].xyz[0];
		 x1= cube[cube_temp[0].nbr[1].at(i)].xyz[1]; 
		 
		 y0= cube[cube_temp[0].nbr[1].at(i)].xyz[2];
		 y1= cube[cube_temp[0].nbr[1].at(i)].xyz[3]; 
		 
		 xc=0.5*(x0+x1);
		 yc=0.5*(z0+z1);
		 
		 if((xc-xc_id)*(xc-xc_id)+(yc-yc_id)*(yc-yc_id)<r)
		 {
			idx=i;
			r=(xc-xc_id)*(xc-xc_id)+(yc-yc_id)*(yc-yc_id); 
		 }		 		 
	 }
	  
	  nbr_id=cube_temp[0].nbr[1].at(idx);
	  cube[index].nbr[1].push_back(nbr_id);
	  
	   for(i=0;i<cube[nbr_id].nbr[0].size();i++)
	 {
		 if(cube[nbr_id].nbr[0].at(i)==id)
		 {
			cube[nbr_id].nbr[0].at(i)=index;
		 }
	 }     	  	
	  
  }
  
  if(cube_temp[0].nbr[3].size()==1)
  {
	cube[index].nbr[3].push_back(cube_temp[0].nbr[3].at(0)); 
  }
  else
  {
	     xc_id= cube[index].xyz[0]+cube[index].xyz[1];		 
		 zc_id= cube[index].xyz[4]+cube[index].xyz[5];  
		 
		 r=1000.0;
		
	 for(i=0;i<cube_temp[0].nbr[3].size();i++)
	 {
		 
		 x0= cube[cube_temp[0].nbr[3].at(i)].xyz[0];
		 x1= cube[cube_temp[0].nbr[3].at(i)].xyz[1]; 
		 
		 z0= cube[cube_temp[0].nbr[3].at(i)].xyz[4];
		 z1= cube[cube_temp[0].nbr[3].at(i)].xyz[5]; 
		 
		 xc=0.5*(x0+x1);
		 zc=0.5*(z0+z1);
		 
		 if((xc-xc_id)*(xc-xc_id)+(zc-zc_id)*(zc-zc_id)<r)
		 {
			idx=i;
			r=(xc-xc_id)*(xc-xc_id)+(zc-zc_id)*(zc-zc_id); 
		 }		 		 
	 }
	 
	    
	   nbr_id=cube_temp[0].nbr[3].at(idx);
	  cube[index].nbr[3].push_back(nbr_id);
	  
	   for(i=0;i<cube[nbr_id].nbr[2].size();i++)
	 {
		 if(cube[nbr_id].nbr[2].at(i)==id)
		 {
			cube[nbr_id].nbr[2].at(i)=index;
		 }
	 }     	  	
	  
  }
  
  if(cube_temp[0].nbr[5].size()==1)
  {
	cube[index].nbr[5].push_back(cube_temp[0].nbr[5].at(0)); 
  }
  else
  {
	  
	     yc_id= cube[index].xyz[2]+cube[index].xyz[3];
		 zc_id= cube[index].xyz[4]+cube[index].xyz[5];  
		 
		 r=1000.0;
		 
	   for(i=0;i<cube_temp[0].nbr[5].size();i++)
	 {
		 
		 y0= cube[cube_temp[0].nbr[5].at(i)].xyz[0];
		 y1= cube[cube_temp[0].nbr[5].at(i)].xyz[1]; 
		 
		 z0= cube[cube_temp[0].nbr[5].at(i)].xyz[5];
		 z1= cube[cube_temp[0].nbr[5].at(i)].xyz[5]; 
		 
		 yc=0.5*(y0+y1);
		 zc=0.5*(z0+z1);
		 
		 if((yc-yc_id)*(yc-yc_id)+(zc-zc_id)*(zc-zc_id)<r)
		 {
			idx=i;
			r=(yc-yc_id)*(yc-yc_id)+(zc-zc_id)*(zc-zc_id); 
		 }		 		 
	 }	 
	  nbr_id=cube_temp[0].nbr[5].at(idx);
	  cube[index].nbr[5].push_back(nbr_id);
	  
	   for(i=0;i<cube[nbr_id].nbr[4].size();i++)
	 {
		 if(cube[nbr_id].nbr[4].at(i)==id)
		 {
			cube[nbr_id].nbr[4].at(i)=index;
		 }
	 }     	  	 
  }
  
  
//======================================================================

	 //kid 7 (nquad+6)th element is located at 1,1,0 

//======================================================================
  index++;
 
  cube[index].xyz[0]=cube_temp[0].xyz[0];  	
  cube[index].xyz[1]=xmid;
  
  cube[index].xyz[2]=ymid;
  cube[index].xyz[3]=cube_temp[0].xyz[3];
  
  cube[index].xyz[4]=zmid;
  cube[index].xyz[5]=cube_temp[0].xyz[5];
  
  cube[index].nbr[0].push_back(nquad+2); 
  cube[index].nbr[2].push_back(nquad+3); 
  cube[index].nbr[5].push_back(nquad+5); 
   
  if(cube_temp[0].nbr[1].size()==1)
  {
	cube[index].nbr[1].push_back(cube_temp[0].nbr[1].at(0)); 
  }
  else
  {
	  xc_id= cube[index].xyz[1]+cube[index].xyz[2];
		 yc_id= cube[index].xyz[2]+cube[index].xyz[3];  
		 
		 r=1000.0;
		 
	    for(i=0;i<cube_temp[0].nbr[1].size();i++)
	 {
		 
		 x0= cube[cube_temp[0].nbr[1].at(i)].xyz[0];
		 x1= cube[cube_temp[0].nbr[1].at(i)].xyz[1]; 
		 
		 y0= cube[cube_temp[0].nbr[1].at(i)].xyz[2];
		 y1= cube[cube_temp[0].nbr[1].at(i)].xyz[3]; 
		 
		 xc=0.5*(x0+x1);
		 yc=0.5*(z0+z1);
		 
		 if((xc-xc_id)*(xc-xc_id)+(yc-yc_id)*(yc-yc_id)<r)
		 {
			idx=i;
			r=(xc-xc_id)*(xc-xc_id)+(yc-yc_id)*(yc-yc_id); 
		 }		 		 
	 }
	   
	   nbr_id=cube_temp[0].nbr[1].at(idx);
	  cube[index].nbr[1].push_back(nbr_id);
	  
	   for(i=0;i<cube[nbr_id].nbr[0].size();i++)
	 {
		 if(cube[nbr_id].nbr[0].at(i)==id)
		 {
			cube[nbr_id].nbr[0].at(i)=index;
		 }
	 }     	  	
  }
  
  if(cube_temp[0].nbr[3].size()==1)
  {
	cube[index].nbr[3].push_back(cube_temp[0].nbr[3].at(0)); 
  }
  else
  {
	     xc_id= cube[index].xyz[0]+cube[index].xyz[1];		 
		 zc_id= cube[index].xyz[4]+cube[index].xyz[5];  
		 
		 r=1000.0;
		
	 for(i=0;i<cube_temp[0].nbr[3].size();i++)
	 {
		 
		 x0= cube[cube_temp[0].nbr[3].at(i)].xyz[0];
		 x1= cube[cube_temp[0].nbr[3].at(i)].xyz[1]; 
		 
		 z0= cube[cube_temp[0].nbr[3].at(i)].xyz[4];
		 z1= cube[cube_temp[0].nbr[3].at(i)].xyz[5]; 
		 
		 xc=0.5*(x0+x1);
		 zc=0.5*(z0+z1);
		 
		 if((xc-xc_id)*(xc-xc_id)+(zc-zc_id)*(zc-zc_id)<r)
		 {
			idx=i;
			r=(xc-xc_id)*(xc-xc_id)+(zc-zc_id)*(zc-zc_id); 
		 }		 		 
	 }
	 	  
	   nbr_id=cube_temp[0].nbr[3].at(idx);
	  cube[index].nbr[3].push_back(nbr_id);
	  
	   for(i=0;i<cube[nbr_id].nbr[2].size();i++)
	 {
		 if(cube[nbr_id].nbr[2].at(i)==id)
		 {
			cube[nbr_id].nbr[2].at(i)=index;
		 }
	 }     	  	
  }
  
  if(cube_temp[0].nbr[4].size()==1)
  {
	cube[index].nbr[4].push_back(cube_temp[0].nbr[4].at(0)); 
  }
  else
  {
	    yc_id= cube[index].xyz[2]+cube[index].xyz[3];
		 zc_id= cube[index].xyz[4]+cube[index].xyz[5];  
		 
		 r=1000.0;
		
	 for(i=0;i<cube_temp[0].nbr[4].size();i++)
	 {
		 y0= cube[cube_temp[0].nbr[4].at(i)].xyz[0];
		 y1= cube[cube_temp[0].nbr[4].at(i)].xyz[1]; 
		 
		 z0= cube[cube_temp[0].nbr[4].at(i)].xyz[4];
		 z1= cube[cube_temp[0].nbr[4].at(i)].xyz[4]; 
		 
		 yc=0.4*(y0+y1);
		 zc=0.4*(z0+z1);
		 
		 if((yc-yc_id)*(yc-yc_id)+(zc-zc_id)*(zc-zc_id)<r)
		 {
			idx=i;
			r=(yc-yc_id)*(yc-yc_id)+(zc-zc_id)*(zc-zc_id); 
		 }		 		 
	 }
	 
	 
	  nbr_id=cube_temp[0].nbr[4].at(idx);
	  cube[index].nbr[4].push_back(nbr_id);
	  
	   for(i=0;i<cube[nbr_id].nbr[5].size();i++)
	 {
		 if(cube[nbr_id].nbr[5].at(i)==id)
		 {
			cube[nbr_id].nbr[5].at(i)=index;
		 }
	 }     	  	
	 
  }
  
//======================================================================

// Erase cube_temp's nbrs for future use

//======================================================================

for(i=0;i<6;i++)
{
	while(cube_temp[0].nbr[i].size()!=0)
	{
	cube_temp[0].nbr[i].pop_back();
    }
    
    //cube_temp[0].nbr[i].clear();
}
	
}

//======================================================================









