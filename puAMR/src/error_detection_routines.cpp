//#include "/home/jhasbestan/GEM_AMR/src/include/typedefs.h"
#include "typedefs.h"

//======================================================================
//
//  This routine checks to see if elemnt b is a neighbor of element b, 
//        then  element b should be a neighbor of element a
//
//======================================================================

void elem_conn_check(Cube& cube)
{
	 int elem_id,nbr_id,a;
	unsigned int i,j,k,l1,m;
	printf("*******************************************************************\n");
	printf("\n");
	printf("\n"); 
	printf("                   Processor Local Element Connectivity Check ... \n"); 
	printf("\n");
	printf("\n"); 
	
	for(i=0;i<cube.size();i++)
	{
	  for(j=0;j<6;j++)
		{
		 for(k=0;k<cube[i].nbr[j].size();k++)
		 {
			nbr_id=cube[i].nbr[j].at(k);
			  
		    a=0;
		// printf("==============nbr_id= %d\n",nbr_id); 
		 if(nbr_id>-1)
		 {		  
			// printf("==============nbr_id= %d\n",nbr_id); 
		  for(l1=0;l1<6;l1++)
		  {
			for(m=0;m<cube[nbr_id].nbr[l1].size();m++)
			{
				if(cube[nbr_id].nbr[l1].at(m)==(int)i)
				{
					a=1;
					//printf("==============a= %d",a);
			    }
			}
			if(a==1)
			{
				break;
			}    
		  }
		  
		  if(a==0)
		  {
			  printf("Error in connectivity, 'elem_conn_check' function exited\n");
			  printf("\n");
	          printf("*******************************************************************\n"); 
			  exit(0);
		  } 
	  }
		  }
		  
		  
		}
		
	}
	
	
	printf("\n"); 
	printf("               Processor Local Element Connectivity Check Successful \n"); 
	printf("\n");
	printf("*******************************************************************\n"); 
	
	
}
//======================================================================

void CheckConnectivityConsistency(Cube &cube)
{

for(unsigned int i=0;i<cube.size();i++)
{
	for(unsigned int j=0;j<6;j++)
{
	if(cube[i].nbr[j].size()!=1 && cube[i].nbr[j].size()!=4)
{
	printf(ANSI_COLOR_RED"inconsistent four to one balance\n" ANSI_COLOR_RESET) ;
	exit(0);
}
}
}

int count;
for(unsigned int i=0;i<cube.size();i++)
{
			//count same tag
  for(unsigned int j=0;j<6;j++)
{
	count=0;
	
	if(cube[i].nonlocal_nbr.size()!=0)
	{
          //count nume neighbors
          for(unsigned int k=0;k<cube[i].nonlocal_nbr.size();k++)
          {
			  if(cube[i].nonlocal_nbr.at(k).face_tag==j)
			  {
				  count++;
			  }
		  }	 
	  	  
	}
	if(count==4)
	{
	printf("count %d\n",count);
    }
  }
  
  if(count!=0 && count!=1 && count!=4)
  {
	 printf(ANSI_COLOR_RED"inconsistent four to one balance\n" ANSI_COLOR_RESET) ; 
  }
  
  if(count!=0 && count!=1)
  {
	  //printf("%d\n",count);
  }
}





}

void GetCubeCoordLimits(Cube &cube,const Center_coords &XYZ,double *ancestor_length,unsigned int cube_id,double *xyz);
void CalculateArea(double *xyz, double *normal_sum);
void CheckWaterTight(int my_rank,Cube &cube,const Center_coords  &XYZ,double *ancestor_length,Normal *normals)
{
	
double xyz[6];

double normal_sum[6]={0};	
//double sum[3]={0.0,0.0,0.0};
unsigned int nbr_id;
unsigned int  co_face_no[6]={1,0,3,2,5,4};

int idx[4]={0,3,4,7};
unsigned int elem_id;
	
double *sum=NULL;
sum=(double*)calloc(3*cube.size(),sizeof(double));	
	
for(unsigned int i=0;i<cube.size();i++)

//for(unsigned int i=0;i<8;i++)
{
	//unsigned int i=0;
	/*
	   sum[0]=0.0;
		sum[1]=0.0;
		sum[2]=0.0;
		*/
	//elem_id=idx[i];
	elem_id=i;	
		
	//if(cube[elem_id].nbr[0].at(0)>=0 && cube[elem_id].nbr[1].at(0)>=0 && cube[elem_id].nbr[2].at(0)>=0 && cube[elem_id].nbr[3].at(0)>=0 && cube[elem_id].nbr[4].at(0)>=0 && cube[elem_id].nbr[5].at(0)>=0)
	{	
		/*
		if(my_rank==0)
		{
			printf("elem_id=%lu\n",elem_id);
		}
		*/
		
	for(int j=0;j<6;j++)
	{	
	if(cube.at(elem_id).nbr[j].at(0)>=0)
	{			
	for(unsigned int k=0;k<cube.at(elem_id).nbr[j].size();k++)
	{	
	 nbr_id=cube.at(elem_id).nbr[j].at(k);
	 		 	
	 normal_sum[0]=normals[nbr_id].face_normal[co_face_no[j]].nx;
	 normal_sum[1]=normals[nbr_id].face_normal[co_face_no[j]].ny;
	 normal_sum[2]=normals[nbr_id].face_normal[co_face_no[j]].nz;
	 
	 
	 // watch out for the multiple neighbors 
	 
	 if(cube.at(nbr_id).nbr[co_face_no[j]].size()==4)
	 {
     normal_sum[0]=0.25*normal_sum[0];
	 normal_sum[1]=0.25*normal_sum[1];
	 normal_sum[2]=0.25*normal_sum[2];
		 
	 }
	 		
	 /*
	 if(my_rank==0)
	 {
	 printf("nbr_id %d %lf %lf %lf\n",nbr_id,normal_sum[0],normal_sum[1],normal_sum[2]);		
     }
     */
	 sum[3*i+0]=sum[3*i+0]+normal_sum[0];
	 sum[3*i+1]=sum[3*i+1]+normal_sum[1];
	 sum[3*i+2]=sum[3*i+2]+normal_sum[2];
			
    }
     }
     /*
     else
     {
		 	 	
	 normal_sum[0]=normals[elem_id].face_normal[j].nx;
	 normal_sum[1]=normals[elem_id].face_normal[j].ny;
	 normal_sum[2]=normals[elem_id].face_normal[j].nz;		
	
	 sum[3*i+0]=sum[3*i+0]-normal_sum[0];
	 sum[3*i+1]=sum[3*i+1]-normal_sum[1];
	 sum[3*i+2]=sum[3*i+2]-normal_sum[2];
	 }
	 */
	
	}
	
  }	
}
	//printf(ANSI_COLOR_GREEN"surface area %d %lf  %lf %lf \n" ANSI_COLOR_RESET,my_rank,sum[0],sum[1],sum[2]);
	/*
	for(unsigned int i=0;i<cube.size();i++)
{
	if(cube[i].nbr[0].at(0)>=0 && cube[i].nbr[1].at(0)>=0 && cube[i].nbr[2].at(0)>=0 && cube[i].nbr[3].at(0)>=0 && cube[i].nbr[4].at(0)>=0 && cube[i].nbr[5].at(0)>=0)
{
	printf(ANSI_COLOR_GREEN "my_rank =%d %d %16.16lf %16.16lf %16.16lf\n"ANSI_COLOR_RESET,my_rank,i,sum[3*i+0],sum[3*i+1],sum[3*i+2]);
}
}
*/ 
int pass=1;
for(unsigned int i=0;i<cube.size();i++)
{
	if(cube[i].nbr[0].at(0)>=0 && cube[i].nbr[1].at(0)>=0 && cube[i].nbr[2].at(0)>=0 && cube[i].nbr[3].at(0)>=0 && cube[i].nbr[4].at(0)>=0 && cube[i].nbr[5].at(0)>=0)
{
	
	if(fabs(sum[3*i])+fabs(sum[3*i+1])+fabs(sum[3*i+2])>1.e-12)
	{
		printf(ANSI_COLOR_RED "inconsistency in local connectivity %16.16lf %16.16lf %16.16lf"ANSI_COLOR_RESET,sum[0],sum[1],sum[2]);
		pass=0;
		exit(0);
	}
}

}
	if(pass)
	{
		printf("*****************************\n");
		printf("Test for Closedness of the Local Volumes Passed Successfully\n");
		printf("*****************************\n");
	}
	
	
	
	
}



void GetCubeCoordLimits(Cube &cube,const Center_coords &XYZ,double *ancestor_length,unsigned int cube_id,double *xyz)
{
	
	  unsigned int level=cube.at(cube_id).level;
	  unsigned int coord_index=cube.at(cube_id).centeroid_index;
	 
	  double dl,denum;
	  denum=pow(2,level);
	  dl=ancestor_length[0]/denum;
	  xyz[0]=XYZ.at(coord_index).x-dl*0.5;
      xyz[1]=XYZ.at(coord_index).x+dl*0.5;
	
	  dl=ancestor_length[1]/denum;
      xyz[2]=XYZ.at(coord_index).y-dl*0.5;
      xyz[3]=XYZ.at(coord_index).y+dl*0.5;
	  
	  dl=ancestor_length[2]/denum;
      xyz[4]=XYZ.at(coord_index).z-dl*0.5;
      xyz[5]=XYZ.at(coord_index).z+dl*0.5; 
	
	  
	
}


void CalculateArea(double *xyz, double *normal_sum)
{
	
	// v1=(xyz[1]-xyz[0])i+0  j+0 k 
	// v2=0 i+(xyz[3]-xyz[2]) j +0 k
	// v3=0 i+0 j +(xyz[5]-xyz[4]) k 
	double n1,n2,n3;
	n1=(xyz[1]-xyz[0])*(xyz[3]-xyz[2]);	
	n2=(xyz[1]-xyz[0])*(xyz[5]-xyz[4]);	
	n3=(xyz[3]-xyz[2])*(xyz[5]-xyz[4]);	

	normal_sum[0]=-n1;
	normal_sum[1]=n1;
	
	normal_sum[2]=-n2;
	normal_sum[3]=n2;
	
	normal_sum[4]=-n3;
	normal_sum[5]=n3;
		
}

void CalculateNormals(int my_rank,Cube &cube,const Center_coords &XYZ,double *ancestor_length,Normal *normals)
{
	
	double normal_sum[6];
	double xyz[6];
	
	for(unsigned int i=0;i<cube.size();i++)
	{
		
		GetCubeCoordLimits(cube,XYZ,ancestor_length,i,xyz);
		CalculateArea(xyz,normal_sum);		
		
		normals[i].face_normal[0].nx=0.0;
		normals[i].face_normal[0].ny=0.0;
		normals[i].face_normal[0].nz=normal_sum[0];
		
		normals[i].face_normal[1].nx=0.0;
        normals[i].face_normal[1].ny=0.0;
        normals[i].face_normal[1].nz=normal_sum[1];
        
		normals[i].face_normal[2].nx=0.0;
		normals[i].face_normal[2].ny=normal_sum[2];
		normals[i].face_normal[2].nz=0.0;
		
		normals[i].face_normal[3].nx=0.0;		
		normals[i].face_normal[3].ny=normal_sum[3];
		normals[i].face_normal[3].nz=0.0;
				
		normals[i].face_normal[4].nx=normal_sum[4];
		normals[i].face_normal[4].ny=0.0;
        normals[i].face_normal[4].nz=0.0;
        
		normals[i].face_normal[5].nx=normal_sum[5];
		normals[i].face_normal[5].ny=0.0;
        normals[i].face_normal[5].nz=0.0;
		
	}
	/*
	if(my_rank==0)
	{
	for(unsigned int i=0;i<cube.size();i++)
	{
		printf("cube%d\n",i);
		
	for(int j=0;j<6;j++)
	{
		printf("%d %lf  %lf %lf \n",j,normals[i].face_normal[j].nx,normals[i].face_normal[j].ny,normals[i].face_normal[j].nz);
	}
    }
  }
	*/
}

//===============================================================
//
//                    CRS format construction
//
//===============================================================

void ConstructCRSFormat(int my_rank,Cube &cube,unsigned int **ia,unsigned int **ja, MPI_Comm Comm,int offset, int ncube_total,int np)
{

unsigned int size=cube.size();

(*ia)=(unsigned int*)calloc((size+1),sizeof(unsigned int));

unsigned int mysize;

for(unsigned int i=0;i<cube.size();i++)
{
mysize=0;

 for(unsigned int j=0;j<6;j++)
{
if(cube.at(i).nbr[j].at(0)>=0)
{
mysize=mysize+cube.at(i).nbr[j].size();
}
}
mysize=mysize+cube.at(i).nonlocal_nbr.size();
(*ia)[i+1]=mysize;

}

for(unsigned int i=0;i<size;i++)
{
(*ia)[i+1]=(*ia)[i+1]+(*ia)[i];
}


unsigned int size2;

(*ja)=(unsigned int*)calloc((*ia)[size],sizeof(unsigned int));

unsigned int *vtx_slist=NULL;
unsigned int *vtx_rlist=NULL;

vtx_slist=new unsigned int[np+1];
vtx_rlist=new unsigned int[np+1];

for(int i=0;i<np;i++)
{
vtx_slist[i]=(*ia)[size];
}

MPI_Alltoall(vtx_slist, 1, MPI_UNSIGNED, vtx_rlist, 1, MPI_UNSIGNED, Comm);

for(int i=0;i<np;i++)
{
vtx_slist[i+1]=vtx_slist[i+1]+vtx_slist[i];
}

printf("ia size %d\n",(*ia)[size]);








for(int i=0;i<np;i++)
{
vtx_slist[i]=offset;
}

MPI_Alltoall(vtx_slist, 1, MPI_UNSIGNED, vtx_rlist, 1, MPI_UNSIGNED, Comm);

vtx_rlist[np]=ncube_total;

printf("ia size %d\n",(*ia)[size]);


unsigned int count=0;

for(unsigned int i=0;i<cube.size();i++)
{
 for(unsigned int j=0;j<6;j++)
{
if(cube.at(i).nbr[j].at(0)>=0)
{
// add locals
for(unsigned int k=0;k<cube.at(i).nbr[j].size();k++)
{
  (*ja)[count]=cube.at(i).nbr[j].at(k)+offset;
count++;
}
}
}
// add nonlocals
for(unsigned int k=0;k<cube.at(i).nonlocal_nbr.size();k++)
{
 (*ja)[i+1]=cube.at(i).nonlocal_nbr.at(k).elem_id+vtx_rlist[cube.at(i).nonlocal_nbr.at(k).proc_id];
count++;
}

}

//
//
// gather all the data in one processor
//
//

for(int i=0;i<np;i++)
{
vtx_slist[i]=(*ia)[size];
}

MPI_Alltoall(vtx_slist, 1, MPI_UNSIGNED, vtx_rlist, 1, MPI_UNSIGNED, Comm);


unsigned int max_size=0;

for(int i=0;i<np;i++)
{
max_size+=vtx_rlist[i];
}

printf("%lu\n",max_size);

//if(my_rank==0)
{

unsigned int* iat=NULL;
unsigned int* jat=NULL;
iat=(unsigned int*)calloc((ncube_total+1),sizeof(unsigned int)); 
jat=(unsigned int*)calloc(max_size,sizeof(unsigned int)); 

}









delete[] vtx_slist;
delete[] vtx_rlist;


}

/*
unsigned int level,coord_index;
double Xa,Xb,Ya,Yb,Za,Zb;
 
  double denum;

      level=cube.at(0).level;
      coord_index=cube.at(0).centeroid_index;
      denum=pow(2,level);
      dx=ancestor_length[0]/denum;
      dy=ancestor_length[1]/denum;
      dz=ancestor_length[2]/denum;
      Xa=XYZ.at(coord_index).x-dx*0.5;
      Xb=XYZ.at(coord_index).x+dx*0.5;
      Ya=XYZ.at(coord_index).y-dy*0.5;
      Yb=XYZ.at(coord_index).y+dy*0.5;
      Za=XYZ.at(coord_index).z-dz*0.5;
      Zb=XYZ.at(coord_index).z+dz*0.5;


 
printf("CENTER %lf %lf %lf\n",XYZ.at(0).x,XYZ.at(0).y,XYZ.at(0).z);

printf("dxyz1 %lf %lf %lf denum =%lf\n",dx,dy,dz,denum);

printf("coords %lf %lf %lf %lf %lf %lf\n",Xa,Xb,Ya,Yb,Za,Zb);
*/
