#include "typedefs.h"
#define nvar 1


//======================================================================
//
//                Compact Version of Adding Kids
//
//======================================================================
void check_coordinates(Cube&cube);

void find_adjacent_face(Cube& cube,Cube& cube_temp,int id, int tag, int face_id, int *neigbor_id);

void AMR_Refine_Cube(Cube& cube,int id,int L, int M, int N)
{
 unsigned int i,j,k,l1;	
 unsigned int nquad=cube.size();
 unsigned int index=nquad;	

// define aa temporary cube and save current object to a temporary object 
 // note that adaptation level increase by one 
 Cube cube_temp;
 cube_temp.push_back(cube_data());
 cube_temp[0].level=cube[id].level;
 
 for(i=0;i<6;i++){
	cube_temp[0].xyz[i]=cube[id].xyz[i];
    }
	
  for(i=0;i<6;i++){
	 for(j=0;j<cube[id].nbr[i].size();j++){
	cube_temp[0].nbr[i].push_back(cube[id].nbr[i].at(j));
       }
   }
   unsigned int q_size=L*M*N*nvar;
/*   cube_temp[0].q=(double*)malloc(q_size*sizeof(double));
   
   for(j=0;j<q_size;j++)
		{
	     cube_temp[0].q[i]=1.0;
       }
   */
// add 7 elements to the number of cubes
    
 for(i=0;i<7;i++) {
  cube.push_back(cube_data());
 } 
 
 //  set adaptation level,
 cube[id].level= cube_temp[0].level+1;
 
  for(j=0;j<7;j++)
		{
	     cube[nquad+j].level=cube_temp[0].level+1;
	     //cube[nquad+j].q=(double*)malloc(q_size*sizeof(double));     
       }
 
 // eliminate the old neighbor data from the parent element since now 
 //parent will be ith child (kind of like each element pointinmg to itself )
 
 for(i=0;i<6;i++)
 {
 cube[id].nbr[i].clear();
 }
//======================================================================
// set the coordinate bounds for each cube

double xmin=cube_temp[0].xyz[0];  	
double ymin=cube_temp[0].xyz[2];
double zmin=cube_temp[0].xyz[4];

double xmax=cube_temp[0].xyz[1];  	
double ymax=cube_temp[0].xyz[3];
double zmax=cube_temp[0].xyz[5];

double xmid=(xmin+xmax)*0.5;
double ymid=(ymin+ymax)*0.5;
double zmid=(zmin+zmax)*0.5;
 
// Making use of mappings helps of write one code and execute separately 
// for each child cube, map_index refers to element number

unsigned int map_index[8];
map_index[0]=id;

for(i=1;i<8;i++)
{
	map_index[i]=nquad+(i-1);
	//printf("map index %d\n",map_index[i]);
}

double X[3]={xmin,xmid,xmax};
double Y[3]={ymin,ymid,ymax};
double Z[3]={zmin,zmid,zmax};

unsigned int Xmin_idx[8]={0,1,1,0,0,1,1,0};
unsigned int Xmax_idx[8]={1,2,2,1,1,2,2,1};

unsigned int Ymin_idx[8]={0,0,1,1,0,0,1,1};
unsigned int Ymax_idx[8]={1,1,2,2,1,1,2,2};

unsigned int Zmin_idx[8]={0,0,0,0,1,1,1,1};
unsigned int Zmax_idx[8]={1,1,1,1,2,2,2,2};

for(i=0;i<8;i++)
{
  cube[map_index[i]].xyz[0]=X[Xmin_idx[i]];  	
  cube[map_index[i]].xyz[1]=X[Xmax_idx[i]];
  
  cube[map_index[i]].xyz[2]=Y[Ymin_idx[i]];
  cube[map_index[i]].xyz[3]=Y[Ymax_idx[i]];
  
  cube[map_index[i]].xyz[4]=Z[Zmin_idx[i]];
  cube[map_index[i]].xyz[5]=Z[Zmax_idx[i]];
 
}
//======================================================================
// first loop over the neighbors of the cube_temp, at each edge that 
  // the cube_temp has only one neighbor add the new elements
  // remove the id from the neighborhood, this adds to the computaion 6 addition of vectors
  // but makes coding easier and more readable, not a good idea I prefer to overwrite instead of removing and adding
   
 // int face_no[6]={0,1,2,3,4,5};
 int co_face_no[6]={1,0,3,2,5,4};
    
 for(i=0;i<6;i++)
  {
	for(j=0;j<cube_temp[0].nbr[i].size();j++)
	{
	if(cube_temp[0].nbr[i].at(j)>-1)
   {
		for(k=0;k<cube[cube_temp[0].nbr[i].at(j)].nbr[co_face_no[i]].size();k++)
		{
			if(cube[cube_temp[0].nbr[i].at(j)].nbr[co_face_no[i]].at(k)==id)
			{
				cube[cube_temp[0].nbr[i].at(j)].nbr[co_face_no[i]].erase(cube[cube_temp[0].nbr[i].at(j)].nbr[co_face_no[i]].begin()+k);
			}
	    }
	}     
  }
}

// insert all the new elements to as neighbors to neighboring elements

/*
for(i=0;i<cube.size();i++)
{
	printf("============================ %d\n",i);
	
	for(j=0;j<6;j++)
	{
		for(k=0;k<cube[i].nbr[j].size();k++)
		{
			printf("%d\t",cube[i].nbr[j].at(k));
		}
		printf("\n");
	}
	
}
*/

//======================================================================
// Assign interior faces for each element 
//======================================================================
 
unsigned int I[24]={1,3,5,1,3,4,1,2,4,1,2,5,0,3,5,0,3,4,0,2,4,0,2,5}; 

unsigned int nbr_index[24]={nquad+3,nquad+2,nquad,nquad+4,nquad+1,id,nquad+5,nquad,nquad+2,nquad+6,id,nquad+1,id,nquad+6,nquad+4,nquad,nquad+5,nquad+3,nquad+1,nquad+4,nquad+6,nquad+2,nquad+3,nquad+5};

//unsigned int nbr_index[24]={map_index[4],map_index[3],map_index[1],map_index[5],map_index[2],map_index[0],map_index[6],map_index[1],map_index[3],map_index[7],map_index[0],map_index[2],map_index[0],map_index[7],map_index[5],map_index[1],map_index[6],map_index[4],map_index[2],map_index[5],map_index[7],map_index[3],map_index[4],map_index[6]};

for(i=0;i<8;i++)
{ 
	for(j=0;j<3;j++)
	{
		cube[map_index[i]].nbr[I[3*i+j]].push_back(nbr_index[3*i+j]);
		//printf(" %d %d\t",I[3*i+j],nbr_index[3*i+j]);
	}	
}


  // update the kid neighbor, neigbors z=0, y=0 and x=0 dont change
  // update the other 3 neighbors
  // sit in the middle of the element, 
  // 0 >> bottom 1>> top 
  // 2 >> right  3>>left
  // 4 >> back   5>> front
  // some of the neghbirs will be inherited from the parent element
  
  unsigned int elem_id;
  
   // now there are only 3 faces exposed, the other three were assigned before
   // same mapping for upper elements only the dfirst index changes
  
  I[0]=0,I[1]=2,I[2]=4;
  I[3]=0,I[4]=2,I[5]=5;
  I[6]=0,I[7]=3,I[8]=5;
  I[9]=0,I[10]=3,I[11]=4;
  I[12]=1,I[13]=2,I[14]=4;
  I[15]=1,I[16]=2,I[17]=5;
  I[18]=1,I[19]=3,I[20]=5;
  I[21]=1,I[22]=3,I[23]=4;   
  
  int idx;
  int a,b;
  
  // index2 is generated to eliminate 3 if's
  
  int index2[6]={0,0,1,1,2,2};
  // bottom and top planes both faces need xy coordinates (index2[0] and index2[1]) >>> tag[0]
  //left and right need xz (index2[2] and index2[3]) >>> tag[1]
  // back and front faces need yz (index2[4] and index2[5]) >>> tag[2]
  
 for(i=0;i<8;i++)
  {
	  elem_id=map_index[i];
	  for(j=0;j<3;j++)
	  {
		   a=index2[I[3*i+j]];
		   b=I[3*i+j];		   
		   find_adjacent_face(cube,cube_temp,elem_id,a,b,&idx);
		  // printf("elem_id %d  face id =%d idx=%d\n",elem_id,b,idx);
		   cube[elem_id].nbr[b].push_back(idx);
	  
	  // update the neighbor if neighbor exists
	  	  if(idx>-1)
	  	  {	   
		    cube[idx].nbr[co_face_no[b]].push_back(elem_id);	       
	      }
	      
	      }          
    }	  
  
  
//======================================================================

// Erase cube_temp's nbrs for future use

//======================================================================

//free(cube_temp[0].q);

for(i=0;i<6;i++)
{
	while(cube_temp[0].nbr[i].size()!=0)
	{
	cube_temp[0].nbr[i].pop_back();
    }
}
	
}


//======================================================================
void find_adjacent_face(Cube& cube,Cube& cube_temp,int id, int tag, int face_id, int *neigbor_id)
{
	 unsigned int i,idx;
	 unsigned int I[2];
	 unsigned int J[2];
	 
	switch(tag)
	{
		case 0:
		 I[0]=0;
		 I[1]=1;
		 J[0]=2;
		 J[1]=3;
		 break;
		case 1:
		 I[0]=0;
		 I[1]=1;
		 J[0]=4;
		 J[1]=5;
		 break;
		case 2:
		 I[0]=2;
		 I[1]=3;
		 J[0]=4;
		 J[1]=5;
		 break;
		
	}
	double xc_id,yc_id,r;
	double x0,x1,y0,y1,xc,yc;

	     xc_id=0.5*(cube[id].xyz[I[0]]+cube[id].xyz[I[1]]);		 
		 yc_id= 0.5*(cube[id].xyz[J[0]]+cube[id].xyz[J[1]]);  
		 	 
		
		
		//printf("coordinates current element =%lf %lf \n",xc_id,yc_id);
		if(face_id<0)
		{
		printf("???????????????????????????????????%d\n",face_id);
	}
	
	#if(DEBUG)
	if(cube_temp[0].nbr[face_id].size()>4)
	{
		printf("Failure in satistying 4:1 ration\n");
		printf("Error  Code Excited in find_adjacent_face function\n");
		exit(0);
	}
	#endif
	
	 if(cube_temp[0].nbr[face_id].size()==1)	
		{
			*neigbor_id=cube_temp[0].nbr[face_id].at(0);
         }
	else
	{		 
		r=1000.0;
	 for(i=0;i<cube_temp[0].nbr[face_id].size();i++)
	 {
		 x0= cube[cube_temp[0].nbr[face_id].at(i)].xyz[I[0]];
		 x1= cube[cube_temp[0].nbr[face_id].at(i)].xyz[I[1]]; 
		 
		 y0= cube[cube_temp[0].nbr[face_id].at(i)].xyz[J[0]];
		 y1= cube[cube_temp[0].nbr[face_id].at(i)].xyz[J[1]]; 
		 
		 xc=0.5*(x0+x1);
		 yc=0.5*(y0+y1);
		 //  printf("distance=%lf \n",sqrt(r));
		 if((xc-xc_id)*(xc-xc_id)+(yc-yc_id)*(yc-yc_id)<r)
		 {
			idx=i;
			r=(xc-xc_id)*(xc-xc_id)+(yc-yc_id)*(yc-yc_id); 
		 }		 		 
	 }
	 // printf("distance_min= %lf nbr_size=%d\n",sqrt(r),cube_temp[0].nbr[face_id].size());
	 *neigbor_id=cube_temp[0].nbr[face_id].at(idx);
	 
	 #if(DEBUG)
	 if(sqrt(r)>1.e-6)
	 {
		 printf("the error in face centers might be large\n");
		 sleep(1);
	 }
	 #endif
     }
	 
	 //printf("nbr_id in adjacent face function %d\n",cube_temp[0].nbr[face_id].at(idx));  
	 
 }
//======================================================================


// std::sort (myvector.begin(), myvector.end(), myfunction)


//======================================================================


// enforces the 
void check_list(Vector_Int& refine_list,int elem_id, int *T);

void Amr_four_to_one_enforce(Cube& cube,Vector_Int& refine_list,int L,int M, int N)
{
	
	// dynamically check and seee what happens to the 
	// element if I refine the element, if level offset is more than one adapt that one
	// as well
	
	unsigned int elem_id,i,j,k,l1;
	int nbr_id;
	Vector_Int adapt_level; 
		
	unsigned int list_size;	
		
int a=1;
int T;
 int istart=0;	
 int iend=refine_list.size();
//printf("istart %d iend %d\n",istart,iend);

int count=0;	
	
#if(1)
 while(a==1)
 {		
	a=0;
	
	printf("istart %d iend %d\n",istart,iend);
			
	 for(i=istart;i<iend;i++)
	 {
		elem_id=refine_list.at(i);
		//printf("===================elem_id %d\n",elem_id);
		
		for(j=0;j<6;j++)
		{
			for(k=0;k<cube[elem_id].nbr[j].size();k++)
			{
				nbr_id=cube[elem_id].nbr[j].at(k);
				//exclude the boundaries 
				//printf("===========%d\n",nbr_id);
				
				if(nbr_id>-1)
				{
					//printf("===========%d %d %d\n",nbr_id,cube[nbr_id].level,cube[elem_id].level);
				if(cube[(unsigned)nbr_id].level<cube[elem_id].level)
				{
					// refine it
					//
					//refine_list.push_back(nbr_id);
					//printf("nbr_id=%d\n",nbr_id);
					
					check_list(refine_list,nbr_id,&T);
					// check the original list, if this, set is as -1  
				  	if(T==0)
				  	{
							refine_list.push_back(nbr_id);
							//printf("nbr_id=%d\n",nbr_id);
							a=1;
					}
				}
			 }		
			}	
		}
		
	}
	
	/*
	for(i=istart;i<iend;i++)
	{
		elem_id=(unsigned)(refine_list.at(i));
		//printf("%d\n",elem_id);
		AMR_Refine_Cube(cube,elem_id,L,M,N);
	 }
	 */
	 
	 if(a==1)
    {
	  istart=iend;
	  iend=refine_list.size();  
    }
	
	// printf("istart %d iend %d\n",istart,iend);
	 
	 /*
	 for(i=istart;i<iend;i++)
	{
		printf("refine_list=%d\n",refine_list.at(i));
	 }
	 */
 }
#endif

//std::sort (refine_list.begin(), refine_list.end(), compare_level);

unsigned int level_max=1;

for(i=0;i<refine_list.size();i++)
{
	if(cube[refine_list.at(i)].level>0)
	{
		level_max=cube[refine_list.at(i)].level;
	}
}
//printf("level_max=%d\n",level_max);

for(i=0;i<=level_max;i++)
{
	for(j=0;j<refine_list.size();j++)
	{
		//printf("-------------------------%d %d\n",refine_list.at(j),cube[refine_list.at(j)].level);
		if(cube[refine_list.at(j)].level==i)
		{
		elem_id=(unsigned)(refine_list.at(j));
		//printf("-------------------------%d\n",elem_id);
		AMR_Refine_Cube(cube,elem_id,L,M,N);
		refine_list.at(j)=-1;
	    }
	}
	
}	

//printf("size %d\n",refine_list.size());

 while(refine_list.size()!=0)
 {
	 refine_list.pop_back();
 }

#if(DEBUG)
    printf("*******************************************************************\n");
    printf("\n");
	printf("                    4 to 1 Enforcing Check ... \n"); 
	printf("\n");
	printf("\n"); 

for(i=0;i<cube.size();i++)
{
	//printf("hi\n");
	for(j=0;j<6;j++)
	{
		if(cube[i].nbr[j].size()>4)
		{
			printf("4:1 ration was not satsified\n");
			printf("Element at fault is %d\n",i);
			printf("\n");
	        printf("*******************************************************************\n"); 
			exit(0);
		}
	}
}

   
	printf("\n"); 
	printf("                    4 to 1 Enforcing Successful\n"); 
	printf("\n");
	printf("*******************************************************************\n"); 
	
#endif


}
//======================================================================

void check_list(Vector_Int& refine_list,int elem_id, int *T)
{
	int i;
	
	*T=0;
	
	for(i=0;i<refine_list.size();i++)
	{
		if(elem_id==refine_list.at(i))
		{
			*T=1;
		}
	}	
}


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
	printf("                   Element Connectivity Check ... \n"); 
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
				if(cube[nbr_id].nbr[l1].at(m)==i)
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
	printf("               Element Connectivity Check Successful \n"); 
	printf("\n");
	printf("*******************************************************************\n"); 
	
	
}
//======================================================================
void check_coordinates(Cube&cube)
{
	int i;
	
	printf("*******************************************************************\n");
	printf("\n");
	printf("\n"); 
	printf("                   Coordinate Bounds Check ... \n"); 
	
	for(i=0;i<cube.size();i++)
	{
		if(cube[i].xyz[0]>cube[i].xyz[1] || cube[i].xyz[2]>cube[i].xyz[3] || cube[i].xyz[4]>cube[i].xyz[5])
		{
			printf("error in corrdinate bounds in function check_coordinates\n");
			exit(0);
		}
	}
	
	
	printf("\n"); 
	printf("               Coordinate Bounds Check Successful \n"); 
	printf("\n");
	printf("*******************************************************************\n"); 
	
}


//======================================================================
// n incules the number of points in each direction
void sphere(Cube& cube, Vector_Int& refine_list,double C)
{

double xc,yc,zc;
unsigned int i;

for(i=0;i<cube.size();i++)
 {
	xc=0.5*(cube[i].xyz[0]+cube[i].xyz[1]);
	yc=0.5*(cube[i].xyz[2]+cube[i].xyz[3]);
	zc=0.5*(cube[i].xyz[4]+cube[i].xyz[5]);
	
	if(xc*xc+yc*yc+zc*zc<C)
	//if(xc*xc+yc*yc<C)
	{
		refine_list.push_back(i);
	}
 }	
		
}

//======================================================================
//======================================================================
void circle(double **x, double **y,int n)
{

double xc,yc,zc;
unsigned int i;

double theta;
(*x)=(double*)malloc(n*sizeof(double));
(*y)=(double*)malloc(n*sizeof(double));

double r=0.5;
double pi=4.*atan(1.0);
for(i=0;i<n;i++)
 {
	theta=2.*pi/(n-1)*i;
	(*x)[i]=r*cos(theta);
	(*y)[i]=r*sin(theta);
 }	
		
}

//======================================================================

void identify(Cube& cube, Vector_Int& refine_list,double *x,double *y,int n)
{

double xc,yc,zc;
unsigned int i,j;
int a,b;

for(i=0;i<cube.size();i++)
 {
	 for(j=0;j<n;j++)
	 {
	a=x[j]> cube[i].xyz[0] && x[j]< cube[i].xyz[1];
	b=y[j]> cube[i].xyz[2] && y[j]< cube[i].xyz[3];
	if(a && b)
	{
		refine_list.push_back(i);
		break;
	}
   }
 }	
		
}


//======================================================================


//======================================================================
//
//                     Remove Kids
//
//======================================================================

void AMR_Derefine_Cube(Cube& cube,int parent_id)
{
	unsigned int nquad=cube.size();	
	Cube cube_temp;
	cube_temp.push_back(cube_data());
	unsigned int i,j,k,l1,elem_id;
	int nbr_id;
	
	unsigned int map_index[8];
    map_index[0]=parent_id;

//for(i=1;i<8;i++)
{
	//map_index[i]=parent_id+(i-1);
	i=1;
	map_index[i]=8;
	map_index[i+1]=9;
	map_index[i+2]=10;
	map_index[i+3]=11;
	map_index[i+4]=12;
	map_index[i+5]=13;
	map_index[i+6]=14;
	
	
}

for(i=0;i<8;i++)
{
	if(i!=0)
	{
		map_index[i]=map_index[i]+7;
	}
	printf("map index %d\n",map_index[i]);
}


	int co_face_no[6]={1,0,3,2,5,4};
	
	
	
	
	unsigned int I[24];
	
  I[0]=0,I[1]=2,I[2]=4;
  I[3]=0,I[4]=2,I[5]=5;
  I[6]=0,I[7]=3,I[8]=5;
  I[9]=0,I[10]=3,I[11]=4;
  I[12]=1,I[13]=2,I[14]=4;
  I[15]=1,I[16]=2,I[17]=5;
  I[18]=1,I[19]=3,I[20]=5;
  I[21]=1,I[22]=3,I[23]=4;  
  
 // update xyz for parent_id'th element
// indices 0, 1,3   are the same

  cube[parent_id].xyz[1]=cube[map_index[1]].xyz[1];
  cube[parent_id].xyz[3]=cube[map_index[2]].xyz[3];	
  cube[parent_id].xyz[5]=cube[map_index[4]].xyz[5]; 

// saves all the neighbors of the refined element to a temporary object
	int bol=1;

	for(i=0;i<8;i++)
	{
		 elem_id=map_index[i];
		 
		for(j=0;j<3;j++)
	 {
		 		  
		 for(k=0;k<cube[elem_id].nbr[I[3*i+j]].size();k++)
		 {
			 bol=1;
			 for(l1=0;l1<cube_temp[0].nbr[I[3*i+j]].size();l1++)
			 {
				 if(cube_temp[0].nbr[I[3*i+j]].at(l1)==cube[elem_id].nbr[I[3*i+j]].at(k))
				 {
					 bol=0;					 
					 break;
				 }
			 }
			 if(bol==1)
			 {
			 cube_temp[0].nbr[I[3*i+j]].push_back(cube[elem_id].nbr[I[3*i+j]].at(k));
		     }
		 //printf("bol=%d\n",bol);
		 //cube_temp[0].nbr[I[3*i+j]].push_back(cube[elem_id].nbr[I[3*i+j]].at(k));
		 }
   }
 }
 
 // reduce the parent adaptation level by one due to derefinement 
 
  cube[parent_id].level= cube[parent_id].level-1;
	/*
	for(i=0;i<6;i++)
	{
		for(j=0;j<cube_temp[0].nbr[i].size();j++)
		{
			printf("%d\t",cube_temp[0].nbr[i].at(j));
		}
		printf("\n");
	}
		*/
//======================================================================	
// remove the eliminated elements from the neighbors neighborhood data
		
	for(i=0;i<8;i++)
	{
		 elem_id=map_index[i];
		 
		for(j=0;j<3;j++)
	 {
		 for(k=0;k<cube_temp[0].nbr[I[3*i+j]].size();k++)
		 {
			 	 printf("nbr_id %d\n",cube_temp[0].nbr[I[3*i+j]].at(k));
			 	 nbr_id=cube_temp[0].nbr[I[3*i+j]].at(k);		 			 
			 if(nbr_id>-1)
			 {			 		 	
			 	printf("nbr_id = %d =%d\n",nbr_id,cube[nbr_id].nbr[co_face_no[I[3*i+j]]].size());	
			 	
			for(l1=0;l1<cube[nbr_id].nbr[co_face_no[I[3*i+j]]].size();l1++)
			{					 
			  if(cube[nbr_id].nbr[co_face_no[I[3*i+j]]].at(l1)==elem_id)
			  {
		       cube[nbr_id].nbr[co_face_no[I[3*i+j]]].erase(cube[nbr_id].nbr[co_face_no[I[3*i+j]]].begin()+l1);
			  }
		    }
		 }
	    }
	 }
	}
	
//======================================================================	
	
	// add parent_id to the neighbors neighborlood data

	for(i=0;i<6;i++)
	{		 
		 for(j=0;j<cube_temp[0].nbr[i].size();j++)
		 {
			 nbr_id=cube_temp[0].nbr[i].at(j);
			 if(nbr_id>-1)
			 {
			 //printf("nbr_id = %d =%d\n",nbr_id,cube[nbr_id].nbr[co_face_no[i]].size());	
		     cube[nbr_id].nbr[co_face_no[i]].push_back(parent_id);
		    }
		 }
	}
	
// remove the neighborhood connectivity

for(j=0;j<6;j++)
{
	while(cube[map_index[0]].nbr[j].size()>0)
	{
     cube[parent_id].nbr[j].pop_back();
    }
}

for(i=0;i<6;i++)
{
	for(j=0;j<cube_temp[0].nbr[i].size();j++)
	{
     cube[parent_id].nbr[i].push_back(cube_temp[0].nbr[i].at(j));
    }
}

// remove the kids from the connectivity list
for(i=1;i<8;i++)
{
  cube.erase(cube.begin()+map_index[i]);
}

// update the original

for(i=0;i<6;i++)
{
	while(cube_temp[0].nbr[i].size()!=0)
	{
	cube_temp[0].nbr[i].pop_back();
    }
}

//printf("nquad original=  %d number of quads after removal =%d\n",nquad,cube.size());
// now the numbering is going to change depending on which element is derefined first

int max=map_index[7];

for(i=0;i<cube.size();i++)
{
	for(j=0;j<6;j++)
	{
		for(k=0;k<cube[i].nbr[j].size();k++)
		{
			if(cube[i].nbr[j].at(k)>max)
			{
				cube[i].nbr[j].at(k)=cube[i].nbr[j].at(k)-7;
			}
		}
		
	}
}




	
}


//======================================================================

  /*
  cube[id].nbr[1].push_back(nquad+3);
  cube[id].nbr[3].push_back(nquad+2);  
  cube[id].nbr[5].push_back(nquad);

// kid 1
  cube[index].nbr[1].push_back(nquad+4); 
  cube[index].nbr[3].push_back(nquad+1); 
  cube[index].nbr[4].push_back(id); 
// kid 2
  cube[index].nbr[1].push_back(nquad+5); 
  cube[index].nbr[2].push_back(nquad); 
  cube[index].nbr[4].push_back(nquad+2); 
  //kid 3
   cube[index].nbr[1].push_back(nquad+6); 
  cube[index].nbr[2].push_back(id); 
  cube[index].nbr[5].push_back(nquad+1); 
    //kid 4
  cube[index].nbr[0].push_back(id); 
  cube[index].nbr[3].push_back(nquad+6); 
  cube[index].nbr[5].push_back(nquad+4); 
// kid 5
  cube[index].nbr[0].push_back(nquad); 
  cube[index].nbr[3].push_back(nquad+5); 
  cube[index].nbr[4].push_back(nquad+3); 

//kid 6
  cube[index].nbr[0].push_back(nquad+1); 
  cube[index].nbr[2].push_back(nquad+4); 
  cube[index].nbr[4].push_back(nquad+6);
  // kid 7
  cube[index].nbr[0].push_back(nquad+2); 
  cube[index].nbr[2].push_back(nquad+3); 
  cube[index].nbr[5].push_back(nquad+5); 
*/  
