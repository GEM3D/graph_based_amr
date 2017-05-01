/** \mainpage
 *  
 *   "Paralle implementation of AMR for immersed boundary simulations" 
 * 
 *            <br> Part of NSF project: GEM3D 
 *                          
 *
 *             Required Libraries: 
 *               HDF5, 
 *               MPI,
 *               Zoltan, 
 *               CMAKE <br>
 *
 *             Usage:
 *              progName input/geometry *.stl <run.txt
 *              (run needs 4 integers, the first three specify the process topology
                and the last one is the level of adaptation)
 *
 * @authors     Jaber J. Hasbestan, Inanc Senocak <br>
 *  PI:                 Inanc Senocak
 * @date       30 May 2017
 * \details
 * Copyright (c) 2017
 * Mechanical and Bio-medical Engineering Department
 * Boise State University
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */





#include "typedefs.h"
//#include "/home/jhasbestan/GEM_AMR/src/include/typedefs.h"
//#include "${CMAKE_HOME_DIRECTORY}/src/include/typedefs.h"

// define write 0 for points only (Fastest Version, good for debug)
// define write 1 for unstructured representation
// define write 2 for Multi-Block (Slowest Version)
#define WRITE 1
/* ZOLTAN =1 turns partitioning on set to zero to turn it off*/
#define ZOLTAN_ON 0
/* if you set the value to 1 to use geometirc based partitioning */
#define ZOLTAN_GEOMETRIC_PARTITION 0


int main(int argcs, char* pArgs[])
{

  MPI_Init(&argcs,&pArgs);

  static int nface=6;
  static int nchild=8;
  
  unsigned int i,j,k,l;
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  int Dim[3]={0,0,0};
  MPI_Comm Comm_cart;

  // xyz of the mother cube 
  double xyz_ancestor[6];
  xyz_ancestor[0]=-1.0;
  xyz_ancestor[1]=1.0;

  xyz_ancestor[2]=-1.0;
  xyz_ancestor[3]=1.0;

  xyz_ancestor[4]=-1.0;
  xyz_ancestor[5]=1.0;


 unsigned int L=2;
 unsigned int M=L;
 unsigned int N=L;

int max_level;
if(my_rank==0)
{
remove("soln/Pxdmf3d.h5");
remove("soln/Pxdmf3d.xmf");

	printf("Enter the number of processors at each direction\n");
	scanf("%d",&(Dim[0]));
	scanf("%d",&(Dim[1]));
	scanf("%d",&(Dim[2]));
	printf("%d %d %d\n",Dim[0],Dim[1],Dim[2]);
	printf("Enter the number of levels you wish to adapt\n");
	scanf("%d",&max_level);
	/*printf("Define the outer Cube that encloses the entire domain\n");
	
	<scanf("%d",&(xyz_ancestor[0]));
	scanf("%d",&(xyz_ancestor[1]));
	scanf("%d",&(xyz_ancestor[2]));
	scanf("%d",&(xyz_ancestor[3]));
	scanf("%d",&(xyz_ancestor[4]));
	scanf("%d",&(xyz_ancestor[5]));
	*/
}



MPI_Bcast(Dim,3,MPI_INT,0,MPI_COMM_WORLD);
MPI_Bcast(&max_level,1,MPI_INT,0,MPI_COMM_WORLD);

int npx=Dim[0];
int npy=Dim[1];
int npz=Dim[2];

int np;
MPI_Comm_size(MPI_COMM_WORLD,&np); 

if(my_rank==0)
{
printf("total number of procs=%d\n",np);
}

if(my_rank==0)
{
if(np!=(Dim[0]*Dim[1]*Dim[2]))
{
	printf(ANSI_COLOR_RED "number of processors does not match user provided topology\n" ANSI_COLOR_RESET);
	printf(ANSI_COLOR_RED "Please be consistent in assigning\n" ANSI_COLOR_RESET);
	exit(0);
	MPI_Finalize();
}
}



 const int periods[3]={0,0,0};

 MPI_Cart_create(MPI_COMM_WORLD,3,Dim,periods,1,&Comm_cart);

 int coords[3];
 int maxdims=3;

 if(my_rank==0){
   for(int i=0;i<Dim[0]*Dim[1]*Dim[2];i++){
     MPI_Cart_coords(Comm_cart,i,maxdims,coords);
    	printf("%d %d %d\n",coords[0],coords[1],coords[2]);
   }
 }

 // reserve your initial guess of the vector
 Cube cube;
 cube.reserve(500000);
 cube.push_back(cube_data());
// cube[0].nbr[0].reserve(4);
// dx dy and dz are the original size of the single cube

 double dx=xyz_ancestor[1]-xyz_ancestor[0];
 double dy=xyz_ancestor[3]-xyz_ancestor[2];
 double dz=xyz_ancestor[5]-xyz_ancestor[4];
 
// double ancestor_length[3]={dx,dy,dz};

 MPI_Cart_coords(Comm_cart,my_rank,maxdims,coords);


 cube.at(0).level=0;
 cube.at(0).centeroid_index=0;
 
 Center_coords XYZ;
 XYZ.reserve(500000);
 XYZ.push_back(Center_Coords());
 // calculate the center of the parent 
 double xyz[6];
 
 xyz[0]=xyz_ancestor[0]+dx/Dim[0]*coords[0];
 xyz[1]=xyz_ancestor[0]+dx/Dim[0]*(coords[0]+1);
 xyz[2]=xyz_ancestor[2]+dy/Dim[1]*coords[1];
 xyz[3]=xyz_ancestor[2]+dy/Dim[1]*(coords[1]+1);
 xyz[4]=xyz_ancestor[4]+dz/Dim[2]*coords[2];
 xyz[5]=xyz_ancestor[4]+dz/Dim[2]*(coords[2]+1);

//printf("dxyz0 %lf %lf %lf\n",dx,dy,dz);
//printf("endpoints %lf %lf %lf %lf %lf %lf\n",xyz[0],xyz[1],xyz[2],xyz[3],xyz[4],xyz[5]);


double ancestor_length[3]={xyz[1]-xyz[0],xyz[3]-xyz[2],xyz[5]-xyz[4]};

 XYZ.at(0).x=(xyz[0]+xyz[1])*0.5;
 XYZ.at(0).y=(xyz[2]+xyz[3])*0.5;
 XYZ.at(0).z=(xyz[4]+xyz[5])*0.5;

// for debugremove later
//


 // use cartesian coords to construct neighbors
 int coord_nbr[3];
 int nbr_rank;

 // construct the initial nbr connectivity 

 Cart_proc_con cart_proc_con;
 
 Cartesian_Proc_Conn(Dim,coords,Comm_cart,&cart_proc_con);
  
 for(i=0;i<6;i++)
   {
     printf("mapping my_rank %d %d\n",my_rank,cart_proc_con.nbr[i].at(0));
   }
  
 /*
 printf(" my_rank=%d  coordes %d %d %d %d %d %d\n",my_rank,coords[0],coords[1],coords[2],Dim[0],Dim[1],Dim[2]);
*/
/*
 printf("X my_rank %d %lf %lf\n",my_rank,cube[0].xyz[0],cube[0].xyz[1]);
 printf("Y my_rank %d %lf %lf\n",my_rank,cube[0].xyz[2],cube[0].xyz[3]);
 printf("Z my_rank %d %lf %lf\n",my_rank,cube[0].xyz[4],cube[0].xyz[5]);
 // check the volume	
 if(DEBUG)
   {	
     printf("volume[%d]= %lf\n",my_rank,fabs((cube[0].xyz[0]-cube[0].xyz[1])*(cube[0].xyz[2]-cube[0].xyz[3])*(cube[0].xyz[4]-cube[0].xyz[5])));
   }

 printf("size of cube %d\n",cube.size());
*/
 // every processor need to no his offset for parallel I/O, therefore, 
 // it should now that trying to avoid coolective operation for extreme scale computing
 int send=0;

 MPI_Group World_Group;
 MPI_Comm_group(MPI_COMM_WORLD, &World_Group);
 MPI_Group Wait_Group;
 MPI_Status status;
 int source;
 int destination;
 Vector_Unint ranlist;
 //===================================================================
 //
 //                 Now update the neighborhood 
 //
 //===================================================================
 // each proc has one block, its neighbor is either another proc or geometry boundary
 // in either case, all are negative numbers

 for(int i=0;i<nface;i++)
   {
     cube[0].nbr[i].push_back(-i-1);
   }


 Vector_Int bface_id;
 int cube_id=0;
 find_boundary_face(cube,cube_id,bface_id); 
 
 for(int j=0;j<nface;j++)
   {
     printf("my_rank %d before modification %d \n",my_rank,cube[0].nbr[j].at(0));
   }
 
 int index=0;

 for(j=0;j<bface_id.size();j++)
   {
//     printf("-------my_rank %d bface_id %d \n",my_rank,bface_id.at(j));
     
     if(cart_proc_con.nbr[j].at(0)>=0)
       {
	 //	 cube[0].nonlocal_nbr.proc_id.push_back(cart_proc_con.nbr[j].at(0));
	 cube[0].nonlocal_nbr.push_back(Nonlocal_Nbr());

	 cube[0].nonlocal_nbr[index].proc_id=cart_proc_con.nbr[j].at(0);
	 cube[0].nonlocal_nbr[index].face_tag=bface_id.at(j);
	 // no adaptation so far therefore elem_id is 0 at this stage
	 cube[0].nonlocal_nbr[index].elem_id=0;
	 index++;
       }
   }
 
 printf("------------------------------------------------------\n");
 
/*
//if(my_rank==1)
{
 for(j=0;j<nface;j++)
   {
     printf("nbr_id %d \n",cube[0].nbr[j].at(0));
   }

 for(i=0;i<cube[0].nonlocal_nbr.size();i++)
   {
     printf("************nbr_id %d face tag=%d proc_id=%d\n",cube[0].nonlocal_nbr[i].elem_id,cube[0].nonlocal_nbr[i].face_tag,cube[0].nonlocal_nbr[i].proc_id);
   }
}
*/
 //MPI_Barrier(MPI_COMM_WORLD);
 /*
 for(j=0;j<6;j++)
   {
     printf("my_rank %d after modification %d\n",my_rank,cube[0].nbr[j].at(0));
   }
 printf("------------------------------------------------------\n");
 */
//AMR_Refine_Cube(my_rank,MPI_COMM_WORLD,cube,0,L,M,N); 

//void AmrRefineParallel(const int my_rank,const MPI_Comm comm,Cube& cube,const Vector_Int &ordered_refine_list,L,M,N);

Vector_Unint refine_list;
Vector_Unint ordered_refine_list;

refine_list.reserve(9000);
ordered_refine_list.reserve(10000);

#if(0)
if(my_rank==0)
{
	for(unsigned int k=0;k<cube.size();k++)
	{
printf("elem_id =%d elem_level=%d\n",k,cube[k].level);
 /*
 for(j=0;j<nface;j++)
   {
     printf("%d \t",cube[k].nbr[j].at(0));
   }
printf("\n");

 for(i=0;i<cube[k].nonlocal_nbr.size();i++)
   {
     printf("************ %d nbr_id %d face tag=%d proc_id=%d\n",k,cube[k].nonlocal_nbr[i].elem_id,cube[k].nonlocal_nbr[i].face_tag,cube[k].nonlocal_nbr[i].proc_id);
   } 
*/
}

}


printf("**************************************************************\n");
printf("**************************************************************\n");
#endif

MPI_Barrier(MPI_COMM_WORLD);

//printf("+++++++++++++++++++++my_rank =%d ordered list size=%lu\n",my_rank,ordered_refine_list.size());

//======================================================================
//                            debugging
//======================================================================
// debugging for q values in hdf5, no refinement
#if(0)
int rank_id=1;

//if(my_rank==rank_id)
{
	refine_list.push_back(0);
	
}

//refine_list.push_back(0);

AmrNewList(my_rank,MPI_COMM_WORLD,cube,refine_list,ordered_refine_list);

AmrRefineParallel(my_rank,MPI_COMM_WORLD,cube,ordered_refine_list,L,M,N,XYZ,ancestor_length);


if(my_rank==rank_id)
{
	//refine_list.pop_back();
	for(int i=0;i<8;i++)
	{
	refine_list.push_back(i);
    }	
}

AmrNewList(my_rank,MPI_COMM_WORLD,cube,refine_list,ordered_refine_list);
AmrRefineParallel(my_rank,MPI_COMM_WORLD,cube,ordered_refine_list,L,M,N,XYZ,ancestor_length);


//if(my_rank==rank_id)
{
	//refine_list.pop_back();
	ran_list(cube,ranlist,8);
	for(int i=0;i<12;i++)
	{
	refine_list.push_back(ranlist.at(i));
    }	
}

AmrNewList(my_rank,MPI_COMM_WORLD,cube,refine_list,ordered_refine_list);
AmrRefineParallel(my_rank,MPI_COMM_WORLD,cube,ordered_refine_list,L,M,N,XYZ,ancestor_length);

for(unsigned int i=0;i<cube.size();i++)
{
	cube[i].q=my_rank;
}


{
	//refine_list.pop_back();
	ran_list(cube,ranlist,8);
	for(int i=0;i<8;i++)
	{
	refine_list.push_back(ranlist.at(i));
    }	
}

AmrNewList(my_rank,MPI_COMM_WORLD,cube,refine_list,ordered_refine_list);
AmrRefineParallel(my_rank,MPI_COMM_WORLD,cube,ordered_refine_list,L,M,N,XYZ,ancestor_length);

{
	//refine_list.pop_back();
	ran_list(cube,ranlist,8);
	for(int i=0;i<8;i++)
	{
	refine_list.push_back(ranlist.at(i));
    }	
}

AmrNewList(my_rank,MPI_COMM_WORLD,cube,refine_list,ordered_refine_list);
AmrRefineParallel(my_rank,MPI_COMM_WORLD,cube,ordered_refine_list,L,M,N,XYZ,ancestor_length);

{
	//refine_list.pop_back();
	//ran_list(cube,ranlist,8);
	for(unsigned int i=0;i<cube.size();i++)
	{
	//refine_list.push_back(ranlist.at(i));
	refine_list.push_back(i);
    }	
}

AmrNewList(my_rank,MPI_COMM_WORLD,cube,refine_list,ordered_refine_list);
AmrRefineParallel(my_rank,MPI_COMM_WORLD,cube,ordered_refine_list,L,M,N,XYZ,ancestor_length);


{
	//refine_list.pop_back();
	//ran_list(cube,ranlist,8);
	for(unsigned int i=0;i<cube.size();i++)
	{
	//refine_list.push_back(ranlist.at(i));
	refine_list.push_back(i);
    }	
}

AmrNewList(my_rank,MPI_COMM_WORLD,cube,refine_list,ordered_refine_list);
AmrRefineParallel(my_rank,MPI_COMM_WORLD,cube,ordered_refine_list,L,M,N,XYZ,ancestor_length);

{
	//refine_list.pop_back();
	//ran_list(cube,ranlist,8);
	for(unsigned int i=0;i<100;i++)
	{
	//refine_list.push_back(ranlist.at(i));
	refine_list.push_back(i);
    }	
}
/*
AmrNewList(my_rank,MPI_COMM_WORLD,cube,refine_list,ordered_refine_list);
AmrRefineParallel(my_rank,MPI_COMM_WORLD,cube,ordered_refine_list,L,M,N,XYZ,ancestor_length);

*/


#endif






#if(0)
// another evel
if(my_rank==rank_id)
{
	//refine_list.pop_back();
	
/*	ran_list(cube,ranlist,20);
	
	for(int i=0;i<ranlist.size();i++)
{
 refine_list.push_back(ranlist.at(i));
}
*/
	for(int i=0;i<(int)cube.size();i++)
	{
	refine_list.push_back(i);
    }	

//refine_list.push_back(1);
}

AmrNewList(my_rank,MPI_COMM_WORLD,cube,refine_list,ordered_refine_list);
AmrRefineParallel(my_rank,MPI_COMM_WORLD,cube,ordered_refine_list,L,M,N);
//======================================================================

if(my_rank==rank_id)
{
	//refine_list.pop_back();
	
	/*ran_list(cube,ranlist,5);
	
	for(int i=0;i<ranlist.size();i++)
{
 refine_list.push_back(ranlist.at(i));
}
*/
	for(unsigned int i=0;i<cube.size();i++)
	{
	refine_list.push_back(i);
    }	

//refine_list.push_back(1);
}

AmrNewList(my_rank,MPI_COMM_WORLD,cube,refine_list,ordered_refine_list);
AmrRefineParallel(my_rank,MPI_COMM_WORLD,cube,ordered_refine_list,L,M,N);
//======================================================================

// assign your 
#endif
double *geom_xyz=NULL;
int geom_nn;
// airfoil test case
//ReadInGeom(&geom_xyz,&geom_nn);

//double test;
//Pow(1,3,&test);
//printf("TESTSSSSSSSSSSSS%lf\n",test);


ReadSTLGeom(argcs,pArgs,&geom_xyz,&geom_nn,xyz);


Normal *normals=NULL;
 
  
 clock_t beginnewlist,endnewlist,beginrefine,endrefine;
 clock_t start,end;
 
 start=clock();

for(int i=0;i<max_level;i++)
{

IsInsideSolid(cube,XYZ,refine_list,geom_xyz,ancestor_length,geom_nn);
/*
for(unsigned int i=0;i<refine_list.size();i++)
{
	printf("??????????????? %u\n",refine_list.at(i));
}
*/

//beginnewlist=clock();

AmrNewList(my_rank,MPI_COMM_WORLD,cube,refine_list,ordered_refine_list);

//printf("AMR_NEW_List Successful\n");

//clear initial list

refine_list.clear();
//refine_list.shrink_to_fit();

//endnewlist=clock();

//beginrefine=clock();

AmrRefineParallel(my_rank,MPI_COMM_WORLD,cube,ordered_refine_list,L,M,N,XYZ,ancestor_length);

// clear new list
ordered_refine_list.clear();
//ordered_refine_list.shrink_to_fit();


//endrefine=clock();

//elem_conn_check(cube);

//CheckConnectivityConsistency(cube);
/*
normals=(Normal *)malloc(cube.size()*sizeof(Normal));	
CalculateNormals(my_rank,cube,XYZ,ancestor_length,normals);
CheckWaterTight(my_rank,cube,XYZ,ancestor_length,normals);
free(normals);
*/

printf("level=%d\n",i);
//fflush(stdout);
}
end=clock();

printf("Run time without I/O %16.16lf\n",double(end-start)/CLOCKS_PER_SEC);

//======================================================================
#if(0)
if(my_rank==0)
{
	for(unsigned int k=0;k<cube.size();k++)
	{
//k=0;
  printf("elem_id =%d elem_level=%d\n",k,cube[k].level);

 for( int j=0;j<nface;j++)
   {
    for(unsigned int l=0;l<cube[k].nbr[j].size();l++)
    {
     printf(ANSI_COLOR_GREEN "%d \t" ANSI_COLOR_RESET,cube[k].nbr[j].at(l));
    }
    printf("\n");
   }  
   
  printf("\n===================\n");

 for(unsigned int i=0;i<cube[k].nonlocal_nbr.size();i++)
   {
     printf("************ %d nbr_id %d face tag=%d proc_id=%d\n",k,cube[k].nonlocal_nbr[i].elem_id,cube[k].nonlocal_nbr[i].face_tag,cube[k].nonlocal_nbr[i].proc_id);
   } 

printf("\n===================\n");
}
}

#endif


//======================================================================
// initialize offset since this does not get initialized for rank 0
int ncube=0,offset=0;
int mycube=cube.size();
// send and recieve should be changed to win_create to mkae it more efficient
int rma=1; offset=mycube;
// printf("mryank=%d offset =%d cube.size=%d mycube %d\n",my_rank,offset,cube.size(),mycube);

//  using active synchronization

 MPI_Group group,sub_group;
 int group_size;
 int *ranks=NULL;
 MPI_Win win;
 int new_rank;

 ranks=(int*)malloc(np*sizeof(int));
 int *off1=NULL;
 int ncube_total=0;
 calculate_processor_offset_for_IO(my_rank,mycube,1,np,&offset,ranks,off1);
 Broadcast(my_rank,mycube,0,np,offset,ranks,off1,&ncube_total);
 printf("my_rank=%d total number of cubes=%d offset=%d\n",my_rank,ncube_total,offset);
 MPI_Free_mem(off1);
 free(ranks);
 
for(unsigned int i=0;i<cube.size();i++)
{
	//cube[i].q=offset+i;
	//cube[i].q=cube[i].level;
	cube[i].q=my_rank*my_rank;
}



/*
if(cube.size()>2764)
{
	cube.at(2764).q=10000.;
}
*/
 
 /*
 if(my_rank==(np-1))
   {
     ncube_total=offset+cube.size();
   }

  MPI_Bcast(&ncube_total,1,MPI_INT,np-1,MPI_COMM_WORLD);
*/
/*
unsigned int *ia=NULL;
unsigned int *ja=NULL;
ConstructCRSFormat(my_rank,cube, &ia, &ja,MPI_COMM_WORLD,offset,ncube_total,np);
*/





/*======================================================================

                         prepare for Zoltan

  ======================================================================
*/ 
if(ZOLTAN_ON)
{
  float ver;  
  int rc;
   
 rc = Zoltan_Initialize(argcs, pArgs, &ver);

  if (rc != ZOLTAN_OK){
    printf(ANSI_COLOR_RED "Problem in Initialzing Zoltan\n" ANSI_COLOR_RESET);
    MPI_Finalize();
    exit(0);
  }
 struct Zoltan_Struct *zz=NULL;
 zz = Zoltan_Create(MPI_COMM_WORLD);
 // I define this sturcuture to mask sending all these arrays separately
 // and send only one structure instead
 
 Zoltan_Out zoltan_out;

if(ZOLTAN_GEOMETRIC_PARTITION)
{
   ZoltanGeometricPartitioner(argcs,pArgs,my_rank,np,MPI_COMM_WORLD,cube,ncube_total,offset,1,zz,XYZ,ancestor_length,&zoltan_out);
}
else 
{
   ZoltanGraphPartitioner(argcs,pArgs,my_rank,np,MPI_COMM_WORLD,cube,ncube_total,offset,0,zz,XYZ,ancestor_length,&zoltan_out);
}
/* for visualization: each processor sets the elements to be exported cube[i].q= rank_of_destination */
//if(my_rank==0)
{
	for(int i=0;i<zoltan_out.numExport;i++)
	{
		//printf("-------------%d %d\n",zoltan_out.importGlobalGids[i],zoltan_out.importToPart[i]);
		//printf("-------------%d %d\n",zoltan_out.exportLocalGids[i],zoltan_out.exportToPart[i]);
		cube[zoltan_out.exportLocalGids[i]].q=pow(zoltan_out.exportToPart[i],2);
	}	
}

//printf("my_rank-------------%d %u  %d\n",my_rank,cube.size(),cube.size()-zoltan_out.numExport+zoltan_out.numImport);
/* clean up the memory used by zoltan */

Zoltan_LB_Free_Part(&zoltan_out.importGlobalGids,&zoltan_out.importLocalGids,&zoltan_out.importProcs,&zoltan_out.importToPart);
Zoltan_LB_Free_Part(&zoltan_out.exportGlobalGids,&zoltan_out.exportLocalGids,&zoltan_out.exportProcs,&zoltan_out.exportToPart);

Zoltan_Destroy(&zz);
}
// =====================================================================
//                      prepare for ParMETIS
//======================================================================

#if(0)
//printf("***************************** %lu %d\n",ndims,sizeof(idx_t));
idx_t *part=NULL;
part=(idx_t*)malloc(cube.size()*sizeof(idx_t));
CallParMetis(cube,np,offset,ncube_total,my_rank,MPI_COMM_WORLD,part);
//if(my_rank==2)
{
	
	for(int i=0;i<cube.size();i++)
	{
	printf(ANSI_COLOR_RED "my_rank =%d i=%d id= %d\n" ANSI_COLOR_RESET,my_rank,i,part[i]);
    }
}

#endif
//======================================================================
// write routines
//======================================================================
start=clock();

if(WRITE==0)
{
 WriteHdf5ParallelPolyvertex(cube,my_rank,npx,npy,npz,MPI_COMM_WORLD,MPI_INFO_NULL,offset,ncube_total,XYZ);
}
else if(WRITE==1)
{
	 Center_coords XYZ2;
  XYZ2.reserve(4*cube.size());
WriteHdf5ParallelUnstructured(cube,my_rank,npx,npy,npz,MPI_COMM_WORLD,MPI_INFO_NULL,offset,ncube_total,XYZ, ancestor_length,XYZ2);
}
else if(WRITE==2)
{
 	 	
WriteHdf5ParallelSpatialCollection(cube,my_rank,L,M,N,npx,npy,npz,MPI_COMM_WORLD,Comm_cart,MPI_INFO_NULL,offset,ncube_total,XYZ,ancestor_length);
}
else
{
	printf(ANSI_COLOR_RED "Invalid Number for Selecting the Visualization Method\n" ANSI_COLOR_RESET);
	exit(0);
}
end=clock();

printf("Run time for I/O %16.16lf\n",double(end-start)/CLOCKS_PER_SEC);


  printf("my_rank=%d offset=%d cube_size=%lu\n",my_rank,offset,cube.size());


/*
MPI_Cart_coords(Comm_cart,my_rank,np,coords);	
   //MPI_Barrier(MPI_COMM_WORLD);
   printf("Coords %d %d %d\n",coords[0],coords[1],coords[2]);
*/
/*
for(i=0;i<6;i++)
{
	cube[0].nbr[i].push_back(-i-1);
	printf("%d\n",cube[0].nbr[i].at(0));
}
*/

double *X=NULL;
double *Y=NULL;
double *Z=NULL;



int q_size=L*M*N;



double *x=NULL;
double *y=NULL;
unsigned int n=200;


double c=1.;
int l1,l2;



int id;

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

if(my_rank==0)
{
 std::unique_ptr<int[]> ptr(new int[10]);
 
 for(int i=0;i<10;i++)
 {
	 ptr[i]=i;
	 printf("%d\n",ptr[i]);
 }
	
}

// clean up 
/*
for(i=0;i<cube.size();i++)
{
	free(cube[i].q);
}
*/




Cube cube_dealloc;

cube.swap(cube_dealloc);

Center_coords XYZ_dealloc;
XYZ.swap(XYZ_dealloc);

free(X);
free(Y);
free(Z);

//delete[] geom_xyz;
free(geom_xyz);

MPI_Finalize();

}

//======================================================================
// generates the processor connectivity  o that we can easily access without 
// extra confusion by MPI_Cart_shift and MPI_Cart_rank 

void Cartesian_Proc_Conn(const int *Dim,const int *coords,const MPI_Comm Comm_cart, Cart_proc_con *cart_proc_con)
{
  int nbr_rank;
  int coord_nbr[3];

 if(coords[0]==0)
   {
     cart_proc_con->nbr[4].push_back(-5);
   }
 else
   {
     coord_nbr[0]=coords[0]-1; 
     coord_nbr[1]=coords[1]; 
     coord_nbr[2]=coords[2]; 
    
     MPI_Cart_rank(Comm_cart,coord_nbr, &nbr_rank);
     printf("nbr_rank=%d\n",nbr_rank);
     cart_proc_con->nbr[4].push_back(nbr_rank);   
   }

 if(coords[1]==0)
   {
     cart_proc_con->nbr[2].push_back(-3);
   }
 else
   {
     coord_nbr[0]=coords[0]; 
     coord_nbr[1]=coords[1]-1; 
     coord_nbr[2]=coords[2]; 
    
     MPI_Cart_rank(Comm_cart,coord_nbr, &nbr_rank);
     cart_proc_con->nbr[2].push_back(nbr_rank);   
   }

 if(coords[2]==0)
   {
     (*cart_proc_con).nbr[0].push_back(-1);
   }
 else
   {

     coord_nbr[0]=coords[0]; 
     coord_nbr[1]=coords[1]; 
     coord_nbr[2]=coords[2]-1; 
    
     MPI_Cart_rank(Comm_cart,coord_nbr, &nbr_rank);
     cart_proc_con->nbr[0].push_back(nbr_rank);   
   }

 // whatch out the Dim starts from 1 not zero, in contrast to coords

 if(coords[0]==(Dim[0]-1))
   {
     cart_proc_con->nbr[5].push_back(-6);
   }
 else
   {

     coord_nbr[0]=coords[0]+1; 
     coord_nbr[1]=coords[1]; 
     coord_nbr[2]=coords[2]; 
    
     MPI_Cart_rank(Comm_cart,coord_nbr, &nbr_rank);

     cart_proc_con->nbr[5].push_back(nbr_rank);   
   }

 if(coords[1]==(Dim[1]-1))
   {
     cart_proc_con->nbr[3].push_back(-4);
   }
 else
   {

     coord_nbr[0]=coords[0]; 
     coord_nbr[1]=coords[1]+1; 
     coord_nbr[2]=coords[2]; 
     MPI_Cart_rank(Comm_cart,coord_nbr, &nbr_rank);

     cart_proc_con->nbr[3].push_back(nbr_rank);  

 
   }

 if(coords[2]==(Dim[2]-1))
   {
     cart_proc_con->nbr[1].push_back(-2);
   }
 else
   {
     coord_nbr[0]=coords[0]; 
     coord_nbr[1]=coords[1]; 
     coord_nbr[2]=coords[2]+1; 
    
     MPI_Cart_rank(Comm_cart,coord_nbr, &nbr_rank);

     cart_proc_con->nbr[1].push_back(nbr_rank);   
   }
}

//===============================================================

void find_boundary_face(Cube &cube,int cube_id,Vector_Int &bface_id)
 {  

   unsigned int i,j;
 
     for(i=0;i<6;i++)
       {
	 for(j=0;j<cube[cube_id].nbr[i].size();j++)
	   {
	     if(cube[cube_id].nbr[i].at(j)<0)
	       {
		 bface_id.push_back(i);
	       }
	   }	 
       } 
 }

//======================================================================
