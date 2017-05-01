#include <stdlib.h>
//#include "/home/jhasbestan/GEM_AMR/src/include/typedefs.h"
#include "typedefs.h"
#define TOL "1.1"

/*
   
   Two Different Data Structures are Used for Graph and Mesh Data
         Only hypergraph can handle nonsymmetric graphs    
 
 ----------------------------------------------------------------------
                  GEOMETRIC (Coordinate-Based) of Zoltan
  ----------------------------------------------------------------------
   
   set method=1 for HSFC: Hilbert Space Filling Curve 
   set method=2 for  RCB: Recursive Coordinate Bisection 
   set method=3 for  RIB: Recursive Inertial Bisection 
   
    need to use coordinate transformation to impose the principal peorcessor topology
    direction, otherwise it this routine starts from z=0 plane, this 
    might not necessarily be the best option as it might require too many
    elements to be moved unnecessarily, just rotate the xyz directions if you
    have more processor in x and y directions. e.g. 3 1 1 is needs transformation
    but 1 1 3 is ok, default is that nz>=ny and nz >=nx  
 -----------------------------------------------------------------------
                       HyperGraph of Zoltan    
 -----------------------------------------------------------------------
   
   For HyperGraph in Zoltan
   set method=0
 -----------------------------------------------------------------------
                            ParMETIS
 -----------------------------------------------------------------------
   
   set method=1 for PARTKWAY
   set method=2 for PARTGEOMKWAY (Hybrid Method uses geom+graph)
   set method=3 for ADEAPTIVEREPART
   set method=4 for PARTGEOM
   the *c in GRAPH_DATA is only assigned when needed 
   Default method is ADAPTIVEREPART as set per Zoltan
    
   #define TOL "1.1" defines the imbalance tolerance for 
   Graph and Hypergraph  methods
   
 -----------------------------------------------------------------------
*/ 
 

typedef struct{
  ZOLTAN_ID_TYPE numGlobalPoints;
  ZOLTAN_ID_TYPE numMyPoints;
  ZOLTAN_ID_PTR myGlobalIDs;
  double *c;
} MESH_DATA;

typedef struct{
  ZOLTAN_ID_TYPE numMyVertices; /* total vertices in in my partition */
 // ZOLTAN_ID_TYPE numAllNbors;   /* total number of neighbors of my vertices */
  ZOLTAN_ID_PTR vertexGID;    /* global ID of each of my vertices */
  ZOLTAN_ID_PTR nbrIndex;    /* nborIndex[i] is location of start of neighbors for vertex i */
  ZOLTAN_ID_PTR nbrGID;      /* nborGIDs[nborIndex[i]] is first neighbor of vertex i */
  int *nbrProc;     /* process owning each nbor in nborGID */
  float *c; /* this one is for the center of the elements*/
} GRAPH_DATA;


// refer to manual for ZOLTAN_ID_TYPE and ZOLTAN_ID_PTR
// note that ZOLTAN_ID_PTR=*ZOLTAN_ID_TYPE

static int get_number_of_objects(void *data, int *ierr);

static int get_num_geometry(void *data, int *ierr);

static void get_object_list(void *data, int sizeGID, int sizeLID,
            ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                  int wgt_dim, float *obj_wgts, int *ierr);
                  
static void get_geometry_list(void *data, int sizeGID, int sizeLID,
                      int num_obj,
             ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
             int num_dim, double *geom_vec, int *ierr);
             
           
void ZoltanGeometricPartitioner(int argcs, char* pArgs[],int my_rank,int np,MPI_Comm Comm,Cube &cube,int ncube_total,int offset,int method,struct Zoltan_Struct *zz,Center_coords &XYZ, double ancestor_length[3],Zoltan_Out *zoltan_out)
{

  float ver;  
  int rc;
  
  /* General parameters */

#if(1)

  Zoltan_Set_Param(zz, "DEBUG_LEVEL", "1");

if(method==1)
{
 Zoltan_Set_Param(zz, "LB_METHOD", "HSFC");
}
else if(method==2)
{
    Zoltan_Set_Param(zz, "LB_METHOD", "RCB");
}
else if(method==3)
{
    Zoltan_Set_Param(zz, "LB_METHOD", "RIB");
}
else
{
printf("Specified Method Needs a Graph\n");
printf(ANSI_COLOR_RED "ERROR:Zoltan error \n" ANSI_COLOR_RESET);
 MPI_Finalize();
exit(0);
}
  //  Zoltan_Set_Param(zz, "LB_METHOD", "BLOCK");
  Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1"); 
  Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");
  //Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "0");
  Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");


	MESH_DATA myMesh;
//	myMesh=(MESH_DATA*)malloc(1.0*sizeof(MESH_DATA));
  	
  	myMesh.numGlobalPoints=ncube_total;
   	myMesh.numMyPoints=cube.size();
	myMesh.myGlobalIDs=NULL;
	myMesh.myGlobalIDs=(ZOLTAN_ID_PTR)malloc(cube.size()*sizeof(ZOLTAN_ID_TYPE));
	
		
	for(unsigned int i=0;i<cube.size();i++)
	{
	  myMesh.myGlobalIDs[i]=offset+i;
	}
   
   myMesh.c=NULL;
   myMesh.c=(double*)malloc(3*myMesh.numMyPoints*sizeof(double));
   
   // due to modification in structure
   
        
   
   
   for(unsigned int i=0;i<myMesh.numMyPoints;i++)
   {
	/*
	myMesh.c[3*i]=0.5*(cube[i].xyz[0]+cube[i].xyz[1]); 
	myMesh.c[3*i+1]=0.5*(cube[i].xyz[2]+cube[i].xyz[3]); 
	myMesh.c[3*i+2]=0.5*(cube[i].xyz[4]+cube[i].xyz[5]); 
	*/
    myMesh.c[3*i]=XYZ.at(i).x; 
	myMesh.c[3*i+1]=XYZ.at(i).y; 
	myMesh.c[3*i+2]=XYZ.at(i).z; 

	 
	//printf("my_rank %d %lf\n",my_rank,myMesh.c[i]);  
   }
   
   // printf("my_rank =%d number of my points %d glob_id global_id %d %d \n",my_rank,myMesh->numMyPoints,myMesh->myGlobalIDs[0],myMesh->myGlobalIDs[1]);
  
  
  /* Query functions, to provide geometry to Zoltan */

    Zoltan_Set_Num_Obj_Fn(zz, get_number_of_objects, &myMesh);

    Zoltan_Set_Obj_List_Fn(zz, get_object_list, &myMesh);
    
    Zoltan_Set_Num_Geom_Fn(zz, get_num_geometry, &myMesh);
    
    Zoltan_Set_Geom_Multi_Fn(zz, get_geometry_list, &myMesh);
    
    
  
  
  	rc = Zoltan_LB_Partition(zz, /* input (all remaining fields are output) */
        &(zoltan_out->changes),        /* 1 if partitioning was changed, 0 otherwise */ 
        &(zoltan_out->numGidEntries),  /* Number of integers used for a global ID */
        &(zoltan_out->numLidEntries),  /* Number of integers used for a local ID */
        &(zoltan_out->numImport),      /* Number of vertices to be sent to me */
        &(zoltan_out->importGlobalGids),  /* Global IDs of vertices to be sent to me */
        &(zoltan_out->importLocalGids),   /* Local IDs of vertices to be sent to me */
        &(zoltan_out->importProcs),    /* Process rank for source of each incoming vertex */
        &(zoltan_out->importToPart),   /* New partition for each incoming vertex */
        &(zoltan_out->numExport),      /* Number of vertices I must send to other processes*/
        &(zoltan_out->exportGlobalGids),  /* Global IDs of the vertices I must send */
        &(zoltan_out->exportLocalGids),   /* Local IDs of the vertices I must send */
        &(zoltan_out->exportProcs),    /* Process to which I send each of the vertices */
        &(zoltan_out->exportToPart));  /* Partition to which each vertex will belong */

#endif
if (rc != ZOLTAN_OK){
    printf(ANSI_COLOR_RED "Zoltan Partition Returned Error...\n" ANSI_COLOR_RESET);
    MPI_Finalize();
    exit(0);
  }

//printf( "my_rank %d changes %d numExport %d numImport %d\n",my_rank,zoltan_out->changes,zoltan_out->numExport,zoltan_out->numImport);
 
 
}
// for passing void and type conversion inside a function and how it works see The C++ Programming Language by Bjarne S. section 7.2
static int get_number_of_objects(void *data, int *ierr)
{
  MESH_DATA *mesh= (MESH_DATA *)data;
  *ierr = ZOLTAN_OK;
  //printf("number of my points %d\n",mesh->numMyPoints);
  return mesh->numMyPoints;
  
}

static void get_object_list(void *data, int sizeGID, int sizeLID,
            ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                  int wgt_dim, float *obj_wgts, int *ierr)
{

  MESH_DATA *mesh= (MESH_DATA *)data;
  *ierr = ZOLTAN_OK;

 //printf("?????????????????\n \n %d \n",mesh->numMyPoints);

   for (unsigned int i=0; i<mesh->numMyPoints; i++)
   {
    globalID[i] = mesh->myGlobalIDs[i];
    localID[i] = i;
    //printf("local global %d %d\n",localID[i] ,globalID[i]);
  }
  
  
}

//======================================================================

static int get_num_geometry(void *data, int *ierr)
{
  *ierr = ZOLTAN_OK;
  return 3;
}

//======================================================================

static void get_geometry_list(void *data, int sizeGID, int sizeLID,
                      int num_obj,
             ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
             int num_dim, double *geom_vec, int *ierr)
{
	
int i;

  MESH_DATA *mesh= (MESH_DATA *)data;

/*
  if ( (sizeGID != 1) || (sizeLID != 1) || (num_dim != 2)){
    *ierr = ZOLTAN_FATAL;
    return;
  }
*/

  *ierr = ZOLTAN_OK;

  for (i=0;  i < num_obj ; i++){
    geom_vec[3*i] = (double)mesh->c[3*i];
    geom_vec[3*i + 1] =(double)mesh->c[3*i+1];
    geom_vec[3*i + 2] = (double)mesh->c[3*i+2];
  }

  
}


//======================================================================
//
//                        Graph Partitioning
//
//======================================================================



static int get_number_of_vertices(void *data, int *ierr);

static void get_vertex_list(void *data, int sizeGID, int sizeLID,
            ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                  int wgt_dim, float *obj_wgts, int *ierr);
                  
static void get_num_edges_list(void *data, int sizeGID, int sizeLID,
                      int num_obj,
             ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
             int *numEdges, int *ierr);

static void get_edge_list(void *data, int sizeGID, int sizeLID,
        int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
        int *num_edges,
        ZOLTAN_ID_PTR nbrGID, int *nbrProc,
        int wgt_dim, float *ewgts, int *ierr);


void GetOffSetToGloballyRename(MPI_Comm Comm,int np,const Cube &cube,GRAPH_DATA *myGraph,int ncube_total,int offset);


static void get_geometry_list2(void *data, int sizeGID, int sizeLID,
                      int num_obj,
             ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
             int num_dim, double *geom_vec, int *ierr);   
            
void ZoltanGraphPartitioner(int argcs, char* pArgs[],int my_rank,int np,MPI_Comm Comm,const Cube &cube,int ncube_total,int offset,int method,struct Zoltan_Struct *zz,Center_coords &XYZ, double ancestor_length[3],Zoltan_Out *zoltan_out)
{

  
  int rc;
  
  
  /* General parameters */

  Zoltan_Set_Param(zz, "IMBALANCE_TOL", TOL);
  Zoltan_Set_Param(zz, "DEBUG_LEVEL", "1");
  Zoltan_Set_Param(zz, "LB_METHOD", "GRAPH");
  Zoltan_Set_Param(zz, "LB_APPROACH", "REPARTITION");
  
  //Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");
  

  
 if(method!=0)
 { 
  Zoltan_Set_Param(zz,"GRAPH_PACKAGE","PARMETIS");
   if(method==1)
  {
    Zoltan_Set_Param(zz, "PARMETIS_METHOD", "PARTKWAY");
  }
  else if(method==2)
  {
  Zoltan_Set_Param(zz, "PARMETIS_METHOD", "PARTGEOMKWAY");
  }
  else if(method==3)
  {
  Zoltan_Set_Param(zz, "PARMETIS_METHOD", "ADAPTIVEREPART");
   }
   else if(method==4)
   {
	   Zoltan_Set_Param(zz, "PARMETIS_METHOD", "PARTGEOM");
   }
   else
   {
	   
	   printf(ANSI_COLOR_RED "Undefined Method for GRAPH_PARTITIONING \n" ANSI_COLOR_RESET);
	   MPI_Finalize();
	   exit(0);
	    }
	}
  // Zoltan_Set_Param(zz, "PARMETIS_METHOD", "REFINEKWAY");  
  
 
   // just to debug and show that PARTGEOM produces wrong result
   // need to uncomment line  if(method%2==0) along with this
   //      Zoltan_Set_Param(zz, "LB_METHOD", "HSFC");
  //Zoltan_Set_Param(zz, "PARMETIS_METHOD", "PARTGEOM");
  
  Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1"); 
  Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");
  Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");

	GRAPH_DATA myGraph;

  	myGraph.numMyVertices=cube.size();	
   	//myGraph.vertexGID=NULL;
   	//myGraph.vertexGID=(ZOLTAN_ID_PTR)malloc(myGraph.numMyVertices*sizeof(ZOLTAN_ID_TYPE));
   	
   	myGraph.vertexGID=new ZOLTAN_ID_TYPE[myGraph.numMyVertices];

	for(unsigned int i=0;i<cube.size();i++)
	{
	  myGraph.vertexGID[i]=offset+i;
	}
	
	// fill out the data for the Graph	   
   // printf("my_rank =%d number of my points %d glob_id global_id %d %d \n",my_rank,myMesh->numMyPoints,myMesh->myGlobalIDs[0],myMesh->myGlobalIDs[1]);
  
    /* Query functions, to provide geometry to Zoltan */
    
    Zoltan_Set_Num_Geom_Fn(zz, get_num_geometry, &myGraph);
  // this is to assign geometirc data required by methods 2 and 4
  // the methods 1 and 3 do not need any geometirc data 
    if(method%2==0)
    {
    myGraph.c=NULL;
   myGraph.c=(float*)malloc(3*myGraph.numMyVertices*sizeof(float));
   
   for(unsigned int i=0;i<myGraph.numMyVertices;i++)
   {
	/*
	myGraph.c[3*i]=0.5*(cube[i].xyz[0]+cube[i].xyz[1]); 
	myGraph.c[3*i+1]=0.5*(cube[i].xyz[2]+cube[i].xyz[3]); 
	myGraph.c[3*i+2]=0.5*(cube[i].xyz[4]+cube[i].xyz[5]); 
	//printf("my_rank %d %lf\n",my_rank,myGraph.c[i]);  
	*/
	myGraph.c[3*i]=XYZ.at(i).x; 
	myGraph.c[3*i+1]=XYZ.at(i).y; 
	myGraph.c[3*i+2]=XYZ.at(i).z; 
	//printf("my_rank %d %lf\n",my_rank,myGraph.c[i]);
	
   }
   
    Zoltan_Set_Num_Geom_Fn(zz, get_num_geometry, &myGraph);
    Zoltan_Set_Geom_Multi_Fn(zz, get_geometry_list2, &myGraph);
  }
   // printf("my_rank =%d number of my points %d glob_id global_id %d %d \n",my_rank,myMesh->numMyPoints,myMesh->myGlobalIDs[0],myMesh->myGlobalIDs[1]);
    
  /* Query functions, to provide geometry to Zoltan */
   
    Zoltan_Set_Num_Obj_Fn(zz, get_number_of_vertices, &myGraph);
    Zoltan_Set_Obj_List_Fn(zz, get_vertex_list, &myGraph);
    
    //printf("my_rank=%d myGraph.numMyVertices=%d\n",my_rank,myGraph.numMyVertices);
    
       
    GetOffSetToGloballyRename(MPI_COMM_WORLD,np,cube,&myGraph,ncube_total,offset);
      
    Zoltan_Set_Num_Edges_Multi_Fn(zz, get_num_edges_list, &myGraph);
    Zoltan_Set_Edge_List_Multi_Fn(zz,get_edge_list, &myGraph);
    
   
    /*
    if(my_rank==0)
    {
    for(int i=0;i<myGraph.numMyVertices;i++)
     {
		for(int j=myGraph.nbrIndex[i];j<myGraph.nbrIndex[i+1];j++)
		{
			printf("%d\t",myGraph.nbrGID[j]);						
		}
		printf("\n");
	 }
    }
   */
#if(1) 
 
    //Zoltan_Out zoltan_out;
    
  
  	rc = Zoltan_LB_Partition(zz, /* input (all remaining fields are output) */
        &(zoltan_out->changes),        /* 1 if partitioning was changed, 0 otherwise */ 
        &(zoltan_out->numGidEntries),  /* Number of integers used for a global ID */
        &(zoltan_out->numLidEntries),  /* Number of integers used for a local ID */
        &(zoltan_out->numImport),      /* Number of vertices to be sent to me */
        &(zoltan_out->importGlobalGids),  /* Global IDs of vertices to be sent to me */
        &(zoltan_out->importLocalGids),   /* Local IDs of vertices to be sent to me */
        &(zoltan_out->importProcs),    /* Process rank for source of each incoming vertex */
        &(zoltan_out->importToPart),   /* New partition for each incoming vertex */
        &(zoltan_out->numExport),      /* Number of vertices I must send to other processes*/
        &(zoltan_out->exportGlobalGids),  /* Global IDs of the vertices I must send */
        &(zoltan_out->exportLocalGids),   /* Local IDs of the vertices I must send */
        &(zoltan_out->exportProcs),    /* Process to which I send each of the vertices */
        &(zoltan_out->exportToPart));  /* Partition to which each vertex will belong */


if (rc != ZOLTAN_OK){
    printf(ANSI_COLOR_RED "Zoltan Partition Returned Error...\n" ANSI_COLOR_RESET);
    MPI_Finalize();
    exit(0);
  }

//printf( "my_rank %d changes %d numExport %d numImport %d\n",my_rank,zoltan_out->changes,zoltan_out->numExport,zoltan_out->numImport);
                    
#endif

 // clear up memory in use by myGraph

delete[] myGraph.vertexGID;
delete[] myGraph.nbrGID;
delete[] myGraph.nbrIndex;
delete[] myGraph.nbrProc;

 

if(method%2==0)
    {
   
   delete[] myGraph.c;
}

}

//======================================================================

static int get_number_of_vertices(void *data, int *ierr)
{
  GRAPH_DATA *graph = (GRAPH_DATA *)data;
  *ierr = ZOLTAN_OK;
  return graph->numMyVertices;
}

//======================================================================

static void get_vertex_list(void *data, int sizeGID, int sizeLID,
            ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                  int wgt_dim, float *obj_wgts, int *ierr)
{


  GRAPH_DATA *graph = (GRAPH_DATA *)data;
  *ierr = ZOLTAN_OK;

  /* In this example, return the IDs of our vertices, but no weights.
   * Zoltan will assume equally weighted vertices.
   */

  for (unsigned int i=0; i<graph->numMyVertices; i++){
    globalID[i] = graph->vertexGID[i];
    localID[i] = i;
  }
}

//======================================================================
void finproc(int np,unsigned int *vtx_slist,unsigned int nbrGId,int *proc_id);
void GetOffSetToGloballyRename(MPI_Comm Comm,int np,const Cube &cube,GRAPH_DATA *myGraph,int ncube_total,int offset)
{

unsigned int *vtx_slist=NULL;
unsigned int *vtx_rlist=NULL;

vtx_slist=new unsigned int[np+1];
vtx_rlist=new unsigned int[np+1];


for(int i=0;i<np;i++)
{
vtx_slist[i]=offset;
}

MPI_Alltoall(vtx_slist, 1, MPI_UNSIGNED, vtx_rlist, 1, MPI_UNSIGNED, Comm);

vtx_rlist[np]=ncube_total;

  //  myGraph->nbrIndex=NULL;
	//myGraph->nbrIndex=(ZOLTAN_ID_PTR)malloc((cube.size()+1)*sizeof(ZOLTAN_ID_TYPE));
	myGraph->nbrIndex=new ZOLTAN_ID_TYPE[(cube.size()+1)];
	
unsigned int counter=0;
myGraph->nbrIndex[0]=0;

for(unsigned int i=0;i<cube.size();i++)
{
	counter=0;

	for(int j=0;j<6;j++)
	{
		if(cube[i].nbr[j].at(0)>-1)
		{
			counter=counter+cube[i].nbr[j].size();
		}		
	}
	
	counter=counter+cube[i].nonlocal_nbr.size();
	myGraph->nbrIndex[i+1]=myGraph->nbrIndex[i]+counter;
}

	//myGraph->nbrGID=NULL;
	//myGraph->nbrGID=(ZOLTAN_ID_PTR)malloc(myGraph->nbrIndex[cube.size()]*sizeof(ZOLTAN_ID_TYPE));
	myGraph->nbrGID=new ZOLTAN_ID_TYPE[myGraph->nbrIndex[cube.size()]];
	myGraph->nbrProc=new int[myGraph->nbrIndex[cube.size()]];
	//myGraph->nbrProc=NULL;
	//myGraph->nbrProc=(int*)malloc(myGraph->nbrIndex[cube.size()]*sizeof(int));
	
	
	counter=0;
	int elem_id,proc_id;
	
for(unsigned int i=0;i<cube.size();i++)
{	
	for(int j=0;j<6;j++)
	{
		if(cube[i].nbr[j].at(0)>-1)
		{					
			for(unsigned int k=0;k<cube[i].nbr[j].size();k++)
			{
				myGraph->nbrGID[counter]=cube[i].nbr[j].at(k)+offset;
				counter++;
			}
		}		
	}
	
	for(unsigned int k=0;k<cube[i].nonlocal_nbr.size();k++)
	{
		proc_id=cube[i].nonlocal_nbr.at(k).proc_id;
		elem_id=cube[i].nonlocal_nbr.at(k).elem_id;
		myGraph->nbrGID[counter]=elem_id+vtx_rlist[proc_id];
		counter++;
    }
}

   for(unsigned int i=0;i<counter;i++)
   {
	   elem_id=myGraph->nbrGID[i];
	   finproc(np,vtx_rlist,elem_id,&proc_id);
	   myGraph->nbrProc[i]=proc_id;	   
   }

delete[] vtx_rlist;
delete[] vtx_slist;
		
}

void finproc(int np,unsigned int *vtx_rlist,unsigned int nbrGId,int *proc_id)
{
	*proc_id=0;
	
	for(int i=0;i<np;i++)
	{
		if(nbrGId>=vtx_rlist[i+1])
		{
			*proc_id=*proc_id+1;
			//printf("nbrGId=%d pro_id = %d\n",nbrGId,*proc_id);
		}
	}
		
}

//======================================================================

static void get_num_edges_list(void *data, int sizeGID, int sizeLID,
                      int num_obj,
             ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
             int *numEdges, int *ierr)
{
 int i, idx;

  GRAPH_DATA *graph = (GRAPH_DATA *)data;
/*
  if ( (sizeGID != 1) || (sizeLID != 1) || (num_obj != graph->numMyVertices)){
    *ierr = ZOLTAN_FATAL;
    return;
  }
*/
  for (i=0;  i < num_obj ; i++){
    idx = localID[i];
    numEdges[i] = graph->nbrIndex[idx+1] - graph->nbrIndex[idx];
  }

  *ierr = ZOLTAN_OK;
  return;
}


//======================================================================

static void get_edge_list(void *data, int sizeGID, int sizeLID,
        int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
        int *num_edges,
        ZOLTAN_ID_PTR nbrGID, int *nbrProc,
        int wgt_dim, float *ewgts, int *ierr)
{
	
int i, j, from, to;
int *nextProc;
ZOLTAN_ID_TYPE *nextNbr;

  GRAPH_DATA *graph = (GRAPH_DATA *)data;
  *ierr = ZOLTAN_OK;
/*
  if ( (sizeGID != 1) || (sizeLID != 1) || 
       (num_obj != graph->numMyVertices)||
       (wgt_dim != 0)){
    *ierr = ZOLTAN_FATAL;
    return;
  }
*/
  nextNbr = nbrGID;
  nextProc = nbrProc;

  for (i=0; i < num_obj; i++)
  {
    to = graph->nbrIndex[localID[i]+1];
    from = graph->nbrIndex[localID[i]];
    if ((to - from) != num_edges[i]){
      *ierr = ZOLTAN_FATAL;
      return;
    }

    for (j=from; j < to; j++){

      *nextNbr++ = graph->nbrGID[j];
      *nextProc++ = graph->nbrProc[j];
    }
  }
  return;
}




static void get_geometry_list2(void *data, int sizeGID, int sizeLID,
                      int num_obj,
             ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
             int num_dim, double *geom_vec, int *ierr)
{
	
int i;

  GRAPH_DATA *mesh= (GRAPH_DATA *)data;


  if ( (sizeGID != 1) || (sizeLID != 1) || (num_dim != 3)){
    *ierr = ZOLTAN_FATAL;
    return;
  }


  *ierr = ZOLTAN_OK;

  for (i=0;  i < num_obj ; i++){
    geom_vec[3*i] = (double)mesh->c[3*i];
    geom_vec[3*i + 1] =(double)mesh->c[3*i+1];
    geom_vec[3*i + 2] = (double)mesh->c[3*i+2];
    //printf("Geom=%lf %lf %lf\n",geom_vec[3*i ],geom_vec[3*i + 1],geom_vec[3*i + 2]);
  }

  
}




