//#include "/home/jhasbestan/GEM_AMR/src/include/typedefs.h"
#include "typedefs.h"

void Save(const int id,const Cube& cube,Cube& cube_temp)
{
  cube_temp[0].level=cube[id].level;
	
  cube_temp[0].centeroid_index=cube[id].centeroid_index;
     
  for(unsigned int i=0;i<6;i++)
    {
      for(unsigned int j=0;j<cube[id].nbr[i].size();j++)
	{
	  cube_temp[0].nbr[i].push_back(cube[id].nbr[i].at(j));
	}
      //printf("%d\n",cube[0].nbr[i].at(0));
    }

// for now all I need to know is that the element is a boundary element
// in this program I assume that an element is a boundary element if
// the  cube[id].nonlocal_nbr.size()!=0
// inherit only the face_tag and proc_id
// these elements will be filled out after communication

/* save the newly added tuplet to the structure needed for parallel implementation
   this is in order for the boundary element to know
   on which processor to which element it is connected to
 */

	
	for(unsigned int i=0;i<cube[id].nonlocal_nbr.size();i++)	
	{
	  cube_temp[0].nonlocal_nbr.push_back(Nonlocal_Nbr());
	  cube_temp[0].nonlocal_nbr[i].elem_id=cube[id].nonlocal_nbr[i].elem_id;
	  cube_temp[0].nonlocal_nbr[i].proc_id=cube[id].nonlocal_nbr[i].proc_id;
	  cube_temp[0].nonlocal_nbr[i].face_tag=cube[id].nonlocal_nbr[i].face_tag;
    }
	
	
}

//======================================================================
void EmptyTempCube(Cube& cube_temp)
{
	for(unsigned int i=0;i<6;i++)
	{
	while(cube_temp[0].nbr[i].size()!=0)
	{
	cube_temp[0].nbr[i].pop_back();
     }
    }
    
  while(cube_temp[0].nonlocal_nbr.size()!=0)
	{
	cube_temp[0].nonlocal_nbr.pop_back();
     }
	
	
}
//======================================================================
/*
 * 
 *  Assign coordinates of children according to the parent coordinates
 * 
 * 
 */
//======================================================================
 
void ConstructChildrenCoords(const int id,const Cube &cube_temp,const unsigned int *map_index, Cube &cube,Center_coords& XYZ,double ancestor_length[3])
{

 unsigned int level,coord_index;
 double idenum;
 double dx,dy,dz;
 
 level=cube_temp.at(0).level;
 coord_index=cube_temp.at(0).centeroid_index;
 
 //idenum=1./pow(2,level);
 
 TwoPowN(level,&idenum);
 idenum=1./idenum;
	  
 dx=ancestor_length[0]*idenum;
 dy=ancestor_length[1]*idenum;
 dz=ancestor_length[2]*idenum;
	  
 double xmin=XYZ.at(coord_index).x-dx*0.5; 	
 double ymin=XYZ.at(coord_index).y-dy*0.5;
 double zmin=XYZ.at(coord_index).z-dz*0.5;

 double xmax=XYZ.at(coord_index).x+dx*0.5;
 double ymax=XYZ.at(coord_index).y+dy*0.5;
 double zmax=XYZ.at(coord_index).z+dz*0.5;  

 double xmid=(xmin+xmax)*0.5;
 double ymid=(ymin+ymax)*0.5;
 double zmid=(zmin+zmax)*0.5;
  
  
double X[3]={xmin,xmid,xmax};
double Y[3]={ymin,ymid,ymax};
double Z[3]={zmin,zmid,zmax};

static unsigned int Xmin_idx[8]={0,1,1,0,0,1,1,0};
static unsigned int Xmax_idx[8]={1,2,2,1,1,2,2,1};

static unsigned int Ymin_idx[8]={0,0,1,1,0,0,1,1};
static unsigned int Ymax_idx[8]={1,1,2,2,1,1,2,2};

static unsigned int Zmin_idx[8]={0,0,0,0,1,1,1,1};
static unsigned int Zmax_idx[8]={1,1,1,1,2,2,2,2};

//printf("xyz min max=%lf %lf %lf %lf %lf %lf\n",xmin,xmax,ymin,ymax,zmin,zmax);

for(unsigned i=0;i<8;i++)
{
 cube.at(map_index[i]).centeroid_index=map_index[i];
 
 XYZ.at(map_index[i]).x=(X[Xmin_idx[i]]+X[Xmax_idx[i]])*0.5;
 XYZ.at(map_index[i]).y=(Y[Ymin_idx[i]]+Y[Ymax_idx[i]])*0.5;
 XYZ.at(map_index[i]).z=(Z[Zmin_idx[i]]+Z[Zmax_idx[i]])*0.5;
 
// printf("%lf %lf %lf\n",XYZ.at(map_index[i]).x,XYZ.at(map_index[i]).y,XYZ.at(map_index[i]).z);
 
}

}
//======================================================================
/*
 * 
 *           Remove the parent element from the connectivity 
 * 
 */
//======================================================================
 
void RemoveParentFromConn(const int id,const Cube &cube_temp,const unsigned int* co_face_no,Cube &cube)
{
 for(unsigned int i=0;i<6;i++) 
{
   for(unsigned int j=0;j<cube_temp[0].nbr[i].size();j++)
     {
     if(cube_temp[0].nbr[i].at(j)>-1) 
       {
       for(unsigned int k=0;k<cube[cube_temp[0].nbr[i].at(j)].nbr[co_face_no[i]].size();k++)
	 {
	 if(cube[cube_temp[0].nbr[i].at(j)].nbr[co_face_no[i]].at(k)==id)
	   {		  
	     cube[cube_temp[0].nbr[i].at(j)].nbr[co_face_no[i]].erase(cube[cube_temp[0].nbr[i].at(j)].nbr[co_face_no[i]].begin()+k);
	     k=k-1;
	   }
	 }
       }     
     }
   }
}

//======================================================================
// Assign interior faces for each element, no change in here for parallel
// implementation 
//======================================================================

void AssignInteriorElemConn(const unsigned int id,
	const unsigned int nquad,const unsigned int *map_index,
	const unsigned int *I,Cube &cube)
{

unsigned int nbr_index[24]={nquad+3,nquad+2,nquad,nquad+4,nquad+1,id,nquad+5,nquad,nquad+2,nquad+6,id,nquad+1,id,nquad+6,nquad+4,nquad,nquad+5,nquad+3,nquad+1,nquad+4,nquad+6,nquad+2,nquad+3,nquad+5};

for(unsigned int i=0;i<8;i++)
{ 
	for(unsigned int j=0;j<3;j++)
	{
		cube[map_index[i]].nbr[I[3*i+j]].push_back(nbr_index[3*i+j]);
		//printf(" %d %d\t",I[3*i+j],nbr_index[3*i+j]);
	}	
}

	
}

//======================================================================
/*
 *
 *   extrior face assignemnt using mapping eliminates search and one 
 *   communication
 *  
 */ 
//=====================================================================
//=====================================================================
 // since the nonlocal_nbr stores all all nbrs in one container first we need to 
 // look for the nbrs that have the same face_tag
 
 void CollectSameTag(int my_rank,const Cube &cube_temp,const unsigned int face_id, Vector_Int &temp_elem,Vector_Int &temp_proc)
 {
	  for(unsigned int k=0;k<cube_temp[0].nonlocal_nbr.size();k++)
	        {
				if(cube_temp[0].nonlocal_nbr.at(k).face_tag==face_id)
				{
					/*
					if(my_rank==1)
					{
					printf("inside %d\n",face_id);
				    }
				    */ 
					temp_elem.push_back(cube_temp[0].nonlocal_nbr.at(k).elem_id);
					temp_proc.push_back(cube_temp[0].nonlocal_nbr.at(k).proc_id);
				}
			}
}

//======================================================================
void EmptyContainer(Vector_Int &V)
{
  
	V.clear();
//	V.shrink_to_fit();
	
}
 
void AssignThisFace(int my_rank,Cube &cube,const unsigned int id,const unsigned face_id,const unsigned int i,const unsigned int j,const Vector_Int &temp_elem, const Vector_Int &temp_proc)
{
	
	unsigned int idx,prc,size;
 // these indices are meticulously chosen to eliminate search
 static unsigned int K[24]={0,1,0,1,0,0,2,1,1,3,0,1,0,3,2,1,2,2,2,3,3,3,2,3};
 
 
	if(temp_elem.size()==1)
	        {

		   size=cube[id].nonlocal_nbr.size();	
	       cube[id].nonlocal_nbr.push_back(Nonlocal_Nbr());
	       cube[id].nonlocal_nbr[size].face_tag=face_id;
	       cube[id].nonlocal_nbr[size].proc_id=temp_proc[0];
		   // assigne for now, if nbr is modified will overwrite this after communication
		   //printf(ANSI_COLOR_GREEN "here\n" ANSI_COLOR_RESET);
		   cube[id].nonlocal_nbr[size].elem_id=temp_elem[0];
		/*
		if(my_rank==1)
		{
			printf("----------------------------%d %d %u\n",face_id,temp_elem[0],cube[id].nonlocal_nbr[0].face_tag);		 
		}
		*/
			}
			
	else if(temp_elem.size()==4)
			{
				// divide the 4 elem nbr's between the newly generated elements
			   	idx=temp_elem.at(K[3*i+j]);
			   	prc=temp_proc.at(K[3*i+j]);  	
			   	size=cube[id].nonlocal_nbr.size();
			   	cube[id].nonlocal_nbr.push_back(Nonlocal_Nbr());
			   	cube[id].nonlocal_nbr[size].face_tag=face_id;
			   	cube[id].nonlocal_nbr[size].proc_id=prc;
			   	cube[id].nonlocal_nbr[size].elem_id=idx;
			   	
			   	//printf(ANSI_COLOR_BLUE "here\n" ANSI_COLOR_RESET);
			   	// updating the nbr on the other processor is due after comminucation
			}
	else
			{
				printf("------------------------------------------%lu\n",temp_elem.size());
				printf(ANSI_COLOR_RED "element has more than 4 neighbors (in processor nbr)\n" ANSI_COLOR_RESET);
			    exit(0);
			}
			
				
}


void AssignExteriorFace(int my_rank,Cube &cube,const Cube &cube_temp,const unsigned int *co_face_no,const unsigned int *map_index,const unsigned int*I) 
 {
 
  int idx;
 unsigned int id;
 unsigned int face_id;
 // these indices are meticulously chosen to eliminate search
 static unsigned int K[24]={0,1,0,1,0,0,2,1,1,3,0,1,0,3,2,1,2,2,2,3,3,3,2,3};		 
 Vector_Int temp_proc,temp_elem;		 
 //printf("coordinates current element =%lf %lf \n",xc_id,yc_id);

for(unsigned int i=0;i<8;i++)
{
	//unsigned int i=1;
	id=map_index[i];
	
	
	for(unsigned int j=0;j<3;j++)
	{	
	//unsigned int j=0;
	 face_id=I[3*i+j];		 

    //printf("nonlocal_nbr.size=%d\n",cube_temp[0].nonlocal_nbr.size());
    /*
    if(my_rank==1)
    {
    printf(ANSI_COLOR_GREEN "elem_id=%d face_id exterior=%d\n" ANSI_COLOR_RESET,id,face_id);
    }
    */
     
	 if(cube_temp[0].nbr[face_id].size()==1) 	
		{
			idx=cube_temp[0].nbr[face_id].at(0);
	        cube[id].nbr[face_id].push_back(idx);
	        
	        //printf("nonlocal_nbr.size=%d\n",cube_temp[0].nonlocal_nbr.size());
	        // this section inherits the nonlocal connectivity from the parent
	        
	        if(cube_temp[0].nonlocal_nbr.size()!=0 /*&& cube_temp[0].nonlocal_nbr[0].face_tag==face_id*/)
	        {
	        
			EmptyContainer(temp_elem);
			EmptyContainer(temp_proc);
	        CollectSameTag(my_rank,cube_temp,face_id,temp_elem,temp_proc);
	        
	        //  
	       #if(0)
	        if(my_rank==1)
	        {
				//printf("elemn_id =%d face_id=%d temp_elem.size=%d\n",id,face_id,temp_elem.size());
	        printf("-----------------\n");
	        for(int l=0;l<temp_elem.size();l++)
	        {
				printf("temp_elem=%d temp_proc=%d\n",temp_elem.at(l),temp_proc.at(l));
			}
			printf("-----------------\n");
		    }    
		    #endif  
	        // first check if we have a nonlocal_nbr at that face
	       // printf("temp_elem.size()=%d\n",temp_elem.size());
	        	        
	        if(temp_elem.size()!=0)
	        {
	        AssignThisFace(my_rank,cube,id,face_id,i,j,temp_elem,temp_proc);
		    }
		    
		    
		     #if(0)
	        if(my_rank==1)
	        {
				//printf("elemn_id =%d face_id=%d temp_elem.size=%d\n",id,face_id,temp_elem.size());
	        printf("-----------------\n");
	        for(int l=0;l<cube[id].nonlocal_nbr.size();l++)
	        {
				printf(" :) elem_id=%d tag_after=%d\n",cube[id].nonlocal_nbr.at(l).elem_id,cube[id].nonlocal_nbr.at(l).face_tag);
			}
			printf("-----------------\n");
		    }    
		    #endif  
		    
		    
		    
	       }	     
		 		 
	     if(cube_temp[0].nbr[face_id].at(0)>-1)
	     {
			 cube[idx].nbr[co_face_no[face_id]].push_back(id);	
		 }
        }
        
		else if(cube_temp[0].nbr[face_id].size()==4)
		{	
			// we have to decompose 
			idx=cube_temp[0].nbr[face_id].at(K[3*i+j]);
			cube[id].nbr[face_id].push_back(idx);		
			cube[idx].nbr[co_face_no[face_id]].push_back(id);
			//printf("hereeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee\n");			
		}
	else
	{
		//printf("------------------------------------------%lu\n",cube_temp[0].nbr[face_id].size());
		printf(ANSI_COLOR_RED "element has inconsistent number of neighbors %lu \n" ANSI_COLOR_RESET,cube_temp[0].nbr[face_id].size());
		exit(0);
	}
    }
    
  
    /*
if(my_rank==1)
{
	for(unsigned int k=0;k<cube.size();k++)
	{
 
 for(j=0;j<nface;j++)
   {
     printf("nbr_id %d \n",cube[0].nbr[j].at(0));
   }

 for(unsigned int i=0;i<cube[k].nonlocal_nbr.size();i++)
   {
     //printf("************ %d nbr_id %d face tag=%d proc_id=%d\n",k,cube[k].nonlocal_nbr[i].elem_id,cube[k].nonlocal_nbr[i].face_tag,cube[k].nonlocal_nbr[i].proc_id);
   }
}

}
*/
    
    
    
    
    
    
    
    
 }
 
 }
 

 //=====================================================================
void AssignCoordsForChildren(const Cube &cube,const unsigned int face_tag,const unsigned int elem_id,Center_coords &XYZ,Message_struct (&Children)[6]); 
void PrepareChilrenForCommunication(const unsigned int id,const unsigned int nquad,const Cube &cube,const Cube &cube_temp,Center_coords &XYZ, Message_struct (&Children)[6], unsigned int *Dest)  
 {
	 unsigned int face_tag;
	// double xc_id,yc_id,r;
	//double x0,x1,y0,y1,xc,yc;
		 
    //printf("%d %d\n",cube_temp[0].nonlocal_nbr.size(),cube_temp[0].nonlocal_nbr[0].elem_id); 
		 
// if(cube[id].nonlocal_nbr.size()!=0)
// regardless of the number of nonlocal neighbors this loop will 
// include the elements on its own face for communication


//unsigned int size=cube_temp[0].nonlocal_nbr.size();

for(unsigned int i=0;i<cube_temp[0].nonlocal_nbr.size();i++)
{
	face_tag=cube_temp[0].nonlocal_nbr.at(i).face_tag;
	Dest[face_tag]=cube_temp[0].nonlocal_nbr.at(i).proc_id;
	//printf("Destination %d\n",cube_temp[0].nonlocal_nbr.at(i).proc_id);
	//printf("%d %d\n",face_tag,Dest[face_tag]);
	switch(face_tag)
	{
	case 0:
	AssignCoordsForChildren(cube,face_tag,id,XYZ,Children);
	AssignCoordsForChildren(cube,face_tag,(nquad),XYZ,Children);
	AssignCoordsForChildren(cube,face_tag,(nquad+1),XYZ,Children);
	AssignCoordsForChildren(cube,face_tag,(nquad+2),XYZ,Children);
	//printf("face_tag %d\n",face_tag);	
	break;
	
	case 1:
	AssignCoordsForChildren(cube,face_tag,(nquad+3),XYZ,Children);
	AssignCoordsForChildren(cube,face_tag,(nquad+4),XYZ,Children);
	AssignCoordsForChildren(cube,face_tag,(nquad+5),XYZ,Children);
	AssignCoordsForChildren(cube,face_tag,(nquad+6),XYZ,Children);
	//printf("face_tag %d\n",face_tag);	
	break;
	
	case 2:
	AssignCoordsForChildren(cube,face_tag,(id),XYZ,Children);
	AssignCoordsForChildren(cube,face_tag,(nquad),XYZ,Children);
	AssignCoordsForChildren(cube,face_tag,(nquad+3),XYZ,Children);
	AssignCoordsForChildren(cube,face_tag,(nquad+4),XYZ,Children);
	//printf("face_tag %d\n",face_tag);	
	break;	
		
	case 3:
	AssignCoordsForChildren(cube,face_tag,(nquad+1),XYZ,Children);
	AssignCoordsForChildren(cube,face_tag,(nquad+2),XYZ,Children);
	AssignCoordsForChildren(cube,face_tag,(nquad+5),XYZ,Children);
	AssignCoordsForChildren(cube,face_tag,(nquad+6),XYZ,Children);
	//printf("face_tag %d\n",face_tag);	
	break;
	
	case 4:
	AssignCoordsForChildren(cube,face_tag,(id),XYZ,Children);
	AssignCoordsForChildren(cube,face_tag,(nquad+2),XYZ,Children);
	AssignCoordsForChildren(cube,face_tag,(nquad+3),XYZ,Children);
	AssignCoordsForChildren(cube,face_tag,(nquad+6),XYZ,Children);
	//printf("face_tag %d\n",face_tag);	
	break;
	
	case 5:
	AssignCoordsForChildren(cube,face_tag,(nquad),XYZ,Children);
	AssignCoordsForChildren(cube,face_tag,(nquad+1),XYZ,Children);
	AssignCoordsForChildren(cube,face_tag,(nquad+4),XYZ,Children);
	AssignCoordsForChildren(cube,face_tag,(nquad+5),XYZ,Children);
	//printf("face_tag %d\n",face_tag);	
	break;			
	}	
}
	
	     

} 
 
void IsInList(Message_struct (&Children)[6],unsigned int face_tag,unsigned int elem_id, int *bol); 
void AssignCoordsForChildren(const Cube &cube,const unsigned int face_tag,const unsigned int elem_id,Center_coords &XYZ,Message_struct (&Children)[6])
{
	/*
	static const unsigned int tag1[12]={0,1,0,1,0,1,0,1,2,3,2,3};
	static const unsigned int tag2[12]={2,3,2,3,4,5,4,5,4,5,4,5};
	*/
	double xc,yc;
	unsigned int size=Children[face_tag].size();
	
	// add only if it is not in the list 
	int bol;
	IsInList(Children,face_tag,elem_id,&bol);
	
	switch(face_tag)
	{
		case 0:
		case 1:
		xc=XYZ.at(elem_id).x;
		yc=XYZ.at(elem_id).y;
		break;
		/*
		case 1:
		xc=XYZ.at(elem_id).x;
		yc=XYZ.at(elem_id).y;
		break;
		*/ 
		case 2:
		case 3:
		xc=XYZ.at(elem_id).x;
		yc=XYZ.at(elem_id).z;
		break;
		/*
		case 3:
		xc=XYZ.at(elem_id).x;
		yc=XYZ.at(elem_id).z;
		break;
		*/ 
		case 4:
		case 5:
		xc=XYZ.at(elem_id).y;
		yc=XYZ.at(elem_id).z;
		break;
		/*
		case 5:
		xc=XYZ.at(elem_id).y;
		yc=XYZ.at(elem_id).z;
		break;
		*/
		
	}
	
	
	
	if(bol==1)
	{
	Children[face_tag].push_back(Message_Struct());		
	Children[face_tag].at(size).elem_id=elem_id;
	// change here for center using face_tag
	// already have  the center coords
	/*
	xc=0.5*(cube.at(elem_id).xyz[tag1[2*face_tag]]+cube.at(elem_id).xyz[tag1[2*face_tag+1]]);		 
    yc=0.5*(cube.at(elem_id).xyz[tag2[2*face_tag]]+cube.at(elem_id).xyz[tag2[2*face_tag+1]]);
    */ 
	Children[face_tag].at(size).c1=xc;
    Children[face_tag].at(size).c2=yc;
    }
    
  //  printf("Children[%d] size=%d elem_id\n",face_tag,Children[face_tag].size(),elem_id);
   #if(0)
   if(size>0)
   {
    printf("Children[%d] size=%d %lf %lf\n",face_tag,Children[face_tag].at(size-1).elem_id,xc,yc);
   }
   #endif
} 
//======================================================================

void IsInList(Message_struct (&Children)[6],unsigned int face_tag,unsigned int elem_id, int *bol)
{
	*bol=1;
	for(unsigned int i=0;i<Children[face_tag].size();i++)
	{
		if(Children[face_tag].at(i).elem_id==elem_id)
		{
			*bol=0;
			 break;
		}
		
	}
	
}


//======================================================================
void CalculateCenter(Center_coords &XYZ,const unsigned int id,const int face_tag,double *xc);

void SearchAndUpdate(int my_rank,Cube &cube, Message_Struct *rma_buff,const int face_tag,const int co_face_tag,const int off_set,double *xc, int elem_id,int peoc_id,double *ancestor_length,Center_coords &XYZ,Message_struct &rcv_list);

void ListAffectedElements(int my_rank, const Cube &cube,const int face_tag,const int co_face_tag,const int off_set,const Message_Struct* rma_buff,Vector_Unint &affected_list,Center_coords &XYZ, double ancestor_length[3]);

void EliminateFromAffectedElements(Cube &cube,Message_Struct* rma_buff,unsigned int face_tag, Vector_Unint &affected_list,int *proc_id);

void UpdateNonlocalNbrs(int my_rank,Cube &cube, Message_Struct* rma_buff,const int off_set,const int face_tag,const int co_face_tag,Center_coords &XYZ,double ancestor_length[3],Vector_Unint &affected_list,Message_struct &rcv_list)
{
	
	unsigned int elem_id;
	double xc[2];
	int proc_id;
	clock_t start,end;

    static const unsigned int tag[12]={0,1,0,1,0,2,0,2,1,2,1,2};
	
	unsigned int coord_index;
	double xyz[3];
	
	ListAffectedElements(my_rank,cube,face_tag,co_face_tag,off_set,rma_buff,affected_list,XYZ,ancestor_length);
	
	EliminateFromAffectedElements(cube,rma_buff,face_tag,affected_list,&proc_id);

    //printf("%d\n",affected_list.size());

    start=clock();
	for(unsigned int i=0;i<affected_list.size();i++)
	{
		elem_id=affected_list.at(i);	
		//printf("%d\n",elem_id);
		CalculateCenter(XYZ,elem_id,face_tag,xc);
	
	/* or use this to avoid function calls here
	     coord_index=cube.at(elem_id).centeroid_index;
	     xyz[0]=XYZ.at(coord_index).x;
	     xyz[1]=XYZ.at(coord_index).y;
	     xyz[2]=XYZ.at(coord_index).z;
	     xc[0]=xyz[tag[2*face_tag]];
	     xc[1]=xyz[tag[2*face_tag+1]];
	  */  	 
	
		//if(my_rank==0)
		{
		// printf(ANSI_COLOR_YELLOW "=======================id =%d center=%lf %lf\n" ANSI_COLOR_RESET,elem_id,xc[0],xc[1]);	
		
		SearchAndUpdate(my_rank,cube,rma_buff,face_tag,co_face_tag,off_set,xc,elem_id,proc_id,ancestor_length,XYZ,rcv_list);			
		}
		//find closest elements 
		
		
		
	}
	
	end=clock();
	
	printf("myrank=%d spent on search and update\n",my_rank,double(end-start)/CLOCKS_PER_SEC);
	
	

}

//======================================================================
int IsInList(unsigned int elem,const int face_tag, const int off_set, const Message_Struct *rma_buff);
void FindInsideRectangle(const Cube& cube,const double *xc,const int face_tag,unsigned int *elem_id,Center_coords &XYZ, double *ancestor_length);
int IsInList(const Vector_Unint &affected_list,unsigned int elem_id);
void ListAffectedElements(int my_rank,const Cube &cube,const int face_tag,const int co_face_tag,const int off_set,const Message_Struct* rma_buff,Vector_Unint &affected_list,Center_coords &XYZ, double ancestor_length[3])
{
	
	unsigned int elem_id;
	/*
	while(affected_list.size()!=0)
	{
		affected_list.pop_back();
	}
	*/
	affected_list.clear();
//	affected_list.shrink_to_fit();
	
// if my face_tag cooredsponds to nbr co_face_tag and the updated e;ement
// in the neighboring processor is in the connectivity of the elements of the current processor 
	/*
	for(unsigned int i=0;i<cube.size();i++)
	{
		//printf("affetcet candidate=%d\n",cube[i].nonlocal_nbr.size());
		
		for(unsigned int j=0;j<cube[i].nonlocal_nbr.size();j++)
		{			
			if(cube[i].nonlocal_nbr.at(j).face_tag==co_face_tag)
			{
				elem_id=cube[i].nonlocal_nbr.at(j).elem_id;
				
				//printf("affetcet candidate=%d\n",elem_id);
				// get previous number of the elment
				
				if(IsInList(elem_id,face_tag,off_set,rma_buff)!=0)
				{
				   affected_list.push_back(i);
			    }
			    
			}
		}
	}
	*/
	double xc[2];
	//printf("SIZE=%u\n",rma_buff[6*off_set+face_tag*off_set+1].elem_id);
	// due to round off error need to specify some kind of tolerance 
	// for purturbatins needed to find the affected elements
	double prtrb=1.e-6;
	double perturbation_x[4]={prtrb,prtrb,-prtrb,-prtrb};
	double perturbation_y[4]={prtrb,-prtrb,prtrb,-prtrb};
	
for(int j=0;j<4;j++)
{
	for(unsigned int i=0;i<rma_buff[6*off_set+face_tag*off_set+1].elem_id-2;i++)
	{
		xc[0]=rma_buff[6*off_set+face_tag*off_set+2+i].c1;
		xc[1]=rma_buff[6*off_set+face_tag*off_set+2+i].c2;
		
		xc[0]=xc[0]+perturbation_x[j];
		xc[1]=xc[1]+perturbation_y[j];
		
		FindInsideRectangle(cube,xc,face_tag,&elem_id,XYZ,ancestor_length);
		//printf("elem_id=%u  center=%lf %lf\n",elem_id,xc[0],xc[1]);
		
		if(IsInList(affected_list,elem_id)==0)
		{
		   
			affected_list.push_back(elem_id);
			
	      //printf("elem_idddddddddd=%u\n",affected_list.at(0));
			
		}
	}
}
			 //printf("elem_idddddddddd=%u\n",affected_list.at(0));
	#if(0)
	if(my_rank==2)
	{
	for(unsigned int i=0;i<affected_list.size();i++)
	{
	 	printf("my_rank=%d ????????????????????  %d\n",my_rank,affected_list.at(i));
	}
    }
    #endif	
	//printf("affetcet size=%d\n",affected_list.size());
	// find the smallest element conecting these guys and add that to list
	
}

//======================================================================

void FindInsideRectangle(const Cube& cube,const double *xc,const int face_tag,unsigned int *elem_id,Center_coords &XYZ, double *ancestor_length)
{
	static const unsigned int tag1[12]={0,1,0,1,0,1,0,1,2,3,2,3};
	static const unsigned int tag2[12]={2,3,2,3,4,5,4,5,4,5,4,5};
	//static const unsigned int tag3[12]={2,3,2,3,4,5,4,5,4,5,4,5};
	
	//*elem_id=0;
	  double xyz[6];
	 //  level,coord_index;
	  
unsigned int level, coord_index;
double dx,dy,dz,idenum;
	
	for(unsigned int i=0;i<cube.size();i++)
	{
		if(cube[i].nbr[face_tag].at(0)==(-face_tag-1))
		{
			
	  level=cube.at(i).level;
	  coord_index=cube.at(i).centeroid_index;
	  //idenum=1./pow(2,level);
	  
	  TwoPowN(level,&idenum);
	  idenum=1./idenum;
	  
	  dx=ancestor_length[0]*idenum;
	  dy=ancestor_length[1]*idenum;
	  dz=ancestor_length[2]*idenum;
      xyz[0]=XYZ.at(coord_index).x-dx*0.5;
      xyz[1]=XYZ.at(coord_index).x+dx*0.5;
      xyz[2]=XYZ.at(coord_index).y-dy*0.5;
      xyz[3]=XYZ.at(coord_index).y+dy*0.5;
      xyz[4]=XYZ.at(coord_index).z-dz*0.5;
      xyz[5]=XYZ.at(coord_index).z+dz*0.5;
					
	//if((xyz[tag1[2*face_tag]]<=xc[0] && xc[0]<=xyz[tag1[2*face_tag+1]]) && (xyz[tag2[2*face_tag]]<=xc[1] && xc[1]<=xyz[tag2[2*face_tag+1]]))  

	if((xyz[tag1[2*face_tag]]<xc[0] && xc[0]<xyz[tag1[2*face_tag+1]]) && (xyz[tag2[2*face_tag]]<xc[1] && xc[1]<xyz[tag2[2*face_tag+1]]))  
		{
			*elem_id=i;
			break;
	    }
	  }
	}	

}
//======================================================================

void EliminateFromAffectedElements(Cube &cube,Message_Struct* rma_buff,unsigned int face_tag, Vector_Unint &affected_list,int *proc_id)
{
	unsigned int elem_id;
	
	
	for(unsigned int i=0;i<affected_list.size();i++)
	{
		elem_id=affected_list.at(i);
		
		for(unsigned int j=0;j<cube[elem_id].nonlocal_nbr.size();j++)
		{
			if(cube[affected_list.at(i)].nonlocal_nbr.at(j).face_tag==face_tag)
			{
				//printf("i=%d eliminated nonlocal nbr %d\n",i,affected_list.at(i));
				*proc_id=cube[elem_id].nonlocal_nbr[j].proc_id;
				cube[elem_id].nonlocal_nbr.erase(cube[elem_id].nonlocal_nbr.begin()+j);
	                        j=j-1; 			
			}			
		}
	}
	
	
	// watch out this erase statement reduces the size by one 
	
	
	// find proc ID
	
	
}

//======================================================================

void CalculateCenter(Center_coords &XYZ,const unsigned int id,const int face_tag,double *xc)
{
	/*	
	static const unsigned int tag1[12]={0,1,0,1,0,1,0,1,2,3,2,3};
	static const unsigned int tag2[12]={2,3,2,3,4,5,4,5,4,5,4,5};

	xc[0]=0.5*(cube.at(id).xyz[tag1[2*face_tag]]+cube.at(id).xyz[tag1[2*face_tag+1]]);		 
    xc[1]=0.5*(cube.at(id).xyz[tag2[2*face_tag]]+cube.at(id).xyz[tag2[2*face_tag+1]]);
    */
    
    	switch(face_tag)
	{
		case 0:
	        case 1:
		xc[0]=XYZ.at(id).x;
		xc[1]=XYZ.at(id).y;
		break;
/*
		case 1:
		xc[0]=XYZ.at(id).x;
		xc[1]=XYZ.at(id).y;
		break;
*/
		case 2:
	        case 3:
		xc[0]=XYZ.at(id).x;
		xc[1]=XYZ.at(id).z;
		break;
/*
		case 3:
		xc[0]=XYZ.at(id).x;
		xc[1]=XYZ.at(id).z;
		break;
*/
		case 4:
		case 5:
		xc[0]=XYZ.at(id).y;
		xc[1]=XYZ.at(id).z;
		break;
/*
		case 5:
		xc[0]=XYZ.at(id).y;
		xc[1]=XYZ.at(id).z;
		break;
*/		
		
	}
    
    
   // printf(ANSI_COLOR_YELLOW "=======================id =%d center=%lf %lf\n" ANSI_COLOR_RESET,id,xc[0],xc[1]);
				 		
}
//======================================================================

int IsInList(unsigned int elem_id,const int face_tag, const int off_set, const Message_Struct *rma_buff)
{
	int bol=0;
		
	for(unsigned int i=0;i<rma_buff[6*off_set+face_tag*off_set+1].elem_id-2;i++)
	{
		if(rma_buff[6*off_set+face_tag*off_set+i+2].elem_id==elem_id)
		{
			bol=1;
			break;
		}		
	}
	return(bol);	
} 
//======================================================================

int IsInList(const Vector_Unint &affected_list,unsigned int elem_id)
{
	int bol=0;
		
	for(unsigned int i=0;i<affected_list.size();i++)
	{
		
		if(affected_list.at(i)==elem_id)
		{
			bol=1;
			break;
		}
		 		
	}
	//printf("bol=%d\n",bol);
	return(bol);	
}

//======================================================================
// possible bugs in terms of face assignment
bool compare(int i,int j) { return (i<j); }
#if(0)
void SearchAndUpdate(int my_rank,Cube &cube, Message_Struct *rma_buff,const int face_tag,const int co_face_tag,const int off_set,double *xc, int elem_id,int proc_id)
{
	
	double r=1000.0,delta[2];
	unsigned int nbr,nbr0;
	int size;
	int count=0;
	unsigned int index=0;
	Vector_Int tmp;
	double tol=1.e-6;
	
	unsigned int BIG=1e7;
	//double BIG=1.0e8;
	//printf("BIG =%u\n",BIG);
	//printf("BIG =%lf\n",BIG);
	 //printf("\n face tag=%d recver rank=%d ==============================index size=%d\n",face_tag,rma_buff[6*off_set+off_set*face_tag].elem_id,rma_buff[6*off_set+off_set*face_tag+1].elem_id);
		
	for(unsigned int i=0;i<rma_buff[6*off_set+off_set*face_tag+1].elem_id-2;i++)
    {
	 //  if(rma_buff[6*off_set+off_set*face_tag+2+i].elem_id!=BIG)
	  {
	        delta[0]=xc[0]-rma_buff[6*off_set+off_set*face_tag+2+i].c1;
			delta[1]=xc[1]-rma_buff[6*off_set+off_set*face_tag+2+i].c2;
			
			if((delta[0]*delta[0])+(delta[1]*delta[1])<r)
			{
				r=(delta[0]*delta[0])+(delta[1]*delta[1]);
				// 						
				index=i;		
		    }			
    }
    }
    
    //printf("index=%d elem_id=%d r=%lf elem_id %d my coord=%lf %lf\n",index,rma_buff[6*off_set+off_set*face_tag+2+index].elem_id,r,elem_id,xc[0],xc[1]);
   //printf("id=%d R=%lf\n",elem_id,r);
    nbr=rma_buff[6*off_set+off_set*face_tag+2+index].elem_id;	
	tmp.push_back(nbr);
 	//rma_buff[6*off_set+off_set*face_tag+2+index].elem_id=BIG;	
	
//printf("id=%d R=%lf nbr=%d\n",elem_id,r,nbr);
    
    //printf("tmp=%d\n",nbr); 				    
	// find the other three neighbors in case it has any
		
	if(r>tol)
	{
	//tmp.pop_back();	
	
	//for(unsigned int i=0;i<rma_buff[6*off_set+off_set*face_tag+1].elem_id-2 && i!=index;i++)
	for(unsigned int i=0;i<rma_buff[6*off_set+off_set*face_tag+1].elem_id-2 ;i++)
    {
		//if(rma_buff[6*off_set+off_set*face_tag+2+i].c1!=BIG)
		//if(rma_buff[6*off_set+off_set*face_tag+2+i].elem_id!=BIG)
		{
	        delta[0]=xc[0]-rma_buff[6*off_set+off_set*face_tag+2+i].c1;
			delta[1]=xc[1]-rma_buff[6*off_set+off_set*face_tag+2+i].c2;
			 						
			
			if(((delta[0]*delta[0])+(delta[1]*delta[1])-r)<tol)
			{
				nbr=rma_buff[6*off_set+off_set*face_tag+2+i].elem_id;
			//printf("R=%16.16lf\n",(delta[0]*delta[0])+(delta[1]*delta[1]));	
							
			
			/*	if(my_rank==10)
				{
				printf("tmp=%d %lf %lf %lf %lf\n",nbr,xc[0],xc[1],rma_buff[6*off_set+off_set*face_tag+2+i].c1,rma_buff[6*off_set+off_set*face_tag+2+i].c2); 						
			   }
			 */ 
			   if(i!=index)
			    {
				tmp.push_back(nbr);
			//	rma_buff[6*off_set+off_set*face_tag+2+i].elem_id=BIG;
				
			    }
		    }
		   } 		
      }		
	}
	
		
	if(tmp.size()!=1 && tmp.size()!=4)
	{
		printf(ANSI_COLOR_RED "Error in Assigning Neighbors, 4:1 balance not preserved in SearchAndUpdate tmp_size= %lu\n" ANSI_COLOR_RESET,tmp.size());
		printf(ANSI_COLOR_RED "my_rank %d elements %lu %lu %d\n" ,my_rank,tmp.at(0),tmp.at(1),face_tag);
		exit(0);
			
	}
	
	//printf(ANSI_COLOR_RED "Error in Assigning Neighbors, 4:1 balance not preserved in SearchAndUpdate\n");
	
	  
	if(tmp.size()==4)
	{
		std::sort (tmp.begin(),tmp.end(), compare);
	}

	 #if(1)
	
	for(unsigned int i=0;i<tmp.size();i++)
    {
		size=cube[elem_id].nonlocal_nbr.size();
		cube[elem_id].nonlocal_nbr.push_back(Nonlocal_Nbr());
		cube[elem_id].nonlocal_nbr[size].elem_id=tmp.at(i);
		cube[elem_id].nonlocal_nbr[size].face_tag=face_tag;
		cube[elem_id].nonlocal_nbr[size].proc_id=proc_id;
		// add proc id and face_tag
	} 
	
	#endif

  
  //printf("size of tmp=%d\n",tmp.size());
   tmp.clear();
  // tmp.shrink_to_fit();
}
#endif
//======================================================================

// modify this, first compare the centers, if the coinside then asssign and eliminate from
// list otherwise construct the dges and look for 4 instead of 1 
bool compareStruct(Message_Struct XYZ1, Message_Struct XYZ2) { return XYZ1.c1 < XYZ2.c1; };

void SearchAndUpdate(int my_rank,Cube &cube, Message_Struct *rma_buff,const int face_tag,const int co_face_tag,const int off_set,double *xc, int elem_id,int proc_id,double *ancestor_length,Center_coords &XYZ,Message_struct &rcv_list)
{
	
	//double r=1000.0,delta[2];
	unsigned int nbr,nbr0;
	int size;
	int count=0;
	unsigned int index=0;
	Vector_Int tmp;
	double tol=1.e-6;
	
		
	unsigned int level, coord_index;
    double dx,dy,dz,idenum;
	double Xc[4];
	double xyz[6];
			
	  level=cube.at(elem_id).level;
	  coord_index=cube.at(elem_id).centeroid_index;
	  //idenum=1./pow(2,level);
	  TwoPowN(level,&idenum);
	  idenum=1./idenum;
	  
	  dx=ancestor_length[0]*idenum;
	  dy=ancestor_length[1]*idenum;
	  dz=ancestor_length[2]*idenum;
      xyz[0]=XYZ.at(coord_index).x-dx*0.5;
      xyz[1]=XYZ.at(coord_index).x+dx*0.5;
      xyz[2]=XYZ.at(coord_index).y-dy*0.5;
      xyz[3]=XYZ.at(coord_index).y+dy*0.5;
      xyz[4]=XYZ.at(coord_index).z-dz*0.5;
      xyz[5]=XYZ.at(coord_index).z+dz*0.5;
      
    static const unsigned int tag1[12]={0,1,0,1,0,1,0,1,2,3,2,3};
	static const unsigned int tag2[12]={2,3,2,3,4,5,4,5,4,5,4,5};	
	
	Xc[0]=xyz[tag1[2*face_tag]];
	Xc[1]=xyz[tag1[2*face_tag+1]];
	Xc[2]=xyz[tag2[2*face_tag]];
	Xc[3]=xyz[tag2[2*face_tag+1]];
		
	
	
	  unsigned int idx=6*off_set+off_set*face_tag+2;
	//double BIG=1.0e8;
	//printf("BIG =%u\n",BIG);
	//printf("BIG =%lf\n",BIG);
	 //printf("\n face tag=%d recver rank=%d ==============================index size=%d\n",face_tag,rma_buff[6*off_set+off_set*face_tag].elem_id,rma_buff[6*off_set+off_set*face_tag+1].elem_id);
#if(CPP==1)
		
	for(unsigned int i=0;i<rma_buff[6*off_set+off_set*face_tag+1].elem_id-2;i++)
    {
	  if(Xc[0]<=rma_buff[idx+i].c1 && Xc[1]>=rma_buff[idx+i].c1 && Xc[2]<=rma_buff[idx+i].c2 && Xc[3]>=rma_buff[idx+i].c2)
	  // if(Xc[0]<rma_buff[idx+i].c1 && Xc[1]>rma_buff[idx+i].c1 && Xc[2]<rma_buff[idx+i].c2 && Xc[3]>rma_buff[idx+i].c2)
	  {
	        nbr=rma_buff[idx+i].elem_id;	
	       // printf("inside if\n");
	        tmp.push_back(nbr);		
      }
    }
#endif
    
   #if(CPP==0)
    std::sort(rcv_list.begin(), rcv_list.end(), &compareStruct);
    
    for(unsigned int i=0;i<rcv_list.size();i++)
    {
		
	  if(fabs(xc[0]-rcv_list[i].c1)<1.e-12 && fabs(xc[1]-rcv_list[i].c2)<1.e-12)
	  {
		   nbr=rcv_list[i].elem_id;	
	       // printf("inside if\n");
	        tmp.push_back(nbr);
	        
	        rcv_list.erase(rcv_list.begin()+i);
	        i--;
	        break;
		  
	  }
	  else if(Xc[0]<=rcv_list[i].c1 && Xc[1]>=rcv_list[i].c1 && Xc[2]<=rcv_list[i].c2 && Xc[3]>=rcv_list[i].c2)
	  // if(Xc[0]<rma_buff[idx+i].c1 && Xc[1]>rma_buff[idx+i].c1 && Xc[2]<rma_buff[idx+i].c2 && Xc[3]>rma_buff[idx+i].c2)
	  {
	        nbr=rcv_list[i].elem_id;	
	       // printf("inside if\n");
	        tmp.push_back(nbr);
	        
	        
	        if(tmp.size()==4)
	        {
			rcv_list.erase(rcv_list.begin()+i);
	        i--;
	        break;
	        
		    }
		    
      }
    }
    #endif
    
    
    //printf("index=%d elem_id=%d r=%lf elem_id %d my coord=%lf %lf\n",index,rma_buff[6*off_set+off_set*face_tag+2+index].elem_id,r,elem_id,xc[0],xc[1]);
   //printf("id=%d R=%lf\n",elem_id,r);
   
	
		/*
	if(tmp.size()!=1 && tmp.size()!=4)
	{
		printf(ANSI_COLOR_RED "Error in Assigning Neighbors, 4:1 balance not preserved in SearchAndUpdate tmp_size= %lu\n" ANSI_COLOR_RESET,tmp.size());
		printf(ANSI_COLOR_RED "my_rank %d elements %lu %lu %d\n" ,my_rank,tmp.at(0),tmp.at(1),face_tag);
		exit(0);
			
	}
	*/
	//printf(ANSI_COLOR_RED "Error in Assigning Neighbors, 4:1 balance not preserved in SearchAndUpdate\n");
	
	  
	if(tmp.size()==4)
	{
		std::sort (tmp.begin(),tmp.end(), compare);
	}
	else if(tmp.size()!=1)
	{
		printf(ANSI_COLOR_RED "Error in Assigning Neighbors, 4:1 balance not preserved in SearchAndUpdate tmp_size= %lu\n" ANSI_COLOR_RESET,tmp.size());
		printf(ANSI_COLOR_RED "my_rank %d elements %lu %lu %d\n" ,my_rank,tmp.at(0),tmp.at(1),face_tag);
		exit(0);
		
	}

	 #if(1)
	
	for(unsigned int i=0;i<tmp.size();i++)
    {
		size=cube[elem_id].nonlocal_nbr.size();
		cube[elem_id].nonlocal_nbr.push_back(Nonlocal_Nbr());
		cube[elem_id].nonlocal_nbr[size].elem_id=tmp.at(i);
		cube[elem_id].nonlocal_nbr[size].face_tag=face_tag;
		cube[elem_id].nonlocal_nbr[size].proc_id=proc_id;
		// add proc id and face_tag
	} 
	
	#endif

  
  //printf("size of tmp=%d\n",tmp.size());
   tmp.clear();
 //  tmp.shrink_to_fit();
}




//======================================================================

//  for later to do some templates 
/*
template<class VECTOR>
void Vector_Clear(VECTOR &V)
{
	while(V.size()!=0)
	{
		V.pop_back();
	}
	
}
*/

 void CountSameTag(int my_rank,const Cube &cube,const unsigned int face_tag, const unsigned int id, int *bol ,unsigned int *nbr, unsigned int *proc_id)
 {
	 int count=0;
	 
	 *bol=1;
	 *nbr=1e6;
	 
	 
	  for(unsigned int k=0;k<cube[id].nonlocal_nbr.size();k++)
	        {
				if(cube[id].nonlocal_nbr.at(k).face_tag==face_tag)
				{
					count++;
					
					*nbr=cube[id].nonlocal_nbr.at(k).elem_id;
					*proc_id=cube[id].nonlocal_nbr.at(k).proc_id;
					 // count>1 means my level is higher than my nbrs level
					if(count>1)
					{
						*bol=0;
						break;
					}
				}
			}
			
			if(count==0)
			{
				*bol=0;
			}
			/*
			if(count==1)
			{
				
			*bol=1;	
				
			}
			*/
			
		//	printf(ANSI_COLOR_RED "bol=%d nbr=%d nonlocal_size=%d count %d\n" ANSI_COLOR_RESET,*bol,*nbr,cube[id].nonlocal_nbr.size(),count);
		
}


//======================================================================
#if(0)
void WriteOut(int my_rank,const Cube &cube)
{
FILE *fo=NULL;
 char buf[12];
 sprintf(buf, "%d.txt", my_rank);
  fo=fopen(buf,"w");
  fprintf(fo,"%lu\n",cube.size());
  
  for(unsigned int i=0;i<cube.size();i++)
  {
	  for(unsigned int j=0;j<6;j++)
	  {
	  fprintf(fo,"%lf\t",cube[i].xyz[j]);
      }
      fprintf(fo,"\n");
     for(unsigned int j=0;j<6;j++)
	  {
		  for(unsigned int k=0;k<cube[i].nbr[j].size();k++)
		{
	  fprintf(fo,"%d\t",cube[i].nbr[j].at(k));
        }
	  fprintf(fo,"\n");
      }
      // 1 is added as a bolean to know that an element is a processor boundary element
      	
      	if(cube[i].nonlocal_nbr.size()==0)
        {
			fprintf(fo,"%d\n",0);
		}
		else
		{
			fprintf(fo,"%lu\n",cube[i].nonlocal_nbr.size());
		}
      	 
      	  for(unsigned int k=0;k<cube[i].nonlocal_nbr.size();k++)
		{
	  fprintf(fo,"%d %d %d\n",cube[i].nonlocal_nbr.at(k).elem_id,cube[i].nonlocal_nbr.at(k).face_tag,cube[i].nonlocal_nbr.at(k).proc_id);
        }
        
        
        
	  //fprintf(fo,"\n");
      
      
  }
  
fclose(fo);

}

//======================================================================

void ReadMesh(int my_rank,Cube &cube)
{
	FILE *fo=NULL;
	int bol,size;
	
	if(my_rank==0)
	{
		fo=fopen("0.txt","r");
	}
	else
	{
		fo=fopen("1.txt","r");
	}
	
	int cube_size;
	
	fscanf(fo,"%d",&cube_size);
	
	
	for(int i=1;i<cube_size;i++)
	{
		cube.push_back(cube_data());
	}
	
	//printf("myrank=%d cubesize=%d %d\n",my_rank,cube_size,cube.size());
//		

	for(int i=0;i<cube_size;i++)
   {
	   while(cube[i].nonlocal_nbr.size()!=0)
	   {
		   cube[i].nonlocal_nbr.pop_back();
	   }
	   for(int j=0;j<6;j++)
	   {
		 while(cube[i].nbr[j].size()!=0)
	   {
		   cube[i].nbr[j].pop_back();
	   }
	      
	   }
   }
	
	for(int i=0;i<cube_size;i++)
   {  
	 //  int i=0;
	 for(int j=0;j<6;j++)
	 {  
	   fscanf(fo,"%lf", &(cube[i].xyz[j]));
	  // printf("myrank=%d %lf\n",my_rank,cube[i].xyz[j]);
	 }
	   
	 for(int j=0;j<6;j++)
	 {  
		cube[i].nbr[j].push_back(0);
	   fscanf(fo,"%d", &(cube[i].nbr[j].at(0)));
	   //printf("myrank=%d %d\n",my_rank,cube[i].nbr[j].at(0));
	 }
	 
	 fscanf(fo,"%d", &(bol));
	 //printf(ANSI_COLOR_GREEN "myrank=%d bol=%d\n" ANSI_COLOR_RESET,my_rank,bol);
	 
	 if(bol!=0)
	 {
		 for(int j=0;j<bol;j++)
		 {
		size=cube[i].nonlocal_nbr.size();	 
		cube[i].nonlocal_nbr.push_back(Nonlocal_Nbr());
		fscanf(fo,"%d %d %d",&(cube[i].nonlocal_nbr[size].elem_id),&(cube[i].nonlocal_nbr[size].face_tag),&(cube[i].nonlocal_nbr[size].proc_id));
		//printf("myrank=%d bol=%d\n",i,j);
	     }
	 }
 }
	 if(cube_size==8)
	 {
		 for(int i=0;i<8;i++)
		 {
			 cube[i].level=1;
		 }
	 }
	 else
	 {
		  cube[0].level=0;
	 }
	
	
	
}

#endif

void IsInList(const unsigned int elem,const Vector_Unint &ranlist,int *bol);

void ran_list(Cube &cube,Vector_Unint &ranlist,unsigned int level)
{

// empty the list first
	 while(ranlist.size()!=0)
 	{
	 ranlist.pop_back();
	 }			
//
// unsigned int size=cube.size();
//
printf(ANSI_COLOR_RED "size=%lu\n" ANSI_COLOR_RESET,cube.size());
//
unsigned int size=cube.size();
//size=130;
unsigned int elem;
unsigned int current_level=0;
int bol;

if(size<level)
{
  printf(ANSI_COLOR_RED "Random selected number is out of range\n" ANSI_COLOR_RESET);
  exit(0);
}
	
 while(current_level<level)
 {
	elem=rand() % size;
	IsInList(elem,ranlist,&bol);
	if(bol==0)
	{	
	ranlist.push_back(elem);
	current_level++;
	}
 }

for(unsigned int i=0;i<ranlist.size();i++)
{
if(ranlist.at(i)>=cube.size())
{
printf("invalid cube number in the list\n");
exit(0);
}
}


}

void IsInList(const unsigned int elem,const Vector_Unint &ranlist,int *bol)
{
*bol=0;
for(unsigned int i=0;i<ranlist.size();i++)
{
if(elem==ranlist.at(i))
{
*bol=1;
break;
}


}

}


void ReadInGeom(double **xyz,int *nn)
{
	
	FILE *fg=NULL;
	//fg=fopen("Geom/0012.txt","r");
	//fg=fopen("Geom/cylinder.txt","r");
	fg=fopen("Geom/0012_dense.txt","r");
	if(fg==NULL)
	{
		printf(ANSI_COLOR_RED "cound not open geometry file\n" ANSI_COLOR_RESET);
		exit(0);
	}
	const int bdim = 132;
char buff[bdim];
int i,b;

//=========================== read number of nodes ===================
  

  fgets(buff,bdim,fg);
  sscanf(buff,"%d",nn);
  printf("\nNumber of points = %d\n",(*nn));
  fflush(stdout);

//========================allocate for coordinates =====================

  (*xyz) = (double*)malloc(3*(*nn)*sizeof(double));
  
//======================== read in coordinates =========================

  for (i=0; i < (*nn); i++)
  {
     fgets(buff,bdim,fg);
     sscanf(buff,"%lf %lf %lf",&((*xyz)[3*i]),&((*xyz)[3*i+1]),&((*xyz)[3*i+2]));
     //printf("i= %d x= %lf, y= %lf\n",i, ((*xyz)[3*i]),((*xyz)[3*i+1]));
     // for 0012
     
     (*xyz)[3*i]=(*xyz)[3*i]-0.5;
     (*xyz)[3*i+1]=(*xyz)[3*i+1];
     /*
     (*xyz)[3*i]=(*xyz)[3*i]*0.01;
     (*xyz)[3*i+1]=(*xyz)[3*i+1]*0.01;
     (*xyz)[3*i+2]=(*xyz)[3*i+2]*0.01;
     */
  }

	
}


void IsInsideSolid(const Cube& cube,const Center_coords &XYZ, Vector_Unint& refine_list,const double *geom_xyz, double *ancestor_length,unsigned int n)
{

double xc,yc,zc;
unsigned int i,j,level,coord_index;
int a,b,c;
double xyz[6],dx,dy,dz,idenum;
const static double half=0.5;

for(i=0;i<cube.size();i++)
 {
	 
	 
	  level=cube.at(i).level;
	  coord_index=cube.at(i).centeroid_index;
	  //idenum=1./pow(2,level);
	  TwoPowN(level,&idenum);
	  idenum=1./idenum;
	  //idenum=1./(double)2<<level;
	  
	  dx=ancestor_length[0]*idenum;
	  dy=ancestor_length[1]*idenum;
	  dz=ancestor_length[2]*idenum;
     
      xyz[0]=XYZ.at(coord_index).x-dx*half;
      xyz[1]=XYZ.at(coord_index).x+dx*half;
      xyz[2]=XYZ.at(coord_index).y-dy*half;
      xyz[3]=XYZ.at(coord_index).y+dy*half;
      xyz[4]=XYZ.at(coord_index).z-dz*half;
      xyz[5]=XYZ.at(coord_index).z+dz*half;
	 
	 for(j=0;j<n;j++)
	 {
	a=geom_xyz[3*j]> xyz[0] && geom_xyz[3*j]< xyz[1];
	b=geom_xyz[3*j+1]> xyz[2] && geom_xyz[3*j+1]< xyz[3];
	c=geom_xyz[3*j+2]> xyz[4] && geom_xyz[3*j+2]< xyz[5];
	
	if(a && b && c)
	{
		refine_list.push_back(i);
		break;
	}
   }
 }	
		
}

//======================================================================

