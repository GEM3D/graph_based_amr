//#include "/home/jhasbestan/GEM_AMR/src/include/typedefs.h"
#include "typedefs.h"
#define nvar 1

void Save(const int id,const Cube& cube,Cube& cube_temp);

void EmptyTempCube(Cube& cube_temp);

void PrepareChilrenForCommunication(const unsigned int id,const unsigned int nquad,const Cube &cube,const Cube &cube_temp,Center_coords &XYZ ,Message_struct (&Children)[6],unsigned int *Dest);

void ConstructChildrenCoords(const int id,const Cube &cube_temp,const unsigned int *map_index, Cube &cube,Center_coords& XYZ,double ancestor_length[3]);

void RemoveParentFromConn(const int id,const Cube &cube_temp,const unsigned int* co_face_no,Cube &cube);

void AssignInteriorElemConn(const unsigned int id,const unsigned int nquad,const unsigned int *map_index,const unsigned int *I,Cube &cube);

void FindAdjacentFace(const Cube& cube,const Cube& cube_temp,const int id,const int tag,const int face_id, int *neigbor_id);

void AssignExteriorFace(int my_rank,Cube &cube,const Cube &cube_temp,const unsigned int *co_face_no,const unsigned int *map_index,const unsigned int*I);

void UpdateNonlocalNbrs(int my_rank,Cube &cube, Message_Struct* rma_buff,const int off_set,const int face_tag,const int co_face_tag,Center_coords &XYZ,double ancestor_length[3],Vector_Unint &affected_list,Message_struct &rcv_list);

void AmrRefineParallel(const int My_rank,const MPI_Comm comm,Cube& cube,const Vector_Unint &ordered_refine_list,const int L,const int M, const int N,Center_coords &XYZ, double ancestor_length[3])
{
  unsigned int i,j,k,l1;	
  unsigned int nquad=cube.size();
  unsigned int index=nquad;	
  unsigned int q_size=L*M*N*nvar;
 
  static const unsigned int co_face_no[6]={1,0,3,2,5,4};
  static const unsigned int I[24]={1,3,5,1,3,4,1,2,4,1,2,5,0,3,5,0,3,4,0,2,4,0,2,5}; 
  unsigned int map_index[8];
 
  int Comm_size;
  MPI_Comm_size(comm, &Comm_size); 
 
  // There are only 3 faces exposed, the other three were assigned before
  // same mapping for upper elements only the first index changes
  //  each three entries in I coorespond to the each exposed face of each newly
  // generated element
  
  static const unsigned int I1[24]={0,2,4,0,2,5,0,3,5,0,3,4,1,2,4,1,2,5,1,3,5,1,3,4};
 
  // elem list is defined to store all the children such that we can use this in communication
  // with other processors, the children is only for boundary elements nbring other process
 
  Message_struct Children[6];
 
  // define MPI_Data type
 
 unsigned int my_rank=(unsigned int)My_rank;
  MPI_Datatype MPI_Message_strct;

  //MPI_Type_contiguous(sizeof(Message_Struct), MPI_BYTE, &MPI_Message_strct);
  /*
    MPI_Datatype oldtype[3] = {MPI_UNSIGNED,MPI_DOUBLE,MPI_DOUBLE};
    int blocklen[3] = {1,1,1};
    MPI_Aint disp[3];
    disp[0] = 0;
    disp[1] = 4;
    disp[2] = 12;
   
    MPI_Type_create_struct(3, blocklen, disp, oldtype, &MPI_Message_strct);
    MPI_Type_commit(&MPI_Message_strct);
    int size1;
    int  MPI_Message_strct_size;
    MPI_Type_size(MPI_Message_strct,&MPI_Message_strct_size);
    int unsigned_size;
    MPI_Type_size(MPI_UNSIGNED,&unsigned_size);
    // size of does not work for MPI_Types
    printf("MPI_Type_Size=%d %d Unsigned =%d MPI_Unsigned=%d\n",sizeof(MPI_Message_strct),MPI_Message_strct_size,sizeof(MPI_UNSIGNED),unsigned_size);
  */
  // do not count on the bitsizes of a structure

  Message_Struct a[2];
  MPI_Datatype oldtype[2] = {MPI_UNSIGNED,MPI_DOUBLE};
  int blocklen[2] = {1,2};
  MPI_Aint disp[2];   
  MPI_Get_address(&a[0].elem_id,&disp[0]);
  MPI_Get_address(&a[0].c1,&disp[1]);
  disp[1]=disp[1]-disp[0];
  disp[0]=0;
    
  //printf("displacement= %d, %d\n",disp[0],disp[1]);
    
  MPI_Type_create_struct(2, blocklen, disp, oldtype, &MPI_Message_strct);
  MPI_Type_commit(&MPI_Message_strct);
  int size1;
  int  MPI_Message_strct_size;
  //MPI_Type_size(MPI_Message_strct,&MPI_Message_strct_size);
  //int unsigned_size;
  //MPI_Type_size(MPI_UNSIGNED,&unsigned_size);
  MPI_Aint S1,S2;
  MPI_Get_address(&a[1],&S1);
  MPI_Get_address(&a[0],&S2);
  MPI_Message_strct_size=S1-S2;
    
  /*
    I1[0]=0,I1[1]=2,I1[2]=4;
    I1[3]=0,I1[4]=2,I1[5]=5;
    I1[6]=0,I1[7]=3,I1[8]=5;
    I1[9]=0,I1[10]=3,I1[11]=4;
    I1[12]=1,I1[13]=2,I1[14]=4;
    I1[15]=1,I1[16]=2,I1[17]=5;
    I1[18]=1,I1[19]=3,I1[20]=5;
    I1[21]=1,I1[22]=3,I1[23]=4;   
  */
 
  // index2 is generated to eliminate 3 if's
  
  static const int index2[6]={0,0,1,1,2,2};
  
  // bottom and top planes both faces need xy coordinates (index2[0] and index2[1]) >>> tag[0]
  //left and right need xz (index2[2] and index2[3]) >>> tag[1]
  // back and front faces need yz (index2[4] and index2[5]) >>> tag[2]
 
  // define aa temporary cube and save current object to a temporary object 
  // note that adaptation level increase by one 
 
  Cube cube_temp;
  cube_temp.push_back(cube_data());
 
  /*Save(0,cube,cube_temp);

  //if(my_rank==0)
  {
  printf("%d=====================================\n",my_rank);
  for(int i=0;i<6;i++)
  {
  for(int j=0;j<cube_temp[0].nbr[i].size();j++)
  {
  printf("%d\t",cube_temp[0].nbr[i].at(j));	 
  }
  printf("\n");	 
  }
  printf("----------------------------------------\n");
  for(int i=0;i<cube_temp[0].nonlocal_nbr.size();i++)
  {
  printf("%d %d %d\t",cube_temp[0].nonlocal_nbr[i].elem_id,cube_temp[0].nonlocal_nbr[i].proc_id,cube_temp[0].nonlocal_nbr[i].face_tag);	 
  }
  printf("\n=====================================\n");	  
  }
  */
clock_t msg_start,msg_end;

msg_start=clock(); 	

  unsigned int id;
  unsigned int comm_size=(unsigned int)Comm_size;
  unsigned int Dest[6]={comm_size,comm_size,comm_size,comm_size,comm_size,comm_size};

  for(unsigned int i=0;i<ordered_refine_list.size();i++)
    {
      id=ordered_refine_list.at(i);
	 
      // Save the cube to be modified 
	 
      Save(id,cube,cube_temp);
	 
#if(0)
      if(my_rank==1)
	{
	  for(int k=0;k<6;k++)
	    {
	      for(unsigned int j=0;j<cube_temp[0].nbr[k].size();j++)
		{
		  printf("%u %d\n",j,cube_temp[0].nbr[k].at(j));
		}
	    }
	  for(unsigned int j=0;j<cube_temp[0].nonlocal_nbr.size();j++)
	    {
	      printf("---------------%u %u %u\n",j,cube_temp[0].nonlocal_nbr.at(j).elem_id,cube_temp[0].nonlocal_nbr.at(j).face_tag);
	    }
	}	
#endif 
      // add 7 elements to the number of cubes
 
      // initial size
      nquad=cube.size(); 
 
      // allocate 8 more 
      for(unsigned int j=0;j<7;j++) 
	{
	  cube.push_back(cube_data());
	  // and need some more memory for XYZ
	  XYZ.push_back(Center_Coords());
	} 
 
      //  set adaptation level, adaptation level might be calcutable given the size of the 
      // initial and current size, due to diviion by 2 at each level
 
      cube[id].level= cube_temp[0].level+1;
      cube[id].centeroid_index= cube_temp[0].centeroid_index;
       
      for(unsigned int j=0;j<7;j++)
	{
	  cube.at(nquad+j).level=cube_temp[0].level+1;
	  cube.at(nquad+j).centeroid_index=nquad+j; 
	}
 
      // eliminate the old neighbor data from the parent element since now 
      //parent will be ith child (kind of like each element pointing to itself )
 
      for(unsigned int j=0;j<6;j++)
	{
	  cube.at(id).nbr[j].clear();
	  cube.at(id).nonlocal_nbr.clear();
	  
	}
 
      // form the map index
 
      map_index[0]=id;
 
      //Children.push_back(map_index[0]);
  
      for(unsigned int j=1;j<8;j++){
		map_index[j]=nquad+(j-1);
      }
	 
      ConstructChildrenCoords(id,cube_temp,map_index,cube,XYZ,ancestor_length);
 
      RemoveParentFromConn(id,cube_temp,co_face_no,cube);
 
      // watch out Interior gets >>>> I
      AssignInteriorElemConn(id,nquad,map_index,I,cube);
 
      // watch out exterior gets >>>> I1
      AssignExteriorFace(my_rank,cube,cube_temp,co_face_no,map_index,I1); 
 
      PrepareChilrenForCommunication(id,nquad,cube,cube_temp,XYZ,Children,Dest);

	  EmptyTempCube(cube_temp);
    }

  //printf("rank=%d\n",my_rank);
  
  
  msg_end=clock();
    
    //if(my_rank==0)
    {
		printf(ANSI_COLOR_RED"my_rank %d time spend before message passing %16.16lf \n" ANSI_COLOR_RESET,my_rank,double(+msg_end-msg_start)/CLOCKS_PER_SEC);
	}

#if(0)
if(my_rank==0)
{
  for(int i=0;i<6;i++)
    {
     // printf("i=%d %d\n",i,Children[i].size());	  
     printf("++++++++++++++\n"); 
      for(int j=0;j<Children[i].size();j++)
	{
	  printf("%d %lf %lf\n",Children[i].at(j).elem_id,Children[i].at(j).c1,Children[i].at(j).c2); 
	}
      printf("++++++++++++++\n"); 
    }
}



  for(int i=0;i<cube.size();i++)
    {
  
      for(int j=0;j<6;j++)
	{
	  printf("%lf\t",cube[i].xyz[j]); 
	}
      printf("\n"); 
    }


  if(my_rank==0)
    {
      printf("----------------------------------------------\n");
      for(int i=0;i<6;i++)
	{
	  //printf("i=%d %d\n",i,Children[i].size());
	  //if(Children[i].size()!=0)
	  {
	    for(int j=0;j<Children[i].size();j++)
	      {
		printf("%d %lf %lf\t",Children[i].at(j).elem_id,Children[i].at(j).c1,Children[i].at(j).c2); 
	      }
	    printf("\n"); 
	  }
	}
      printf("----------------------------------------------\n"); 
    }
#endif

  //=====================================================================
  /*
    communication part 
  */ 
  //=====================================================================
  if(my_rank==0)
    {
      for(unsigned int j=0;j<cube.size();j++)
	{
	  // printf(" cube.size =%d cube[j].nonlocal_nbr.size=%d\n",cube.size(),cube[j].nonlocal_nbr.size());
	  /*
	    for(int k=0;k<6;k++)
	    {
	    //printf("nbr.size=%d\n",cube[j].nbr[k].size());
	    for(i=0;i<cube[j].nbr[k].size();i++)
	    {
	    printf("**********my_rank=%d nbr_indx =%d nbr_id %d \n",j,k,cube[j].nbr[k].at(i));
	    } 
   
   
	    }
	    printf("\n");

 
	    for(i=0;i<cube[j].nonlocal_nbr.size();i++)
	    {
	    printf("**********my_rank=%d cube_id =%d nbr_id %d face tag=%d proc_id=%d\n",my_rank,j,cube[j].nonlocal_nbr[i].elem_id,cube[j].nonlocal_nbr[i].face_tag,cube[j].nonlocal_nbr[i].proc_id);
	    }
	  */
	}

    }
#if(1)
 
 
 
  Message_Struct *rma_buff=NULL; 
  unsigned int my_buff_size,max_buff_size,off_set,destination,face_tag,co_face_tag,displacement1,displacement2,disp_unit;
 
  Vector_Unint affected_list;
 
  MPI_Win win;
 
  my_buff_size=0;
	
  for(unsigned int j=0;j<6;j++)
    {
      if(my_buff_size<Children[j].size())
	{
	  my_buff_size=Children[j].size();
	}
    }
	

	
  MPI_Allreduce(&my_buff_size,&max_buff_size,1,MPI_INT,MPI_MAX,comm);
  
  msg_start=clock();
	
  // need to know tha maximum size for the window for allocation
  // now that we know the worst case scenario (biggest size of matrix)
  // we can allocate the window and perform RMA
 
	
  off_set=max_buff_size+2;
  // there are send and recieves to keep them separate allocate twice as much as you need for now, this can be improved *2  in
  // the following is due to this reason
	
  max_buff_size=6*(off_set)*2;	 
  
  Message_struct rcv_list;
	
  //printf("my_rank=%d max_buff_size%d %lu %lu\n",my_rank,max_buff_size,sizeof(MPI_Message_strct),sizeof(Message_Struct));
	
  printf(ANSI_COLOR_CYAN "my_rank=%d %d\n" ANSI_COLOR_RESET,my_rank,max_buff_size);
	
  MPI_Alloc_mem(MPI_Message_strct_size*max_buff_size, MPI_INFO_NULL,&rma_buff);   
  disp_unit=MPI_Message_strct_size;
  MPI_Win_create(rma_buff,max_buff_size*MPI_Message_strct_size,disp_unit,MPI_INFO_NULL,comm,&win); 
      
  for(unsigned int i=0;i<max_buff_size;i++)
    {
      rma_buff[i].elem_id=0;
      rma_buff[i].c1=0.0;
      rma_buff[i].c2=0.0;
    }
      
      
  // assign something max_buff_size
       		
  for(int i=0;i<6;i++)
    {
      if(Children[i].size()!=0)
	{
		
	  my_buff_size=Children[i].size()+2;
	  destination=Dest[i];
	  rma_buff[i*off_set].elem_id=destination;
	  rma_buff[i*off_set+1].elem_id=my_buff_size;
		
	  //printf("i=%d Dest=%d my_buff_size=%d %d %d offset=%d\n",i,Dest[i],my_buff_size,Children[i].size(),MPI_Message_strct_size*max_buff_size,i*off_set);
		
	  for(unsigned int j=0;j<(my_buff_size-2);j++)
	    {
	      rma_buff[off_set*i+j+2].elem_id=Children[i].at(j).elem_id;
	      rma_buff[off_set*i+j+2].c1=Children[i].at(j).c1;
	      rma_buff[off_set*i+j+2].c2=Children[i].at(j).c2;
	    }	
		
	}
	  
    }
		
		
  //rma_buff[35].c1=1.*my_rank;      
  //printf("my_rank=%d %d\n",my_rank,max_buff_size*sizeof(MPI_Message_strct));     
  // the first entry in rma_buff is the destination and the secomnd entry is 
  // the size of the messaga to be sent
          
  for(int i=0;i<6;i++)
    {
      //i=4;
      my_buff_size=Children[i].size()+2;
      displacement1=(i*off_set);
      displacement2=(6*off_set+co_face_no[i]*off_set);
      destination=Dest[i];
			
      MPI_Win_fence(0, win);
      
      //MPI_Win_lock(MPI_LOCK_EXCLUSIVE,my_rank,0,win);  
      //displacement=(my_rank*off_set*0);
     // printf("displacement1=%d displacement1=%d co_face_no =%d my_buff_size=%d\n",displacement1,displacement2,co_face_no[i],my_buff_size);
      //MPI_Win_lock(MPI_LOCK_EXCLUSIVE,my_rank,0,win);  
      if(Children[i].size()!=0)
	{
		//MPI_Win_lock(MPI_LOCK_SHARED,destination,0,win);  	  
		//MPI_Win_lock(MPI_LOCK_SHARED,my_rank,0,win);  	  
	  // printf("inside here myrank=%d buffsize=%d  dest=%d disp1 %d disp2 %d\n",my_rank,my_buff_size,destination,displacement1,displacement2);
	  //		      MPI_Put(rma_buff,my_buff_size,MPI_Message_strct,destination,displacement,my_buff_size,MPI_Message_strct,win);
	  
	  MPI_Put(rma_buff+displacement1,my_buff_size,MPI_Message_strct,destination,displacement2,my_buff_size,MPI_Message_strct,win);
	
	  //MPI_Win_unlock(destination,win);	
	  // MPI_Win_unlock(my_rank,win);			  		  	      
	  rma_buff[displacement1].elem_id=comm_size;
	  
	  //MPI_Put(rma_buff+30,6,MPI_Message_strct,1,30,6,MPI_Message_strct,win);
	  //MPI_Put(rma_buff,6,MPI_Message_strct,1,30,6,MPI_Message_strct,win);
	}
	       
      MPI_Win_fence(0, win);
      //MPI_Win_unlock(my_rank,win);			  		  
    }
    
   msg_end=clock(); 
   
   {
		printf(ANSI_COLOR_RED"my_rank %d time spend after message passing %16.16lf \n" ANSI_COLOR_RESET,my_rank,double(+msg_end-msg_start)/CLOCKS_PER_SEC);
	}
 msg_start=clock();
 
#if(0)
  if(my_rank==0)
    {             
      for(int k=0;k<max_buff_size;k++)
	{
	  printf("%d %d : %d %lf %lf \n",my_rank,k,rma_buff[k].elem_id,rma_buff[k].c1,rma_buff[k].c2);  
	}
    }		
#endif
  for(int i=0;i<6;i++)
    {
      // displacement1=(i*off_set);
      //displacement2=(6*off_set+i*off_set);	   
				
      //printf("my_rank=%d need to update my nbrs j=%d\n",my_rank,co_face_no[i]);  
				
      if(my_rank==rma_buff[6*off_set+off_set*i].elem_id && rma_buff[6*off_set+i*off_set+1].elem_id!=0)
	{
	  //printf("my_rank=%d need to update my nbrs j=%d\n",my_rank,i);  
	  // assign the changes to the neigboring element
	  //void UpdateNonlocalNbrs(int my_rank,Cube &cube, Message_Struct* rma_buff,const int off_set,const int face_tag,const int co_face_tag,Center_coords &XYZ,Vector_Unint &affected_list)
	  
	   for(unsigned int j=0;j<rma_buff[6*off_set+off_set*i+1].elem_id-2;j++)
	   {
		   rcv_list.push_back(Message_Struct());
		   rcv_list.at(j).elem_id=rma_buff[6*off_set+off_set*i+2+j].elem_id;
		   rcv_list.at(j).c1=rma_buff[6*off_set+off_set*i+2+j].c1;
		   rcv_list.at(j).c2=rma_buff[6*off_set+off_set*i+2+j].c2;
		   
	   }
	  
	   
	   UpdateNonlocalNbrs(my_rank,cube,rma_buff,off_set,i,i,XYZ,ancestor_length,affected_list,rcv_list);
	   
	   rcv_list.clear();
	 //  rcv_list.shrink_to_fit();
	   
	}
			    
    }
/*		   
  if(my_rank==0)
    {
      //printf("here---------------------------------------\n");
      for(unsigned int k=0;k<max_buff_size;k++)
	{
	  //printf("%d %d : %u %lf %lf \n",my_rank,k,rma_buff[k].elem_id,rma_buff[k].c1,rma_buff[k].c2);  
	  //printf("%d %d : %u \n",my_rank,k,rma_buff[k].elem_id);  
	}
    }
*/
#endif

msg_end=clock(); 
   
   {
		printf(ANSI_COLOR_RED"my_rank %d time spend for update %16.16lf \n" ANSI_COLOR_RESET,my_rank,double(+msg_end-msg_start)/CLOCKS_PER_SEC);
	}

  /*
    for(int i=0;i<6;i++)
    {
    while(Children[i].size()!=0)
    {
    Children[i].pop_back();
    }	 
    }
  */
   // the coordinate part is the same as it was for serial 
  //======================================================================

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

    // update the kid neighbor, neigbors z=0, y=0 and x=0 dont change
    // update the other 3 neighbors
    // sit in the middle of the element, 
    // 0 >> bottom 1>> top 
    // 2 >> right  3>>left
    // 4 >> back   5>> front
    // some of the neghbirs will be inherited from the parent element
  
    unsigned int elem_id;  
    int idx;
    int a,b;
    
    // empty the container to eliminate unnecessary push back 
    */
  for(int i=0;i<6;i++)
    {
      
	  Children[i].clear();
	  //Children[i].shrink_to_fit();
	
    }
    
  MPI_Free_mem(rma_buff);
  MPI_Win_free(&win);
  MPI_Type_free(&MPI_Message_strct);
  

}

//======================================================================
//======================================================================

//void check_list(Vector_Unint& refine_list,int elem_id, int *T);

//======================================================================

void check_list(Vector_Unint& refine_list,unsigned int elem_id, int *T)
{
  
	
  *T=0;
	
  for(unsigned int i=0;i<refine_list.size();i++)
    {
      if(elem_id==refine_list.at(i))
	{
	  *T=1;
	}
    }	
}

void check_list_p(Vector_Int ((&send_buff)[6]),int face_id,int elem_id, int *T)
{
  
  *T=0;
  
    for(unsigned int i=0;i<send_buff[face_id].size();i++)
    {
    if(elem_id==send_buff[face_id].at(i))
    {
    *T=1;
    }
    }
   	
}

//======================================================================

void check_level(Cube& cube,int elem_id,unsigned int nbr_elem_id,unsigned int face_tag,int *T)
{
		
  *T=0;
  int nbr;
	
  for(unsigned int i=0;i<6;i++)
    {
	
	//int j=cube[elem_id].nbr[i].size();	
      for(unsigned int j=0;j<cube[elem_id].nbr[i].size();j++)
	{
	  
	 
	 
	  nbr=cube[elem_id].nbr[i].at(j);
	  
	 //  printf(ANSI_COLOR_GREEN "i=%d j=%d elem=%d nbr=%d size=%u %d--------------\n" ANSI_COLOR_RESET,i,j,elem_id,nbr_elem_id,cube[elem_id].nbr[i].size(),cube[elem_id].nbr[i].at(0));
	  
	 
	  
	  // if the nbr's of the elem_id, has the same element as a neighbor at the same face
	  //we have to report this to the neighboring process, note that this is local nbr so it has to include positive elements
          //printf("nbr=%d\n",nbr);
          
          if(nbr>-1)
	    {
	      for(unsigned int k=0;k<cube[nbr].nonlocal_nbr.size();k++)
		{
			
		  if(cube[nbr].nonlocal_nbr.at(k).elem_id==nbr_elem_id && cube[nbr].nonlocal_nbr.at(k).face_tag==face_tag)
		    {
		      *T=1;
		      break;
		    }
		}
	    }
	    
			 
	}
    }
	// printf("T=%d face_tag=%d--------------\n",*T,face_tag);	
}


//======================================================================
/*
 * 
 *  this routine updates the list for every processor by enforcing 4:1 
 *  it will check the neighbors including the nonlocal neighbors and if necessary 
 * 	add them to the list of elements that are to be refined
 *  input for this is ordered_refine_list is the new list 
 *  calculated based on the 4:1 balance criteria
 *  the rest of the arguments are all constants
 * 
 */
//======================================================================
//bool comparez(Center_Coords XYZ1, Center_Coords XYZ2) { return XYZ1.z < XYZ2.z; };

void sort_refine_list(const int my_rank,const Cube& cube, Vector_Unint &refine_list, Vector_Unint &ordered_refine_list); 
void CountSameTag(int my_rank,const Cube &cube,const unsigned int face_tag, const unsigned int id, int *bol ,unsigned int *nbr,unsigned int *proc_id);
void AmrNewList(const int my_rank,const MPI_Comm comm,Cube& cube,Vector_Unint& refine_list,Vector_Unint &ordered_refine_list)
{
  // dynamically check and see what happens to the 
  // element if I refine the current element, if level offset is more than one adapt that one
  // as well
  
  for(unsigned int i=0;i<refine_list.size();i++)
  {
	  if(refine_list.at(i)>=cube.size())
	  {
		  printf("Invalid Elem Number in AmrNewRefine elem_id=%u cube_size=%lu\n",refine_list.at(i),cube.size());
		  exit(0);
	  }
  }
			
  unsigned int elem_id,l1;
  int nbr_id;
  Vector_Int adapt_level; 		
  unsigned int list_size;	
  // need this vector to keep track of which neighbor is going to affected by
  //	refining an boundary element at the current processor
  Vector_Int send_buff[6];
  Vector_Int dest_buff[6];
	
  int Comm_size;
  MPI_Comm_size(comm, &Comm_size);

  static const int co_face_no[6]={1,0,3,2,5,4};
	   
 unsigned int comm_size=(unsigned int)Comm_size;
  unsigned int Dest[6]={comm_size,comm_size,comm_size,comm_size,comm_size,comm_size};
  unsigned int proc_id; 	
		
  int a=1,b=1;
  int T,T1;
  unsigned int istart=0;	
  unsigned int iend=refine_list.size();
  //printf("istart %d iend %d\n",istart,iend);

  int count=0;	
  MPI_Win win;
  int buff_size,nbr_to_add,off_set;
  int destination,recv_destination;
  unsigned int nbr_elem_id,face_id;
  unsigned int my_buff_size,max_buff_size;
  int *rma_buffer=NULL;

 // static const int face_index[6]={-1,-2,-3,-4,-5,-6};
  // this part is perfectly fine for interior to each processor, same as serial no change
  unsigned int face_tag;
  int bol,my_b,my_a;
  unsigned int displacement1,displacement2;
  //printf("istart=%d iend=%d\n",istart,iend);
#if(1)
   while(a==1 || b==1)
  {
    a=0;
    b=0;
    my_b=0;  
    my_a=0;
      		
    for(unsigned int i=istart;i<iend;i++)
      {
	elem_id=refine_list.at(i);
	//printf("istart=%d iend=%d\n",istart,iend);

	for(unsigned int j=0;j<6;j++)
	  {
		  face_tag=j;
		  
	    for(unsigned int k=0;k<cube[elem_id].nbr[j].size();k++)
	      {
		nbr_id=cube[elem_id].nbr[j].at(k);
		//face_tag=face_index[j];
		
		//exclude the boundaries 
		//printf("===========%d\n",nbr_id);
		// excludes boundary elements
		//buff_size=0;
		
		// for locals
		if(nbr_id>-1)
		  {
		    //printf("===========%d %d %d\n",nbr_id,cube[nbr_id].level,cube[elem_id].level);
		    if(cube[(unsigned)nbr_id].level<cube[elem_id].level)
		      {
															
			check_list(refine_list,nbr_id,&T);
			// check the original list, if this, set is as -1  
			if(T==0)
			  {
			    refine_list.push_back(nbr_id);
			    //printf("nbr_id=%d\n",nbr_id);
			    my_a=1;
			  }
		      }
		  }
#if(1)
		// this part accounts for the nonlocal to processor neighbors 
		// if element level< nonloca__nbr_level  
		// buggy
		// for nonlocals
							 
		else if(nbr_id<0 && cube[elem_id].nonlocal_nbr.size()!=0)
		  {		    
			  
		    CountSameTag(my_rank,cube,face_tag,elem_id,&bol,&nbr_elem_id,&proc_id);
		    Dest[j]=proc_id;		    
		   
			    
		    if(bol==1) 
		      {				
	//		  	 printf(ANSI_COLOR_GREEN "elem_id =%d face_tag=%d bol=%d\n" ANSI_COLOR_RESET,elem_id,face_tag,bol);
			// check to see if anyone else from the current lement neighbors
			// have the same nonlocal_nbr, this implies that that element would
			// violate the 4:1 balance if not reported to the nbr'ing processor
			
			check_level(cube,elem_id,nbr_elem_id,face_tag,&T1);
			
			//printf("????????my_rank=%d elem_id=%d nbr_elem_id=%d T1=%d size=%d\n",my_rank,elem_id,nbr_elem_id,T1,cube[elem_id].nonlocal_nbr.size());
			//printf("face_tag=%d\n",face_tag);							
			//printf("???????? T1=%d\n",T1);		
			
			if(T1)	
			  {
			    check_list_p(send_buff,face_tag,nbr_elem_id,&T);	
	//		    printf(ANSI_COLOR_BLUE "T=%d nbr_elem_id=%d\n" ANSI_COLOR_RESET,T,nbr_elem_id);
			#if(1)					
			    if(T==0)
			      {				
				  send_buff[j].push_back(nbr_elem_id);
				
			      }	
			  #endif
			  }
			  
			   
		      }
		   						 		
		  }
#endif
/*
if(my_rank==1)
    {
      printf("----------------------------------------------\n");
      for(int i=0;i<6;i++)
	{
	  //printf("i=%d %d\n",i,Children[i].size());
	  //if(Children[i].size()!=0)
	  {
	    for(int j=0;j<send_buff[i].size();j++)
	      {
		printf("i=%d %d\t",i,send_buff[i].at(j)); 
	      }
	    printf("\n"); 
	  }
	}
      printf("----------------------------------------------\n"); 
    }	
    */   
    
      }
	  }	
      }	
    // each processor calculates its maximum size of send_buff size
    my_buff_size=0;
	
    for(unsigned int j=0;j<6;j++)
      {
	if(my_buff_size<send_buff[j].size())
	  {
	    my_buff_size=send_buff[j].size();
	  }
      }
	
    MPI_Allreduce(&my_buff_size,&max_buff_size,1,MPI_INT,MPI_MAX,comm);
		
    // need to know tha maximum size for the window for allocation
    // now that we know the worst case scenario (biggest size of matrix)
    // we can allocate the window and perform RMA
  
    off_set=max_buff_size+2;
    max_buff_size=6*(off_set)*2;	 
    
   //  printf("my_rank=%d max_buff_size=%d\n",my_rank,max_buff_size);
	
    MPI_Alloc_mem(sizeof(int)*max_buff_size, MPI_INFO_NULL,&rma_buffer); 
    MPI_Win_create(rma_buffer,max_buff_size*sizeof(int),sizeof(int),MPI_INFO_NULL,comm,&win);
    
    
    for(unsigned int j=0;j<max_buff_size;j++)
    {
      rma_buffer[j]=0;
    }
      
    // the first entry in rma_buff is the destination and the secomnd entry is 
    // the size of the messaga to be sent
  
    for(int i=0;i<6;i++)
      {
				 		  
	if(send_buff[i].size()!=0)
	  {
		  //printf("********my_rank=%d i=%d \n",my_rank,i);
	    my_buff_size=send_buff[i].size()+2;
	    rma_buffer[i*off_set]=Dest[i];	
	    rma_buffer[i*off_set+1]=my_buff_size;
		
	    for(unsigned int j=0;j<send_buff[i].size();j++)
	      {
		rma_buffer[off_set*i+j+2]=send_buff[i].at(j);
		//printf("********my_rank=%d i=%d element=%d %d index=%d\n",my_rank,i,send_buff[i].at(j),rma_buffer[off_set*i+j+2],off_set*i+j+2);
	      }  
	  }
      }	
/*
     if(my_rank==1)
     {
     for(unsigned int j=0;j<max_buff_size;j++)
    {
      //printf("%d\n",rma_buffer[j]);
    }
 }
*/
    //write the data to destination
    // change the lock to shared ?????
    
    for(int i=0;i<6;i++)
      {
	//i=4;
	
	my_buff_size=send_buff[i].size()+2;
	displacement1=(i*off_set);
	displacement2=(6*off_set+co_face_no[i]*off_set);
	destination=Dest[i];
	//rma_buffer[6*off_set+off_set*i]=destination;
		  
	//MPI_Win_lock(MPI_LOCK_EXCLUSIVE,my_rank,0,win);
		MPI_Win_fence(0, win);  
	if(send_buff[i].size()!=0)
	  {
	//	MPI_Win_lock(MPI_LOCK_SHARED,destination,0,win);  
	//MPI_Win_lock(MPI_LOCK_EXCLUSIVE,destination,0,win);  
		//MPI_Win_lock(MPI_LOCK_SHARED,my_rank,0,win);  
	  
	    MPI_Put(rma_buffer+displacement1,my_buff_size,MPI_INT,destination,displacement2,my_buff_size,MPI_INT,win);
	  
	  //MPI_Win_unlock(destination,win);
	    //rma_buffer[displacement1]=comm_size;
	  }
	 //MPI_Win_unlock(my_rank,win);
	  MPI_Win_fence(0, win);
     	 
      } 
      
    // give an offset so each guy would not right on window and 
    // then send the data before using it !!
    for(int i=0;i<6;i++)
      {
	if(my_rank==rma_buffer[6*off_set+off_set*i] && rma_buffer[6*off_set+i*off_set+1]!=0)
	  {		     
	    for(int j=0;j<rma_buffer[6*off_set+i*off_set+1]-2;j++)
	      {
		nbr_to_add=rma_buffer[6*off_set+off_set*i+2+j];
		
		//printf(ANSI_COLOR_MAGENTA "elem_to_add %d\n" ANSI_COLOR_RESET,nbr_to_add);			  
		check_list(refine_list,nbr_to_add,&T);
		
		#if(1)
		if(T==0)
		  {
		    refine_list.push_back(nbr_to_add);
		    my_b=1;
		    // might have to overwrite a on the nbr'ing processor as well
		  }
		  #endif
	      }
	  }  
      }
		
		#if(1)  
    // need to make every one repeat the while loop  
    MPI_Allreduce(&my_b,&b,1,MPI_INT,MPI_MAX,comm);	
    MPI_Allreduce(&my_a,&a,1,MPI_INT,MPI_MAX,comm);					  
					  
    for(int i=0;i<6;i++)
      {
	while(send_buff[i].size()!=0)
	  {
	    send_buff[i].pop_back();
	  }
	
      }
	  			 
    if(a==1 || b==1)
      {
	istart=iend;
	iend=refine_list.size();  
      }
      #endif	
      MPI_Free_mem(rma_buffer); 
      MPI_Win_free(&win);
  }
  
  /*
  if(my_rank==1)
    {             
      for(unsigned int k=0;k<max_buff_size;k++)
	{
	 // printf("%d %d : %d \n",my_rank,k,rma_buffer[k]);  
	}
    }	
 */
#endif    

	//printf("MYRANK=%d SIZE OF ORIGINAL LIST %lu\n",my_rank,refine_list.size());
    
    sort_refine_list(my_rank,cube,refine_list,ordered_refine_list);
    
    
   /* 
if(my_rank==0)
{
for(int i=0;i<ordered_refine_list.size();i++)
{
	printf("TTTTTTTTTTTTTTTTTTTTT %d size %d \n",ordered_refine_list.at(i),ordered_refine_list.size());
}
}
*/ 
/*    for(unsigned int i=0;i<ordered_refine_list.size();i++)
    {
		
		//printf("my_rank =%d final list %d level =%d\n",my_rank,ordered_refine_list.at(i),cube.at(ordered_refine_list.at(i)).level);
		//printf("my_rank =%d final list %d level =%d\n",my_rank,refine_list.at(i),cube.at(refine_list.at(i)).level);
     }
*/

/*
     while(refine_list.size()!=0)
    {
      refine_list.pop_back();
      //refine_list.pop_back();
    }
*/
 

}
//======================================================================
/*
 * first find max level and min level in the list 
 */
//====================================================================== 
bool compareLevel(cube_data C1, cube_data C2) { return C1.level < C2.level; };

void sort_refine_list(const int my_rank,const Cube& cube, Vector_Unint &refine_list, Vector_Unint &ordered_refine_list)
{
	unsigned int BIG=1e6;
   unsigned int min_level=BIG;
   unsigned int max_level=0;

#if(CPP==0)
  for(unsigned int i=0;i<refine_list.size();i++)
    {
      if(min_level>cube.at(refine_list.at(i)).level)
	{
	  min_level=cube.at(refine_list.at(i)).level;
	}
      if(max_level<cube.at(refine_list.at(i)).level)
	{
	  max_level=cube.at(refine_list.at(i)).level;
	}
    }
#endif
  
  //printf("max_level=%d\n",cube.at(cont-cube.begin()).level);
  
  #if(CPP==1)
  auto cont=std::minmax_element(cube.begin(),cube.end(),&compareLevel);
  min_level=cube.at(cont.first-cube.begin()).level;
  max_level=cube.at(cont.second-cube.begin()).level;
#endif    
  // printf("max_level min_level=%d\n",min_level,max_level);
    
    
    //printf("max_level=%d min_level=%d input size=%d output_size=%d\n",max_level,min_level,refine_list.size(),ordered_refine_list.size());
	
  // sort the elements based on their level, since we need to refine 
  // the lower levels first to avoid generation of non-consistent elements
  // along the waym since adaptation method assumes that 4:1 will be preserved

// watch out for this bug, sorting has a bug 
	unsigned int elem;	
//	int bol;

  for(unsigned int i=min_level;i<=max_level;i++)
    {
      for(unsigned int j=0;j<refine_list.size();j++)
	{
		elem=refine_list.at(j);
	//	if(refine_list.at(j)!=BIG)
		{
	  if(cube.at(refine_list.at(j)).level==i)
	    {
			//printf(ANSI_COLOR_RED"%d================%d j=%d size=%d\n" ANSI_COLOR_RESET,my_rank,refine_list.at(j),j,ordered_refine_list.size());
			//printf("ELEM %d\n",elem);			
	        ordered_refine_list.push_back(elem);
	        //refine_list.at(j)=BIG;
	        refine_list.erase(refine_list.begin()+j);
	        j--;
	    }
	  }		
	}
    }
     
    if(refine_list.size()!=0)
    {
		printf(ANSI_COLOR_RED "Error in Ordering the List\n" ANSI_COLOR_RESET);
		exit(0);
	}
	
  /*
  for(unsigned int i=min_level;i<=max_level;i++)
    {
      for(unsigned int j=0;j<refine_list.size();j++)
	{
		elem=refine_list.at(j);
		if(refine_list.at(j)!=BIG)
		{
	  if(cube.at(refine_list.at(j)).level==i)
	    {
			//printf(ANSI_COLOR_RED"%d================%d j=%d size=%d\n" ANSI_COLOR_RESET,my_rank,refine_list.at(j),j,ordered_refine_list.size());
			//printf("ELEM %d\n",elem);			
	        ordered_refine_list.push_back(elem);
	        refine_list.at(j)=BIG;
	   //     refine_list.erase(refine_list.begin()+j);
	
	    }
	  }		
	}
    }
   */
  
  
  
    	/*
	if(refine_list.size()!=ordered_refine_list.size())
	{
		printf(ANSI_COLOR_RED "the sizes before and after reordering do not match\n" ANSI_COLOR_RESET);
		exit(0);
	}
	*/
	//printf(ANSI_COLOR_RED "old list size %d new list size %d\n" ANSI_COLOR_RESET,refine_list.size(),ordered_refine_list.size());
  // clear the refine list as we only need the ordered list
  // refine_list.clear();
	
	
}
