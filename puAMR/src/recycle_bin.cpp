#if(0)
void Amr_four_to_one_enforce(int my_rank,MPI_Comm comm,Cube& cube,Vector_Int& refine_list,int L,int M, int N)
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
		AMR_Refine_Cube(my_rank,cube,elem_id,L,M,N);
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
#endif

#if(0)
// now look only for processor boundary elems
a=1;
istart=0;	
iend=refine_list.size();

while(a==1)
 {		
	a=0;	
	 for(unsigned int i=istart;i<iend;i++)
	 {
		elem_id=refine_list.at(i);

		for(unsigned int j=0;j<6;j++)
		{
			for(unsigned int k=0;k<cube[elem_id].nbr[j].size();k++)
			{
				nbr_id=cube[elem_id].nbr[j].at(k);
					
				if(nbr_id>-1 && cube[elem_id].nonlocal_nbr.size()!=0)
				{
					// if my level is higher than the neighbor we are fine so only
					// do this if the neighbor has a lower level
					buff_size=0;
					if(cube[elem_id].nonlocal_nbr.size()==0)
					{
				send_buffer[3*l]=cube[elem_id].nonlocal_nbr[l].proc_id;	
				send_buffer[3*l+1]=cube[elem_id].nonlocal_nbr[l].elem_id;	
				send_buffer[3*l+2]=cube[elem_id].nonlocal_nbr[l].face_tag;	
				destination=send_buffer[3*l];
				buffsize=cube[elem_id].nonlocal_nbr.size()*3;
				}
			   }			
				MPI_Bcast(&recv_destination,1,MPI_INT,my_rank,MPI_COMM_WORLD);				
			// send this elem_id, its level to the destination processor
			   				
				if(my_rank==recv_destination)
				{
				MPI_Sendrecv(&send_buff,buff_size,MPI_INT,destination,0,&recv_buff,buff_size,MPI_INT,my_rank,0,comm,&status);
								
				if(cube[(unsigned)nbr_id].level<cube[elem_id].level)
				{
					// refine it
					//
					//refine_list.push_back(nbr_id);
					//printf("nbr_id=%d\n",nbr_id);
					nbr_id=recv_buff[2];
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
		
	}		 
	 if(a==1)
    {
	  istart=iend;
	  iend=refine_list.size();  
    }
 }

//======================================================================


void AMR_Refine_Cube(const int my_rank,const MPI_Comm comm,Cube& cube,const int id,const int L,const int M, const int N)
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

 /* save the newly added tuplet to the structure needed for parallel implementation
   this is in order for the boundary element to know
   on which processor to which element it is connected to
 */
 
   for(i=0;i<cube[id].nonlocal_nbr.size();i++){
     cube_temp[0].nonlocal_nbr.push_back(Nonlocal_Nbr());
     cube_temp[0].nonlocal_nbr[i].proc_id=cube[id].nonlocal_nbr[i].proc_id; 
     cube_temp[0].nonlocal_nbr[i].face_tag=cube[id].nonlocal_nbr[i].face_tag;
     cube_temp[0].nonlocal_nbr[i].elem_id=cube[id].nonlocal_nbr[i].elem_id;
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
 
 //  set adaptation level, adaptation level might be calcutable given the size of the 
 // initial and current size, due to diviion by 2 at each level
 cube[id].level= cube_temp[0].level+1;
 
 for(j=0;j<7;j++)
   {
     cube[nquad+j].level=cube_temp[0].level+1;
     //cube[nquad+j].q=(double*)malloc(q_size*sizeof(double));     
   }
 
 // eliminate the old neighbor data from the parent element since now 
 //parent will be ith child (kind of like each element pointing to itself )
 
 for(i=0;i<6;i++)
   {
     cube[id].nbr[i].clear();
     cube[id].nonlocal_nbr.clear();
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

for(i=1;i<8;i++){
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

// the coordinate part is the same as it was for serial 
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
	 if(cube[cube_temp[0].nbr[i].at(j)].nbr[co_face_no[i]].at(k)==id){		  
	     cube[cube_temp[0].nbr[i].at(j)].nbr[co_face_no[i]].erase(cube[cube_temp[0].nbr[i].at(j)].nbr[co_face_no[i]].begin()+k);
	   }
	 }
       }     
     }
   }
   
   // tell to the neighboring processor to eliminate the current element form its neghborhood connectivity
   int destination;
   int nbr_id;
   int data[2];
   int rcv_data[2];
   MPI_Status status;
   
   for(i=0;i<cube[id].nonlocal_nbr.size();i++)
   {
	destination=cube[id].nonlocal_nbr[i].proc_id;
	data[0]=my_rank;
	data[1]=cube[id].nonlocal_nbr[i].elem_id;   
	
	MPI_Sendrecv(&data,2,MPI_INT,destination,0,&rcv_data,2, MPI_INT,my_rank,0,comm,&status);

	if(my_rank==destination)
	{
		for(j=0;j<cube[data[1]].nonlocal_nbr.size();j++)
	{			
		if(cube[data[1]].nonlocal_nbr[j].proc_id==data[1])
		{
			cube[data[1]].nonlocal_nbr.erase(cube[data[1]].nonlocal_nbr.begin()+j);
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
// Assign interior faces for each element, no change in here for parallel
// implementation 
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
   // same mapping for upper elements only the first index changes
  //  each three entries in I coorespond to the each exposed face of each newly
  // generated element
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
	  //??????????????????????????????????????????????????????????????????????????????????
	  //find_adjacent_face(cube,cube_temp,elem_id,a,b,&idx);
	  // printf("elem_id %d  face id =%d idx=%d\n",elem_id,b,idx);
	  //cube[elem_id].nbr[b].push_back(idx);
	  
	  // update the neighbor if neighbor exists
	  if(idx>-1)
	    {	   
	      cube[idx].nbr[co_face_no[b]].push_back(elem_id);	       
	    }
	      
	}          
    }	  
//======================================================================

void save(Cube& cube,Cube& cube_temp, int id)
{
  unsigned int i,j;
	
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
}
    
    
    
#endif


      #if(0)
     for(int i=0;i<6;i++)
      {
		//  i=4;
		  my_buff_size=Children[i].size()+2;
		 		  
	  if(Children[i].size()!=0)
		{
		//rma_buffer[i*off_set]=dest_buff[i].at(0);	
		
		destination=Dest[i];
		rma_buff[i*off_set].elem_id=destination;
		rma_buff[i*off_set+1].elem_id=my_buff_size;
		
		printf("i=%d Dest=%d my_buff_size=%d %d %d offset=%d\n",i,Dest[i],my_buff_size,Children[i].size(),MPI_Message_strct_size*max_buff_size,i*off_set);
		
		for(int j=0;j<(my_buff_size-2);j++)
		{
			rma_buff[off_set*i+j+2].elem_id=Children[i].at(j).elem_id;
			//rma_buff[off_set*i+j+2].c1=Children[i].at(j).c1;
			//rma_buff[off_set*i+j+2].c2=Children[i].at(j).c2;
		}	
		 
		 //displacement=(co_face_no[i]*off_set);
		  displacement1=(i*off_set);
		  displacement2=(co_face_no[i]*off_set);
		 
		  printf("--------------------------------------------------\n");
		  //printf("%d %lf %lf\n",rma_buff[i*off_set].elem_id,rma_buff[i*off_set+2].c1,rma_buff[i*off_set+2].c2);
		  //write the data to destination
		  // change the lock to shared ?????
		  /*
		  //MPI_Win_fence(0, win);
		  MPI_Win_lock(MPI_LOCK_EXCLUSIVE,my_rank,0,win);
		 
		  //displacement=(my_rank*off_set*0);
		  printf("displacement=%d co_face_no =%d my_buff_size=%d\n",displacement1,co_face_no[i],my_buff_size);
		  		  
		  //if(my_rank==1)
		  {
		  MPI_Put(rma_buff,my_buff_size,MPI_Message_strct,destination,displacement2,my_buff_size,MPI_Message_strct,win);
		  //MPI_Put(rma_buff,6,MPI_Message_strct,1,30,6,MPI_Message_strct,win);
	       }
		  MPI_Win_unlock(my_rank,win);
		  //MPI_Win_fence(0, win);
          */
		  		  
		  // give an offset so each guy would not right on window and 
		  // then send the data before using it !!
		  #if(1)
		  if(my_rank==rma_buff[co_face_no[i]*off_set].elem_id)
		  {
		   face_tag=i;
		   co_face_tag=co_face_no[i];
		   UpdateNonlocalNbrs(cube,rma_buff,off_set,face_tag,co_face_tag,affected_list);
		  }	
		  #endif 
	  }	
	  
	  
	  
  }
 #endif 
 
 
  // this routine is an older version, it is buggy
// will clean this out soon

#if(0)
  //================================================================================
  //
  //                               write Y 
  //
  //================================================================================

  filespace = H5Screate_simple(RANK, dimsf, NULL); 
  memspace  = H5Screate_simple(RANK, chunk_dims, NULL); 
 
  for(i=0;i<cube.size();i++)
    {  
      
      Xa=cube[i].xyz[0];
      Xb=cube[i].xyz[1];
      Ya=cube[i].xyz[2];
      Yb=cube[i].xyz[3];
      Za=cube[i].xyz[4];
      Zb=cube[i].xyz[5];

      hx=L-1.0;
      hy=M-1.0;
      hz=N-1.0;
      Xh=(Xb-Xa)/(hx);
      Yh=(Yb-Ya)/(hy);
      Zh=(Zb-Za)/(hz);

      
      // generate xyz for every block and write to hdf5
	
      index=0;    
      //printf("%lf %lf %lf %lf %lf %lf\n", Xa,Xb,Ya,Yb,Za,Zb);
    
      for(j=0;j<L;j++)
	{
	  for(k=0;k<M;k++)
	    {
	      for(l=0;l<N;l++)
		{
		  ytemp[index]=Ya+Yh*k;  
		  index++;
		}
	    }
	}	
	 		
      // define the offset, only in the fourth dimension

      offset[3] = my_offset+i;	
 
      // write out X coordinates
      sprintf(str0, "/Y");
  
      /*
       * Create chunked dataset.
       */
      plist_id = H5Pcreate(H5P_DATASET_CREATE);
      H5Pset_chunk(plist_id, RANK, chunk_dims);
      dset_id = H5Dcreate(file_id, str0, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
      H5Pclose(plist_id);
      H5Sclose(filespace);

      filespace = H5Dget_space(dset_id);
      status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, block);

      plist_id = H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    
      status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace,filespace,
			plist_id, ytemp);
 
      H5Dclose(dset_id);
      H5Sclose(filespace);
      H5Pclose(plist_id);
      H5Sclose(memspace);

    }


  //================================================================================
  //
  //                               write Z 
  //
  //================================================================================


  filespace = H5Screate_simple(RANK, dimsf, NULL); 
  memspace  = H5Screate_simple(RANK, chunk_dims, NULL); 
 
  for(i=0;i<cube.size();i++)
    {  
      
      Xa=cube[i].xyz[0];
      Xb=cube[i].xyz[1];
      Ya=cube[i].xyz[2];
      Yb=cube[i].xyz[3];
      Za=cube[i].xyz[4];
      Zb=cube[i].xyz[5];

      hx=L-1.0;
      hy=M-1.0;
      hz=N-1.0;
      Xh=(Xb-Xa)/(hx);
      Yh=(Yb-Ya)/(hy);
      Zh=(Zb-Za)/(hz);

      
      // generate xyz for every block and write to hdf5
	
      index=0;    
      //printf("%lf %lf %lf %lf %lf %lf\n", Xa,Xb,Ya,Yb,Za,Zb);
    
      for(j=0;j<L;j++)
	{
	  for(k=0;k<M;k++)
	    {
	      for(l=0;l<N;l++)
		{
		  ztemp[index]=Za+Zh*l;
		  index++;
		}
	    }
	}	
	 		
      // define the offset, only in the fourth dimension

      offset[3] = my_offset+i;	
 
      // write out X coordinates
      sprintf(str0, "/Z");
  
      /*
       * Create chunked dataset.
       */
      plist_id = H5Pcreate(H5P_DATASET_CREATE);
      H5Pset_chunk(plist_id, RANK, chunk_dims);
      dset_id = H5Dcreate(file_id, str0, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
      H5Pclose(plist_id);
      H5Sclose(filespace);

      filespace = H5Dget_space(dset_id);
      status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, block);

      plist_id = H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    
      status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace,filespace,plist_id, ztemp);
 
      H5Dclose(dset_id);
      H5Sclose(filespace);
      H5Pclose(plist_id);
      H5Sclose(memspace);

    }

#endif







/*

#if(0)
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
#endif



/*
 *  MPI_Bcast(&recv_destination,1,MPI_INT,my_rank,MPI_COMM_WORLD);
			 
 if(my_rank==recv_destination)
 {
 MPI_Sendrecv(&send_buff,buff_size,MPI_INT,destination,0,&recv_buff,buff_size,MPI_INT,my_rank,0,comm,&status);
								
 if(cube[(unsigned)nbr_id].level<cube[elem_id].level)
 {
 // refine it
 //
 //refine_list.push_back(nbr_id);
 //printf("nbr_id=%d\n",nbr_id);
 nbr_id=recv_buff[2];
 check_list(refine_list,nbr_id,&T);
 // check the original list, if this, set is as -1  
 if(T==0)
 {
 refine_list.push_back(nbr_id);
 //printf("nbr_id=%d\n",nbr_id);
 a=1;
 }
 }
								
 }*/
/*
 * this routine checks to see if there is a 4:1 connectivity before we start refining
 * 
 * 
*/


#if(0)
#include <stdlib.h>
#include "zoltan.h"

typedef struct{
  int numGlobalPoints;
  int numMyPoints;
  ZOLTAN_ID_PTR myGlobalIDs;
  double *c;
} MESH_DATA;

static int get_number_of_objects(void *data, int *ierr);

static int get_num_geometry(void *data, int *ierr);

static void get_object_list(void *data, int sizeGID, int sizeLID,
            ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                  int wgt_dim, float *obj_wgts, int *ierr);
                  
static void get_geometry_list(void *data, int sizeGID, int sizeLID,
                      int num_obj,
             ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
             int num_dim, double *geom_vec, int *ierr);
int Zoltan_Partition(int argcs, char* pArgs[],int my_rank,MPI_Comm Comm,const Cube &cube)
//int main(int argcs, char* pArgs[])
{

float ver; 


  MPI_Init(&argcs,&pArgs);
    
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  int np;
  int rc;
  MPI_Comm_size(MPI_COMM_WORLD,&np); 
  

  MPI_Comm cworld;
  MPI_Comm_dup(MPI_COMM_WORLD, &cworld);
  

 rc = Zoltan_Initialize(argcs, pArgs, &ver);

  if (rc != ZOLTAN_OK){
    printf("sorry...\n");
    MPI_Finalize();
    exit(0);
  }
 
 struct Zoltan_Struct *zz=NULL;

 zz = Zoltan_Create(MPI_COMM_WORLD);

  /* General parameters */

#if(1)
  Zoltan_Set_Param(zz, "DEBUG_LEVEL", "1");
  //Zoltan_Set_Param(zz, "LB_METHOD", "HSFC");
    //Zoltan_Set_Param(zz, "LB_METHOD", "RCB");
    Zoltan_Set_Param(zz, "LB_METHOD", "BLOCK");
  Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1"); 
  Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");
  //Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "0");
  Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");


	MESH_DATA myMesh;
//	myMesh=(MESH_DATA*)malloc(1.0*sizeof(MESH_DATA));
  	
  	myMesh.numGlobalPoints=9;
   	myMesh.numMyPoints=2;
	myMesh.myGlobalIDs=NULL;
	myMesh.myGlobalIDs=(ZOLTAN_ID_PTR)malloc(2*sizeof(ZOLTAN_ID_PTR));
	
	myMesh.myGlobalIDs[0]=2*my_rank;
	myMesh.myGlobalIDs[1]=2*my_rank+1;
	
	if(my_rank==2)
	{
	myMesh.numMyPoints=5;
	myMesh.myGlobalIDs=NULL;
	myMesh.myGlobalIDs=(ZOLTAN_ID_PTR)malloc(5*sizeof(ZOLTAN_ID_PTR));
	myMesh.myGlobalIDs[0]=4;
	myMesh.myGlobalIDs[1]=5;
	myMesh.myGlobalIDs[2]=6;
	myMesh.myGlobalIDs[3]=7;
	myMesh.myGlobalIDs[4]=8;
	}
   
   myMesh.c=NULL;
   myMesh.c=(double*)malloc(myMesh.numMyPoints*sizeof(double));
   
   for(int i=0;i<myMesh.numMyPoints;i++)
   {
	myMesh.c[i]=1.*myMesh.myGlobalIDs[i]; 
	printf("my_rank %d %lf\n",my_rank,myMesh.c[i]);  
   }
   
   // printf("my_rank =%d number of my points %d glob_id global_id %d %d \n",my_rank,myMesh->numMyPoints,myMesh->myGlobalIDs[0],myMesh->myGlobalIDs[1]);
  
  
  /* Query functions, to provide geometry to Zoltan */

    Zoltan_Set_Num_Obj_Fn(zz, get_number_of_objects, &myMesh);

    Zoltan_Set_Obj_List_Fn(zz, get_object_list, &myMesh);
    
    Zoltan_Set_Num_Geom_Fn(zz, get_num_geometry, &myMesh);
    
    Zoltan_Set_Geom_Multi_Fn(zz, get_geometry_list, &myMesh);
    
    int changes, numGidEntries, numLidEntries, numImport, numExport;
  ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids; 
  int *importProcs, *importToPart, *exportProcs, *exportToPart;
  int *parts;
  
  
  
	rc = Zoltan_LB_Partition(zz, /* input (all remaining fields are output) */
        &changes,        /* 1 if partitioning was changed, 0 otherwise */ 
        &numGidEntries,  /* Number of integers used for a global ID */
        &numLidEntries,  /* Number of integers used for a local ID */
        &numImport,      /* Number of vertices to be sent to me */
        &importGlobalGids,  /* Global IDs of vertices to be sent to me */
        &importLocalGids,   /* Local IDs of vertices to be sent to me */
        &importProcs,    /* Process rank for source of each incoming vertex */
        &importToPart,   /* New partition for each incoming vertex */
        &numExport,      /* Number of vertices I must send to other processes*/
        &exportGlobalGids,  /* Global IDs of the vertices I must send */
        &exportLocalGids,   /* Local IDs of the vertices I must send */
        &exportProcs,    /* Process to which I send each of the vertices */
        &exportToPart);  /* Partition to which each vertex will belong */
#endif

printf("my_rank %d changes %d numExport %d numImport %d\n",my_rank,changes,numExport,numImport);




Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids, 
                      &importProcs, &importToPart);
  Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids, 
                      &exportProcs, &exportToPart);

  Zoltan_Destroy(&zz);

  /**********************
  ** all done ***********
  **********************/

  MPI_Finalize();

return(0);


}

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

   for (int i=0; i<mesh->numMyPoints; i++)
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
  return 2;
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
    geom_vec[2*i] = (double)mesh->c[i];
    geom_vec[2*i + 1] = 0.0;
  }

  return;
}



#endif


