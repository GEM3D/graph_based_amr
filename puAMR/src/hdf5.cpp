//#include "/home/jhasbestan/GEM_AMR/src/include/typedefs.h"
#include "typedefs.h"
#include "hdf5.h"

//======================================================================
//======================================================================
//
//                    Parallel Version Using Chuncks
//
//======================================================================
//======================================================================


// here rank means dimension of the entire block
#define RANK 4
#define H5FILE_NAME  "soln/Pxdmf3d.h5"
#define XDMF_NAME    "soln/Pxdmf3d.xmf" 
/* 
 * two communicators one stands for world and the other for cartesian
 * 
 */
 
 void XdmfParallelSpatialCollection(Cube &cube,int offset,MPI_Comm comm,MPI_Info info,unsigned int L,unsigned int M,unsigned int N,int ncube_total);

void WriteHdf5ParallelSpatialCollection(Cube& cube,int my_rank,unsigned int L,unsigned int M,unsigned int N, int npx,int npy,int npz,MPI_Comm comm,MPI_Comm com,MPI_Info info,int my_offset,int ncube_total,Center_coords &XYZ, double ancestor_length[3])
{
   
  hid_t       file_id, dset_id;         /* file and dataset identifiers */
  hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
  hsize_t     dimsf[4];                 /* dataset dimensions */
  hsize_t     chunk_dims[4];            /* chunk dimensions */
  hsize_t	count[4];	          /* hyperslab selection parameters */
  hsize_t	block[4];
  hsize_t	offset[4];
  hid_t	plist_id;                 /* property list identifier */
  int         i,j,k,l;
  herr_t	status;
  //int         *data=NULL;    

  /*
   * MPI variables
   */
  int mpi_size, mpi_rank;
  double *xtemp =NULL; 
  double *ytemp =NULL;
  double *ztemp =NULL;
  double *qtemp =NULL;	
  	
  xtemp=new double[L*M*N];
  ytemp=new double[L*M*N];
  ztemp=new double[L*M*N];
  qtemp=new double[L*M*N];
  
  double Xa;
  double Xb;
  double Ya;
  double Yb;
  double Za;
  double Zb;

  double hx;
  double hy;
  double hz;
  double Xh;
  double Yh;
  double Zh;
	
  unsigned int index;

  char str0[50];
  char str1[50];
  char str2[50];
  char str3[50];


  /*
   * Initialize MPI
   */
  
  //MPI_Comm_rank(comm, &mpi_rank);  
  

  /* 
   * Set up file access property list with parallel I/O access
   */

  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, comm, info);

  /*
   * Create a new file collectively and release property list identifier.
   */

  file_id = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);
 

  /*
   * Create the dataspace for the dataset.
   */

   
  dimsf[0] = L;
  dimsf[1] = M;
  dimsf[2] = N;
  dimsf[3] = ncube_total;	

  chunk_dims[0] = L;   
  chunk_dims[1] = M;   
  chunk_dims[2] = N;
  chunk_dims[3] = 1;
  
  /* 
   * Each process defines dataset in memory and writes it to the hyperslab
   * in the file.
   */

  count[0] = 1;
  count[1] = 1 ;
  count[2] = 1;	
  count[3] = 1;	
 
  block[0] = chunk_dims[0];
  block[1] = chunk_dims[1];
  block[2] = chunk_dims[2];		
  block[3] = chunk_dims[3];		

  /*
    for loop such that each processor can write all the blocks
  */ 
  offset[0] = 0;
  offset[1] = 0;
  offset[2] = 0;
  //==============================save coordinates to generate grid==============

 
 
  //================================================================================
  //
  //                               write X
  //
  //================================================================================
 
  filespace = H5Screate_simple(RANK, dimsf, NULL); 
  memspace  = H5Screate_simple(RANK, chunk_dims, NULL); 

  sprintf(str0, "/X");
  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  
  H5Pset_chunk(plist_id, RANK, chunk_dims);
 
  dset_id = H5Dcreate(file_id, str0, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
 #if(1)
 
  H5Pclose(plist_id);
  H5Sclose(filespace);
  filespace = H5Dget_space(dset_id);

  unsigned int level,coord_index;	

  double dx,dy,dz;
  double denum;
	
  for(unsigned int i=0;i<cube.size();i++)
    {  	

	  level=cube.at(i).level;
	  coord_index=cube.at(i).centeroid_index;
	  denum=pow(2,level);
	  dx=ancestor_length[0]/denum;
	  Xa=XYZ.at(coord_index).x-dx*0.5;
      Xb=XYZ.at(coord_index).x+dx*0.5;     
      hx=L-1.0;
      Xh=(Xb-Xa)/(hx);
      
      
      index=0;    
      //printf("%lf %lf %lf %lf %lf %lf\n", Xa,Xb,Ya,Yb,Za,Zb);
      // printf("???????%d %d\n", my_rank,cube.size());
    
      for(unsigned int j=0;j<L;j++)
	{
	  for(unsigned int k=0;k<M;k++)
	    {
	      for(unsigned int l=0;l<N;l++)
		{
		  xtemp[index]=Xa+Xh*j;	      
		  index++;
		}
	    }
	}		 		
      // define the offset, only in the fourth dimension
  
      offset[3] = my_offset+i;	
      // printf("my_rank=%d offset=%d\n",my_rank,offset[3]);
      //status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, block);
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, block);
      plist_id = H5Pcreate(H5P_DATASET_XFER);
      //H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
      H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
      //status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace,filespace,plist_id, xtemp);
      H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace,filespace,plist_id, xtemp);
      H5Pclose(plist_id);
     
      
    }
  //need to close dset_id since the new variable is going to be Y and Z ...
  // no need to open and close memspace and filespace all the time, simply reuse them 

  H5Dclose(dset_id);
  
  //H5Sclose(filespace);
  //H5Sclose(memspace);

 //================================================================================
  //
  //                               write Y 
  //
  //================================================================================

  //filespace = H5Screate_simple(RANK, dimsf, NULL); 
  //memspace  = H5Screate_simple(RANK, chunk_dims, NULL); 

  sprintf(str0, "/Y");
  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(plist_id, RANK, chunk_dims);
  dset_id = H5Dcreate(file_id, str0, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  H5Sclose(filespace);
  filespace = H5Dget_space(dset_id);

  for(unsigned i=0;i<cube.size();i++)
    {  	

      level=cube.at(i).level;
	  coord_index=cube.at(i).centeroid_index;
	  denum=pow(2,level);
	  dy=ancestor_length[1]/denum;
      Ya=XYZ.at(coord_index).y-dy*0.5;
      Yb=XYZ.at(coord_index).y+dy*0.5;
      hy=M-1.0;
      Yh=(Yb-Ya)/(hy);
      index=0;    
      //printf("%lf %lf %lf %lf %lf %lf\n", Xa,Xb,Ya,Yb,Za,Zb);
      // printf("???????%d %d\n", my_rank,cube.size());
    	 
      for(unsigned int j=0;j<L;j++)
	{
	  for(unsigned int k=0;k<M;k++)
	    {
	      for(unsigned int l=0;l<N;l++)
		{
		  ytemp[index]=Ya+Yh*k;  
		  index++;
		}
	    }
	}	
	 	
	 		
      // define the offset, only in the fourth dimension
  
      offset[3] = my_offset+i;	
      //printf("my_rank=%d offset=%d\n",my_rank,offset[3]);
    
      //status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, block);
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, block);
      plist_id = H5Pcreate(H5P_DATASET_XFER);
      //H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
      H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
      //status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace,filespace,plist_id, ytemp);
      H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace,filespace,plist_id, ytemp);
      H5Pclose(plist_id);
     
      
    }

  H5Dclose(dset_id);
 
  // H5Sclose(filespace);
  //H5Sclose(memspace);

 //================================================================================
    //                               write Z 
  //
  //================================================================================

  //filespace = H5Screate_simple(RANK, dimsf, NULL); 
 // memspace  = H5Screate_simple(RANK, chunk_dims, NULL); 

  sprintf(str0, "/Z");
  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(plist_id, RANK, chunk_dims);
  dset_id = H5Dcreate(file_id, str0, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  H5Sclose(filespace);
  filespace = H5Dget_space(dset_id);

  for(unsigned int i=0;i<cube.size();i++)
    {  	

      level=cube.at(i).level;
	  coord_index=cube.at(i).centeroid_index;
	  denum=pow(2,level);
	  dz=ancestor_length[2]/denum;
      Za=XYZ.at(coord_index).z-dz*0.5;
      Zb=XYZ.at(coord_index).z+dz*0.5;      
      hz=N-1.0;
      Zh=(Zb-Za)/(hz);
      index=0;    
      	 
       for(unsigned int j=0;j<L;j++)
	{
	  for(unsigned int k=0;k<M;k++)
	    {
	      for(unsigned int l=0;l<N;l++)
		{
		  ztemp[index]=Za+Zh*l;
		  index++;
		}
	    }
	}	
	 	 		
 // define the offset, only in the fourth dimension
  
      offset[3] = my_offset+i;	
      //printf("my_rank=%d offset=%d\n",my_rank,offset[3]);
    
      //status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, block);
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, block);
      plist_id = H5Pcreate(H5P_DATASET_XFER);
      //H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
	  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);	
      //status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace,filespace,plist_id, ztemp);
      H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace,filespace,plist_id, ztemp);
      H5Pclose(plist_id);
     
      
    }
H5Dclose(dset_id);

//======================================================================
//
//                                write q
//
//======================================================================
// for now the q vector only holds one value for all teh corners 
// this has to modified for flow solver case
#if(1)
  sprintf(str0, "/Q");
  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(plist_id, RANK, chunk_dims);
  dset_id = H5Dcreate(file_id, str0, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  H5Sclose(filespace);
  filespace = H5Dget_space(dset_id);

  for(unsigned int i=0;i<cube.size();i++)
    {  	

     
      index=0;    
      	 
       for(unsigned int j=0;j<L;j++)
	{
	  for(unsigned int k=0;k<M;k++)
	    {
	      for(unsigned int l=0;l<N;l++)
		{
		  qtemp[index]=cube[i].q;
		  index++;
		}
	    }
	}	
	 	 		
 // define the offset, only in the fourth dimension
  
      offset[3] = my_offset+i;	
      //printf("my_rank=%d offset=%d\n",my_rank,offset[3]);
    
      //status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, block);
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, block);
      plist_id = H5Pcreate(H5P_DATASET_XFER);
      //H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
	  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);	
      //status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace,filespace,plist_id, qtemp);
      H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace,filespace,plist_id, qtemp);
      H5Pclose(plist_id);
     
      
    }
H5Dclose(dset_id);


#endif
 
  
//======================================================================  
  H5Sclose(filespace);
  H5Sclose(memspace);
  

#endif
 

  H5Fclose(file_id);
  
  delete[] xtemp;
  delete[] ytemp;
  delete[] ztemp;
  delete[] qtemp;
  
  XdmfParallelSpatialCollection(cube,my_offset,comm,info,L,M,N,ncube_total);
  
}


//========================================================
//
//           write xdmf in parallel
//
//=======================================================
//
//
// 
// This fujnction writes the meta data required by visit or paraview to visualize hdf5 
// this function is hard coded for now, watch that those offset match the ones you use 
//
//

#define offset0 156
//#define offset1 1292 //for L=5
#define offset1 1861
      //#define offset1 1500

//===========================================================

void integer_string(char *strin,int i)
{

      if(i<10)
	{
      sprintf(strin,"00000%d",i);
    }
      else if(i<100)
	{
      sprintf(strin,"0000%d",i);
    }
      else if(i<1000)
	{
      sprintf(strin,"000%d",i);
    }

      else if(i<10000)
	{
      sprintf(strin,"00%d",i);
    }

      else if(i<100000)
	{
      sprintf(strin,"0%d",i);
    }
      else if(i<1000000)
	{
      sprintf(strin,"%d",i);
    }
      else
	{
      printf("not in the range, offset is too big\n");
      exit(0);
    }

    }

//======================================================================
// This MPI_File_Iwrite gives funcky errors, use write for now
#if(1)
void XdmfParallelSpatialCollection(Cube &cube,int offset,MPI_Comm comm,MPI_Info info,unsigned int L,unsigned int M,unsigned int N,int ncube_total)
{
  MPI_File fp; 
  int buf[1000], my_rank,np;  
  MPI_Comm_rank(comm, &my_rank); 
  MPI_File_open(comm, XDMF_NAME, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
  MPI_Comm_size(comm,&np); 
  MPI_Request request; 
  unsigned int i;
  int j;

  MPI_Status status;
  const char* names[] = {"X", "Y", "Z"};
  
  char strL[100];
  char strM[100];
  char strN[100];
  char strNcube[1000];
  char stroff[1000];
  int index;

  if(L<10)
    {
      sprintf(strL,"00%d",L);
      sprintf(strM,"00%d",M);
      sprintf(strN,"00%d",N);
    }
  else if(L<100)
    {
      sprintf(strL,"0%d",L);
      sprintf(strM,"0%d",M);
      sprintf(strN,"0%d",N);
    }
  else
    {
      printf("discretization too big go change xdmf function\n");
      exit(0);
    }


  if(ncube_total<10)
    {
      sprintf(strNcube,"00000%d",ncube_total);
    }
  else if(ncube_total<100)
    {
      sprintf(strNcube,"0000%d",ncube_total);
    }
  else if(ncube_total<1000)
    {
      sprintf(strNcube,"000%d",ncube_total);
    }
  else if(ncube_total<10000)
    {
      sprintf(strNcube,"00%d",ncube_total);
    }
  else if(ncube_total<100000)
    {
      sprintf(strNcube,"0%d",ncube_total);
    }
  else if(ncube_total<1000000)
    {
      sprintf(strNcube,"%d",ncube_total);
    }
  else
    {
      printf("number of cubes larger than 1000000\n");
      exit(0);
    }

 
 char str[1000];
 // a counts the offset for header whic is only written by process rank 0 and  and b the hyperslab part for each cube 
 int a=0,b=0;
  if(my_rank==0)
    {      
      sprintf(str,"<?xml version=\"1.0\" ?>\n");
      MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);     
      a=a+strlen(str);
      sprintf(str,"<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []> \n");
      MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);
      a=a+strlen(str);
      sprintf(str,"<Xdmf Version=\"2.0\">\n");
      MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);
      a=a+strlen(str);
      sprintf(str,"<Domain>\n");
      MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);
       a=a+strlen(str);
      sprintf(str,"<Grid Name=\"AMR\" GridType=\"Collection\" CollectionType=\"Spatial\">\n");
      MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);
       a=a+strlen(str);
      printf("size a=%d\n",a);
      
      for(unsigned int i=0;i<cube.size();i++)
	{
	  // b calculates the offset for each block, need to be set to zero becasue it is inside the loop
	  // but we dont want to accumulate, calculate only for one block
	  b=0;
	  sprintf(str,"   <Grid Name=\"mesh0\" GridType=\"Uniform\">\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	
	  b=b+strlen(str);  
	  sprintf(str,"        <Topology TopologyType=\"3DSMesh\" NumberOfElements=\"%s %s %s\"/>\n",strL,strM,strN);
	  // printf("hyperslab =%d\n",strlen(str));
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
	  b=b+strlen(str);  
	  sprintf(str,"          <Geometry GeometryType=\"X_Y_Z\">  \n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);
	  b=b+strlen(str);  
	  
	  for(j=0;j<3;j++)
	    {
	  sprintf(str,"      <DataItem ItemType=\"HyperSlab\" Dimensions=\"%s %s %s %d\" NumberType=\"Float\" Precision=\"4\"  Type=\"HyperSlab\"> \n",strL,strM,strN,1);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);
	  b=b+strlen(str);  	 
	  sprintf(str,"   	<DataItem Dimensions=\"3 4\" Format=\"XML\" > \n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
	  b=b+strlen(str);  
	  index=offset+i;
	  integer_string(stroff,index);
	  sprintf(str,"   	%d %d %d %s  \n",0,0,0,stroff);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
	  b=b+strlen(str);  
	  sprintf(str,"   	%d %d %d %d  \n",1,1,1,1);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
	  b=b+strlen(str);  
	  sprintf(str,"   	%s %s %s %d  \n",strL,strM,strN,1);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
	  b=b+strlen(str);  
	  sprintf(str,"   	</DataItem>\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	
	  b=b+strlen(str);   
	  sprintf(str,"   	 <DataItem Name=\"%s\" Dimensions=\"%s %s %s %s\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",names[j],strL,strM,strN,strNcube);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
	  b=b+strlen(str);  
	  sprintf(str,"   	 Pxdmf3d.h5:/%s\n",names[j]);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	
	  b=b+strlen(str);  
	  sprintf(str,"         </DataItem>\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
	  b=b+strlen(str);  
	  sprintf(str,"         </DataItem>\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	  
	  b=b+strlen(str);  
	    }
	  sprintf(str,"      </Geometry>   \n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
	  b=b+strlen(str);  
	  // added because of including solution vector Q
	  sprintf(str,"         <Attribute Name=\"Q\" AttributeType=\"Scalar\" Center=\"Node\"> \n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);
	  b=b+strlen(str);
	  sprintf(str,"      <DataItem ItemType=\"HyperSlab\" Dimensions=\"%s %s %s %d\" NumberType=\"Float\" Precision=\"4\"  Type=\"HyperSlab\"> \n",strL,strM,strN,1);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);
	  b=b+strlen(str);
	  sprintf(str,"   	<DataItem Dimensions=\"3 4\" Format=\"XML\" > \n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);
	  b=b+strlen(str);
	  index=offset+i;
	  integer_string(stroff,index);	 
	  sprintf(str,"   	%d %d %d %s  \n",0,0,0,stroff);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
	  b=b+strlen(str);
	  sprintf(str,"   	%d %d %d %d  \n",1,1,1,1);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
	  b=b+strlen(str);
	  sprintf(str,"   	%s %s %s %d  \n",strL,strM,strN,1);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
	  b=b+strlen(str);
	  sprintf(str,"   	</DataItem>\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
	  b=b+strlen(str);
	  sprintf(str,"   	 <DataItem Name=\"Q\" Dimensions=\"%s %s %s %s\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",strL,strM,strN,strNcube);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
	  b=b+strlen(str);
	  sprintf(str,"   	 Pxdmf3d.h5:/Q\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	
	  b=b+strlen(str);
	  sprintf(str,"         </DataItem>\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
	  b=b+strlen(str);
	  sprintf(str,"         </DataItem>\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	  
	  b=b+strlen(str);
	  sprintf(str," </Attribute>\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	
	  b=b+strlen(str);	  
	  sprintf(str,"  </Grid>\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	  
	  b=b+strlen(str);  
	  
	  // put error here 
	  if(a!=offset0 || b!=offset1)
	    {
	      printf("??????????????????????????????????????????\n");
	      printf("Go fix your offset for Xdmf meta data file a=%d offset0=%d b=%d offset1=%d\n",a,offset0,b,offset1);
	      printf("??????????????????????????????????????????\n");
	        exit(0);
	    }
	  // printf("size length9=%d length11=%d\n",strlen("9"),strlen("11"));

	}
      // b=ftell(&fp);
      //    printf("size b=%d\n",b);
      
    }
  
  else{
    
  // other ranks do something else
    // 5 is added because of the header of the file written by  processor 0

      int mpi_offset=0;
      
     for(i=0;i<cube.size();i++)
	{

       mpi_offset=((offset+i)*offset1+offset0);

       mpi_offset=mpi_offset*sizeof(char);
	  //printf("my_rank = %d offset =%d mpi_offset=%d i=%d\n",my_rank,offset,mpi_offset,i);

	   MPI_File_seek(fp,mpi_offset,MPI_SEEK_SET);
	  // MPI_File_set_view(fp,mpi_offset,MPI_CHAR,MPI_CHAR, "native", MPI_INFO_NULL );
	  sprintf(str,"   <Grid Name=\"mesh0\" GridType=\"Uniform\">\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	  
	  sprintf(str,"        <Topology TopologyType=\"3DSMesh\" NumberOfElements=\"%s %s %s\"/>\n",strL,strM,strN);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
	  sprintf(str,"          <Geometry GeometryType=\"X_Y_Z\">  \n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);

	  for(j=0;j<3;j++)
	    {
	  sprintf(str,"      <DataItem ItemType=\"HyperSlab\" Dimensions=\"%s %s %s %d\" NumberType=\"Float\" Precision=\"4\"  Type=\"HyperSlab\"> \n",strL,strM,strN,1);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
	  sprintf(str,"   	<DataItem Dimensions=\"3 4\" Format=\"XML\" > \n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);
	  index=offset+i;
	  integer_string(stroff,index);	 
	  sprintf(str,"   	%d %d %d %s  \n",0,0,0,stroff);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
	  sprintf(str,"   	%d %d %d %d  \n",1,1,1,1);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
	  sprintf(str,"   	%s %s %s %d  \n",strL,strM,strN,1);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
	  sprintf(str,"   	</DataItem>\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
	  sprintf(str,"   	 <DataItem Name=\"%s\" Dimensions=\"%s %s %s %s\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",names[j],strL,strM,strN,strNcube);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
	  sprintf(str,"   	 Pxdmf3d.h5:/%s\n",names[j]);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	
	  sprintf(str,"         </DataItem>\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
	  sprintf(str,"         </DataItem>\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	  
	    }
	  sprintf(str,"      </Geometry>   \n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
	  // insert here the value of qs
	  
	  sprintf(str,"         <Attribute Name=\"Q\" AttributeType=\"Scalar\" Center=\"Node\"> \n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);
	  sprintf(str,"      <DataItem ItemType=\"HyperSlab\" Dimensions=\"%s %s %s %d\" NumberType=\"Float\" Precision=\"4\"  Type=\"HyperSlab\"> \n",strL,strM,strN,1);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);
	  sprintf(str,"   	<DataItem Dimensions=\"3 4\" Format=\"XML\" > \n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);
	  index=offset+i;
	  integer_string(stroff,index);	 
	  sprintf(str,"   	%d %d %d %s  \n",0,0,0,stroff);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
	  sprintf(str,"   	%d %d %d %d  \n",1,1,1,1);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
	  sprintf(str,"   	%s %s %s %d  \n",strL,strM,strN,1);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
	  sprintf(str,"   	</DataItem>\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
	  sprintf(str,"   	 <DataItem Name=\"Q\" Dimensions=\"%s %s %s %s\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",strL,strM,strN,strNcube);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
	  sprintf(str,"   	 Pxdmf3d.h5:/Q\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	
	  sprintf(str,"         </DataItem>\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
	  sprintf(str,"         </DataItem>\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	  
	  sprintf(str," </Attribute>\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	
	  
	  sprintf(str,"  </Grid>\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	  	  
    }
     // add the last part of the xdmf text file to close the arguments
     if(my_rank==(np-1))
       {
      sprintf(str,"  </Grid>\n");
      MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	  
      sprintf(str,"      </Domain>   \n");
      MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
      sprintf(str,"  </Xdmf>\n");
      MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	  
 
    }
    
  }


      MPI_File_close(&fp);

}

#endif

//=======================================================================================
// this function calculates the number of cubes before each process and hence 
// the offset required to write the hdf5
//
//========================================================================================
#if(0)
void calculate_processor_offset_for_IO(int my_rank,int mycube,int RMA,int np,int *offset)
{

  int ncube=0;
  int i;
  MPI_Status status;
  int off_set;
  int source, destination;

  if(RMA==0)
    {
      if(my_rank>0)
	{
	  source=my_rank-1;	
	  //printf("mryand=%d source %d\n",my_rank, source);
	  MPI_Recv(&ncube,1, MPI_INT,source,0,MPI_COMM_WORLD,&status);
	  printf("mryand=%d ncube recieved %d\n",my_rank, ncube);
	  off_set=ncube;
	  (*offset)=off_set;
	}

      if(my_rank<(np-1))
	{
	  destination=my_rank+1;
	  //printf("mryank=%d destination %d\n",my_rank, destination);
	  ncube=ncube+mycube;
	  //printf("mryank=%d ncube %d mycube %d\n",my_rank,ncube,mycube);
	  MPI_Send(&ncube,1, MPI_INT,destination,0,MPI_COMM_WORLD);
	}
    }
  else if(RMA==1)
    {

      MPI_Win win;
      off_set=mycube;
      MPI_Win_create(&off_set,sizeof(int),sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD,&win); 
      MPI_Win_fence(0,win);
      // printf("np=%d\n",np); 
 
      for(i=0;i<np-1;i++)
	{
	  if(my_rank==i)
	    {
	      destination=my_rank+1;
	      MPI_Accumulate(&off_set,1, MPI_INT,destination,0,1,MPI_INT, MPI_SUM, win);
	
	    }
	  MPI_Win_fence(0,win); 
	}
      (*offset)=off_set-mycube;
      MPI_Win_free(&win);
    }
  else
    {
      printf("RMA input is a boolean set it to 1 for true and 0 for false\n");
    }
}
#endif
//==========================================================

void calculate_processor_offset_for_IO(int my_rank,int mycube,int RMA,int np,int *offset,int *ranks, int* off1)
{

  int ncube=0;
  int i,j;
  MPI_Win win;
  MPI_Status status;
  MPI_Group group,sub_group;
  int group_size;
  int off_set;
  int *off=NULL;

  int source, destination;
  if(RMA==0)
    {
//	this is not correct for a general topology
      if(my_rank>0)
	{
	  source=my_rank-1;	
	  //printf("mryand=%d source %d\n",my_rank, source);
	  MPI_Recv(&ncube,1, MPI_INT,source,0,MPI_COMM_WORLD,&status);
	  printf("mryand=%d ncube recieved %d\n",my_rank, ncube);
	  off_set=ncube;
	  (*offset)=off_set;
	}

      if(my_rank<(np-1))
	{
	  destination=my_rank+1;
	  //printf("mryank=%d destination %d\n",my_rank, destination);
	  ncube=ncube+mycube;
	  //printf("mryank=%d ncube %d mycube %d\n",my_rank,ncube,mycube);
	  MPI_Send(&ncube,1, MPI_INT,destination,0,MPI_COMM_WORLD);
	}
 
    }
  else if(RMA==1)
    {
      MPI_Alloc_mem(sizeof(int), MPI_INFO_NULL,&off); 
      off[0]=mycube;
      MPI_Win_create(off,sizeof(int),sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD,&win); 
      MPI_Win_fence(0,win);
      // printf("np=%d\n",np); 
 
      for(i=0;i<np-1;i++)
	{
	  if(my_rank==i)
	    {
	      destination=my_rank+1;
	      MPI_Accumulate(off,1, MPI_INT,destination,0,1,MPI_INT, MPI_SUM, win);
	    }
	  MPI_Win_fence(0,win); 
	}
 
      (*offset)=off[0]-mycube;
      MPI_Free_mem(off);
      MPI_Win_free(&win);
    }
// this section deadlocks on kestrel
#if(0)
  else if(RMA==2)
    {

      MPI_Alloc_mem(sizeof(int), MPI_INFO_NULL,&off); 
      MPI_Win_create(off,sizeof(int),sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD,&win); 
      off[0]=mycube;

      for(i=0;i<np-1;i++)
	{ 
	  MPI_Win_get_group(win,&group);
	  group_size=0;
    
	  for(j=0;j<np;j++)
	    {
	      if(j>=i)
		{
		  ranks[group_size]=j;
		  group_size++;
		}
	    }
          
	  MPI_Group_incl(group,group_size,ranks,&sub_group);
	  MPI_Group_free(&group);

	  if(my_rank>=i)
	    {
	      MPI_Win_start(sub_group,0,win);   
	      MPI_Win_post(sub_group,0,win);   
	    } 
	  if(my_rank==i)
	    {	 
	      destination=my_rank+1;	
	      MPI_Accumulate(off,1, MPI_INT,destination,0,1,MPI_INT, MPI_SUM, win);
	    }
    
	  if(my_rank>=i)
	    {
	      MPI_Win_complete(win);
	      MPI_Win_wait(win);
	    }
	}

      (*offset)=off[0]-mycube;
      MPI_Group_free(&sub_group);
      MPI_Free_mem(off); 
      MPI_Win_free(&win);

    }
#endif
  else
    {
      printf("RMA input is a boolean set it to 1 for true and 0 for false\n");
    }

}

//===================================================================================================
// only need to lock and unlock at the source and destination so lock and unlock are perfect for this

void Broadcast(int my_rank,int mycube,int RMA,int np,int offset,int *ranks,int *off1, int* ncube_total)
{
  int ncube=0;
  int i,j;
  MPI_Win win;
  MPI_Status status;
  MPI_Group group,sub_group;
  int group_size;
  int cube_size;
  int source, destination;      
  int dispatcher;
  dispatcher=np-1;
   
    
  if(RMA==0)
    {
      (*ncube_total)=(offset)+mycube;
      MPI_Bcast(ncube_total,1,MPI_INT,dispatcher,MPI_COMM_WORLD);
    }
  else if(RMA==1)
    {
	  int *off=NULL; 	
      MPI_Alloc_mem(sizeof(int), MPI_INFO_NULL,&off);  
      off[0]=offset+mycube;
      // 
  
      MPI_Win_create(off,sizeof(int),sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD,&win);  
      //MPI_Win_fence(0,win);

      if(my_rank!=(np-1))
	{
	  // lock the access to the target process
	   MPI_Win_lock(MPI_LOCK_EXCLUSIVE,dispatcher,0,win); 
	  //  printf("myrank%d?????????????\n",my_rank);
	  MPI_Get(off,1,MPI_INT,dispatcher,0,1, MPI_INT,win);
	  MPI_Win_unlock(dispatcher,win); 
	}
      // MPI_Win_complete(win);
      //MPI_Win_wait(win);
      
      //MPI_Win_fence(0,win);
      //MPI_Win_fence(0,win);

      (*ncube_total)=off[0];
      
            
      printf("inside my_rank=%d total number of cubes=%d\n",my_rank,off[0]);
      MPI_Win_fence(0,win);
      MPI_Free_mem(off); 
      MPI_Win_free(&win); 
    }
  else
    {
      printf("RMA input is a boolean set it to 1 for true and 0 for false\n");
    }

}

// this routine wries the centeroid of cubes only, no connectivity




void XdmfParallelPolyvertex(const Cube &cube,const int offset,const MPI_Comm comm,const MPI_Info info,const int ncube_total);

void WriteHdf5ParallelPolyvertex(const Cube& cube,const int my_rank,const int npx,const int npy,const int npz,const MPI_Comm comm,const MPI_Info info,const int my_offset,const int ncube_total,const Center_coords &XYZ)
{
   
  hid_t       file_id, dset_id;         /* file and dataset identifiers */
  hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
  hsize_t	count;	          /* hyperslab selection parameters */
  hsize_t	block;
  hsize_t	offset;
  hid_t	plist_id;                 /* property list identifier */
  int         i,j,k,l;
  herr_t	status;
  //int         *data=NULL;    
  hsize_t total_size=ncube_total;
  /*
   * MPI variables
   */
  int mpi_size, mpi_rank;
  double *xtemp =NULL; 
  double *ytemp =NULL;
  double *ztemp =NULL;
  double *qtemp =NULL;
 
  
  xtemp=new double[cube.size()];
  ytemp=new double[cube.size()];
  ztemp=new double[cube.size()];   
  qtemp=new double[cube.size()];
  

  double Xc;
  double Yc;  
  double Zc;
 
  char str0[50];
  char str1[50];
  char str2[50];
  char str3[50];

  /* 
   * Set up file access property list with parallel I/O access
   */


  plist_id = H5Pcreate(H5P_FILE_ACCESS);

  H5Pset_fapl_mpio(plist_id, comm, info);
#if(1)
  /*
   * Create a new file collectively and release property list identifier.
   */

  file_id = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);
 
#endif
  /*
   * Create the dataspace for the dataset.
   */

   
  
  /* 
   * Each process defines dataset in memory and writes it to the hyperslab
   * in the file.
   */

  //==============================save coordinates to generate grid==============

  //================================================================================
  //
  //                               write X
  //
  //================================================================================
   
  count=cube.size();
  offset=my_offset;
  block=1;
     
  filespace = H5Screate_simple(1, &total_size, NULL); 
  memspace  = H5Screate_simple(1, &count, NULL); 


 //================================================================================
  //
  //                               write X 
  //
  //================================================================================
  sprintf(str0, "/X");
  plist_id = H5Pcreate(H5P_DATASET_CREATE);  
  dset_id = H5Dcreate(file_id, str0, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  H5Sclose(filespace);
  
  filespace = H5Dget_space(dset_id);

  unsigned int coord_index;
	
  for(unsigned int i=0;i<cube.size();i++)
    {  	

	  coord_index=cube.at(i).centeroid_index;
	  Xc=XYZ.at(coord_index).x;
        
      xtemp[i]=Xc;
     }
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset, NULL, &count, &block);
      plist_id = H5Pcreate(H5P_DATASET_XFER);
      //H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
      H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
      //status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace,filespace,plist_id, xtemp);
      H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace,filespace,plist_id, xtemp);
      H5Pclose(plist_id);
      
      
 //================================================================================
  //
  //                               write Y 
  //
  //================================================================================
            
  sprintf(str0, "/Y");
  plist_id = H5Pcreate(H5P_DATASET_CREATE);  
  dset_id = H5Dcreate(file_id, str0, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  H5Sclose(filespace);
  
  filespace = H5Dget_space(dset_id);

  
  for(unsigned int i=0;i<cube.size();i++)
    {  	
	  coord_index=cube.at(i).centeroid_index;
      Yc=XYZ.at(coord_index).y;
      ytemp[i]=Yc;
}  
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset, NULL, &count, &block);
      plist_id = H5Pcreate(H5P_DATASET_XFER);
      //H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
      H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
      //status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace,filespace,plist_id, xtemp);
      H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace,filespace,plist_id, ytemp);
      H5Pclose(plist_id);
      
      
 //================================================================================
  //
  //                               write Z
  //
  //================================================================================
      
  sprintf(str0, "/Z");
  plist_id = H5Pcreate(H5P_DATASET_CREATE);  
  dset_id = H5Dcreate(file_id, str0, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  H5Sclose(filespace);
  
  filespace = H5Dget_space(dset_id);

  
  for(unsigned int i=0;i<cube.size();i++)
    {  	
	  coord_index=cube.at(i).centeroid_index;
      Zc=XYZ.at(coord_index).z;
      ztemp[i]=Zc;
}  
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset, NULL, &count, &block);
      plist_id = H5Pcreate(H5P_DATASET_XFER);
      //H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
      H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
      //status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace,filespace,plist_id, xtemp);
      H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace,filespace,plist_id, ztemp);
      H5Pclose(plist_id);
               
      
// ask ray to see what is wrong about closing 
 // H5Fclose(file_id);
 
 //================================================================================
  //
  //                               write X 
  //
  //================================================================================
  sprintf(str0, "/Q");
  plist_id = H5Pcreate(H5P_DATASET_CREATE);  
  dset_id = H5Dcreate(file_id, str0, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  H5Sclose(filespace);
  
  filespace = H5Dget_space(dset_id);

  for(unsigned int i=0;i<cube.size();i++)
    {  	
	   qtemp[i]=cube.at(i).q;
     }
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset, NULL, &count, &block);
      plist_id = H5Pcreate(H5P_DATASET_XFER);
      //H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
      H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
      //status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace,filespace,plist_id, xtemp);
      H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace,filespace,plist_id, qtemp);
      H5Pclose(plist_id);
       
    
      H5Dclose(dset_id);
      H5Sclose(filespace);
      H5Sclose(memspace);
  
  delete[] xtemp;
  delete[] ytemp;
  delete[] ztemp;
  delete[] qtemp;
 
 XdmfParallelPolyvertex(cube,my_offset,comm,info,ncube_total);
 
}

//======================================================================
//
//
//======================================================================

void XdmfParallelPolyvertex(const Cube &cube,const int offset,const MPI_Comm comm,const MPI_Info info,const int ncube_total)
{
  MPI_File fp; 
  int buf[1000], my_rank,np;  
  MPI_Comm_rank(comm, &my_rank); 
  MPI_File_open(comm,XDMF_NAME, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
  MPI_Comm_size(comm,&np); 
  MPI_Request request; 
  unsigned int i;
  int j;

  MPI_Status status;
  const char* names[] = {"X", "Y", "Z"};
  
  char strL[100];
  char strM[100];
  char strN[100];
  char strNcube[1000];
  char stroff[1000];
  int index;

  char str[1000];
 // a counts the offset for header whic is only written by process rank 0 and  and b the hyperslab part for each cube 
 int a=0,b=0;
  if(my_rank==0)
    {      
      sprintf(str,"<?xml version=\"1.0\" ?>\n");
      MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);     
      sprintf(str,"<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []> \n");
      MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);
      sprintf(str,"<Xdmf Version=\"2.0\">\n");
      MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);      
      sprintf(str,"<Domain>\n");
      MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);    
      sprintf(str,"   <Grid Name=\"mesh0\" GridType=\"Uniform\">\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);		   
	  sprintf(str,"        <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"%d\"/>\n",ncube_total);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);		   
	  sprintf(str,"          <Geometry GeometryType=\"X_Y_Z\">  \n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);
	  b=b+strlen(str);  	  
	  sprintf(str,"   	 <DataItem Name=\"X\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",ncube_total);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	
	  sprintf(str,"   	 Pxdmf3d.h5:/%s\n",names[0]);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status); 
	  sprintf(str,"         </DataItem>\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	
	  sprintf(str,"   	 <DataItem Name=\"Y\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",ncube_total);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	
	  sprintf(str,"   	 Pxdmf3d.h5:/%s\n",names[1]);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);  
	  sprintf(str,"         </DataItem>\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	
	  sprintf(str,"   	 <DataItem Name=\"Z\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",ncube_total);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);
	  sprintf(str,"   	 Pxdmf3d.h5:/%s\n",names[2]);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status); 	 
	  sprintf(str,"         </DataItem>\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	
	  sprintf(str,"      </Geometry>   \n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	
      sprintf(str,"         <Attribute Name=\"Q\" AttributeType=\"Scalar\" Center=\"Node\"> \n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);
	  sprintf(str,"   	 <DataItem Name=\"Q\" Dimensions=\"%d \" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",ncube_total);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
	  sprintf(str,"   	 Pxdmf3d.h5:/Q\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	
	  sprintf(str,"         </DataItem>\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
	  sprintf(str," </Attribute>\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	
	  sprintf(str,"  </Grid>\n");  
      MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
      sprintf(str,"      </Domain>   \n");
      MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
      sprintf(str,"  </Xdmf>\n");
      MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	  	  
      
    }
  
  

      MPI_File_close(&fp);

}

//======================================================================
//
//                   Unstructured Way of Showing
//
//======================================================================

void NodeCoords(double *xyz,double *Xa, double *Xb,double *Ya, double *Yb, double *Za, double *Zb,int index);
bool comparex(Center_Coords XYZ1, Center_Coords XYZ2) { return XYZ1.x < XYZ2.x; };
bool comparey(Center_Coords XYZ1, Center_Coords XYZ2) { return XYZ1.y < XYZ2.y; };
bool comparez(Center_Coords XYZ1, Center_Coords XYZ2) { return XYZ1.z < XYZ2.z; };
void IsInList(const Center_coords &XYZ,const double *xyz, int *bol);
void XdmfParallelUnstructured(const Cube &cube,const int offset,const MPI_Comm comm,const MPI_Info info,const int ncube_total,const int nnode_total);
void WriteHdf5ParallelUnstructured(const Cube& cube,const int my_rank,const int npx,const int npy,const int npz,const MPI_Comm comm,const MPI_Info info,const int my_offset,const int ncube_total,const Center_coords &XYZ,const double ancestor_length[3],Center_coords &XYZ2)
{
   
  hid_t       file_id, dset_id;         /* file and dataset identifiers */
  hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
  hsize_t	count;	          /* hyperslab selection parameters */
  hsize_t	block;
  hsize_t	offset;
  hid_t	plist_id;                 /* property list identifier */
  int         i,j,k,l;
  herr_t	status;
  //int         *data=NULL;    
  hsize_t total_size=ncube_total;
  hsize_t total_node_size; 
  /*
   * MPI variables
   */
  int mpi_size, mpi_rank;
  double *xtemp =NULL; 
 /* double *ytemp =NULL;
  double *ztemp =NULL;
  */
  double *qtemp =NULL;
 
  	
  	
  char str0[50];
  char str1[50];
  char str2[50];
  char str3[50];
 double idenum;
 int np=npx*npy*npz;
  /*
   * Initialize MPI
   */
  
  //MPI_Comm_rank(comm, &mpi_rank); 
 
	
  int bol;
  /*
  Center_coords XYZ2;
  XYZ2.reserve(4*cube.size());	 
  */
  double xyz[3];
  double dx,dy,dz;
  double Xa,Xb,Ya,Yb,Za,Zb,denum;
   unsigned int level,coord_index,current_size;	
   
   printf(ANSI_COLOR_GREEN"Inside HDF5\n"  ANSI_COLOR_RESET);
   
   
    
   
  for(unsigned int i=0;i<cube.size();i++)
    {  	

      level=cube.at(i).level;
	  coord_index=cube.at(i).centeroid_index;
	  TwoPowN(level,&idenum);
	  idenum=1./idenum;
	  dx=ancestor_length[0]*idenum;
	  Xa=XYZ.at(coord_index).x-dx*0.5;
      Xb=XYZ.at(coord_index).x+dx*0.5;     
      
      //denum=pow(2,level);
	  dy=ancestor_length[1]*idenum;
      Ya=XYZ.at(coord_index).y-dy*0.5;
      Yb=XYZ.at(coord_index).y+dy*0.5; 
      
      //denum=pow(2,level);
	  dz=ancestor_length[2]*idenum;
      Za=XYZ.at(coord_index).z-dz*0.5;
      Zb=XYZ.at(coord_index).z+dz*0.5;    
      
      
      for(j=0;j<8;j++)
      {
		NodeCoords(xyz,&Xa,&Xb,&Ya,&Yb,&Za,&Zb,j);
	    IsInList(XYZ2,xyz,&bol);
	    
	    if(bol)
	    {
			current_size=XYZ2.size();
			XYZ2.push_back(Center_Coords());
			
			XYZ2.at(current_size).x=xyz[0];
			XYZ2.at(current_size).y=xyz[1];
			XYZ2.at(current_size).z=xyz[2];
			
		}
	}
     }
 printf(ANSI_COLOR_GREEN"AFTER LOOP\n"  ANSI_COLOR_RESET);
/*
if(my_rank==0)
{
	for(i=0;i<XYZ2.size();i++)
	{
		printf("%lf %lf %lf \n",XYZ2.at(i).x,XYZ2.at(i).y,XYZ2.at(i).z);
	}
}
* */
std::sort(XYZ2.begin(), XYZ2.end(), &comparex);

//std::sort(XYZ2.begin(), XYZ2.end(), &comparey);

  /*
if(my_rank==1)
{
	for(i=0;i<XYZ2.size();i++)
	{
		printf("%lf %lf %lf \n",XYZ2.at(i).x,XYZ2.at(i).y,XYZ2.at(i).z);
	}
}  
  
 
if(my_rank==1)
{
Center_Coords t;
t.x=1.0;
 
 //std::vector<Center_Coords>::difference_type cont1
 
 auto cont1=std::lower_bound(XYZ2.begin() , XYZ2.end() ,t,&compare);;
 
 //cont1=std::lower_bound(XYZ2.begin() , XYZ2.end() ,t,&compare);
	   
	   //std::distance(XYZ2.begin(), cont)
  printf("cont %lu\n",cont1-XYZ2.begin());
  
}
*/
  /*
   * Create the dataspace for the dataset.
   */

   
  
   
    MPI_Barrier(MPI_COMM_WORLD);
  /* 
   * Each process defines dataset in memory and writes it to the hyperslab
   * in the file.
   */

  //==============================Nodebased data communication==============

 int mynode=XYZ2.size();
  
 printf(ANSI_COLOR_GREEN"Start of communication\n"  ANSI_COLOR_RESET);
  
 int *ranks=NULL;
  ranks=(int*)malloc(np*sizeof(int));
 int *off1=NULL;
 int my_offset2=0,nnode_total=0;
 calculate_processor_offset_for_IO(my_rank,mynode,1,np,&my_offset2,ranks,off1);
 Broadcast(my_rank,mynode,0,np,my_offset2,ranks,off1,&nnode_total);
 printf(ANSI_COLOR_GREEN"my_rank=%d total number of nodes=%d nodaloffset=%d\n"  ANSI_COLOR_RESET,my_rank,nnode_total,my_offset2);
 MPI_Free_mem(off1);
 free(ranks);

 printf(ANSI_COLOR_GREEN"End of communication\n"  ANSI_COLOR_RESET);
  //================================================================================
  //
  //                               write X
  //
  //================================================================================
 
// use one variable for all to save memory 

  xtemp=new double[mynode];
  //ytemp=new double[mynode];
  //ztemp=new double[mynode];   
  
  
  for(unsigned int i=0;i<XYZ2.size();i++)
  {
	  xtemp[i]=XYZ2.at(i).x;
	//  ytemp[i]=XYZ2.at(i).y;
	//  ztemp[i]=XYZ2.at(i).z;
	//  xtemp[3*i+1]=XYZ2.at(i).y;
	//  xtemp[3*i+2]=XYZ2.at(i).z;
  }
  
  /*
  
  */
  
  /*
  Center_coords XYZ2_temp;
  XYZ2.swap(XYZ2_temp);
  */
  /*
  XYZ2.clear();
  XYZ2.shrink_to_fit();
 
 /* 
  Center_coords XYZ2_temp;
  XYZ2.swap(XYZ2_temp);
  */
    plist_id = H5Pcreate(H5P_FILE_ACCESS);

  H5Pset_fapl_mpio(plist_id, comm, info);
#if(1)
  /*
   * Create a new file collectively and release property list identifier.
   */

  file_id = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);
 
#endif
  
  #if(1)
  
 

//=====================================================================
//
//                              Write Q
//
//=====================================================================

 count=cube.size();
 offset=my_offset;
  block=1;
  
 filespace = H5Screate_simple(1, &total_size, NULL); 
  memspace  = H5Screate_simple(1, &count, NULL);

qtemp=new double[cube.size()];

sprintf(str0, "/Q");
  plist_id = H5Pcreate(H5P_DATASET_CREATE);  
  dset_id = H5Dcreate(file_id, str0, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  H5Sclose(filespace);
  
  filespace = H5Dget_space(dset_id);

  for(unsigned int i=0;i<cube.size();i++)
    {  	
	   qtemp[i]=cube.at(i).q;
     }
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset, NULL, &count, &block);
      plist_id = H5Pcreate(H5P_DATASET_XFER);
      //H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
      H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
      //status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace,filespace,plist_id, xtemp);
      H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace,filespace,plist_id, qtemp);
      H5Pclose(plist_id);
      
 
      H5Dclose(dset_id);
      H5Sclose(filespace);
      H5Sclose(memspace);


delete[] qtemp;
  
 //================================================================================
  //
  //                               write X 
  //
  //================================================================================
  
  offset=my_offset2;
  block=1;
  total_node_size=nnode_total;   
  //count=XYZ2.size();
  count=mynode;
  filespace = H5Screate_simple(1, &total_node_size, NULL); 
  memspace  = H5Screate_simple(1, &count, NULL); 
  
  sprintf(str0, "/X");
  plist_id = H5Pcreate(H5P_DATASET_CREATE);  
  dset_id = H5Dcreate(file_id, str0, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  H5Sclose(filespace);
  
  filespace = H5Dget_space(dset_id);
  
     
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset, NULL, &count, &block);
      plist_id = H5Pcreate(H5P_DATASET_XFER);
      //H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
      H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
      //status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace,filespace,plist_id, xtemp);
      H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace,filespace,plist_id, xtemp);
      H5Pclose(plist_id);
      
      
      for(unsigned int i=0;i<XYZ2.size();i++)
  {
	  
	  xtemp[i]=XYZ2.at(i).y;
	
  }
  
    
      
sprintf(str0, "/Y");
  plist_id = H5Pcreate(H5P_DATASET_CREATE);  
  dset_id = H5Dcreate(file_id, str0, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  H5Sclose(filespace);
  
  filespace = H5Dget_space(dset_id);
  
     
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset, NULL, &count, &block);
      plist_id = H5Pcreate(H5P_DATASET_XFER);
      //H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
      H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
      //status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace,filespace,plist_id, xtemp);
      H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace,filespace,plist_id, xtemp);
      H5Pclose(plist_id);
      
      
      
      for(unsigned int i=0;i<XYZ2.size();i++)
  {
	 
	
	  xtemp[i]=XYZ2.at(i).z;
	
  }
  
 
sprintf(str0, "/Z");
  plist_id = H5Pcreate(H5P_DATASET_CREATE);  
  dset_id = H5Dcreate(file_id, str0, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  H5Sclose(filespace);
  
  filespace = H5Dget_space(dset_id);
  
     
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset, NULL, &count, &block);
      plist_id = H5Pcreate(H5P_DATASET_XFER);
      //H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
      H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
      //status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace,filespace,plist_id, xtemp);
      H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace,filespace,plist_id, xtemp);
      H5Pclose(plist_id);
      
      
// ask ray to see what is wrong about closing 
 // H5Fclose(file_id);
    
      H5Dclose(dset_id);
      H5Sclose(filespace);
      H5Sclose(memspace);
  #endif
  
  
  //====================================================================
  //
  //                        Elem2Node_Connectivity
  //
  //====================================================================
    /*
    int **connectivity=NULL;
    connectivity=(int**)malloc(cube.size()*sizeof(int*));
    // 8 since every cube has eight vertex
    for(unsigned int i=0;i<cube.size();i++)
    {
		connectivity[i]=(int*)malloc(8*sizeof(int));
	}
    */
    
    
    int *connectivity=NULL;   
    //connectivity=(int*)malloc(cube.size()*8*sizeof(int));
    connectivity=new int[cube.size()*8];
    
    
    Center_Coords val;
 
  for(unsigned int i=0;i<cube.size();i++)
    {  	

      level=cube.at(i).level;
	  coord_index=cube.at(i).centeroid_index;
	  //idenum=1./pow(2,level);
	  TwoPowN(level,&idenum);
	  idenum=1./idenum;
	  
	  dx=ancestor_length[0]*idenum;
	  Xa=XYZ.at(coord_index).x-dx*0.5;
      Xb=XYZ.at(coord_index).x+dx*0.5;     
      
      dy=ancestor_length[1]*idenum;
      Ya=XYZ.at(coord_index).y-dy*0.5;
      Yb=XYZ.at(coord_index).y+dy*0.5; 
      
      
	  dz=ancestor_length[2]*idenum;
      Za=XYZ.at(coord_index).z-dz*0.5;
      Zb=XYZ.at(coord_index).z+dz*0.5;   
      
      for(j=0;j<8;j++)
      {
		NodeCoords(xyz,&Xa,&Xb,&Ya,&Yb,&Za,&Zb,j);
	    
	   val.x=xyz[0];
	    
	   #if(CPP==1)
	    
	   auto cont=std::lower_bound(XYZ2.begin() , XYZ2.end() ,val,&comparex);
	   
	   //std::distance(XYZ2.begin(), cont)
	   //printf("cont %lu\n",cont-XYZ2.begin());
	    
	    for(int k=(cont-XYZ2.begin());k<mynode;k++)
	   {
	   //if(xyz[0]==xtemp[k] && xyz[1]==ytemp[k] && xyz[2]==ztemp[k])
	   if(xyz[0]==XYZ2[k].x && xyz[1]==XYZ2[k].y && xyz[2]==XYZ2[k].z)
	   {
	    //connectivity[i][j]=k+my_offset2;
	    connectivity[8*i+j]=k+my_offset2;
	    break;
        }
       }
       #endif
       
       #if(CPP==0)
       for(int k=0;k<mynode;k++)
	   {
	   //if(xyz[0]==xtemp[k] && xyz[1]==ytemp[k] && xyz[2]==ztemp[k])
	   if(xyz[0]==XYZ2[k].x && xyz[1]==XYZ2[k].y && xyz[2]==XYZ2[k].z)
	   {
	    //connectivity[i][j]=k+my_offset2;
	    connectivity[8*i+j]=k+my_offset2;
	    break;
        }
       }
       
       #endif     
	}
  }
  
  
  /*
  for(unsigned int i=0;i<cube.size();i++)
  {
	  for(unsigned int j=0;j<8;j++)
	  {
		    //printf("%d ",connectivity[i][j]);
		   // connectivity[i][j]=1;
		   // printf("%d ",connectivity[i][j]);
		    printf("%d ",connectivity[8*i+j]);
	  }
	  printf("\n");
  }
  */
  //====================================================================
  //
  //                 Write Connectivity to hdf5
  //
  //====================================================================
                  /* hyperslab selection parameters */
  
  
  hsize_t	count2D[2];	         
  hsize_t	block2D[2];
  hsize_t	offset2D[2];
  hsize_t total_element_size2D[2];
    
  
  count2D[0]=cube.size();
  count2D[1]=8;
  total_element_size2D[0]=ncube_total;
  total_element_size2D[1]=8;
  filespace = H5Screate_simple(2, total_element_size2D, NULL); 
  memspace  = H5Screate_simple(2, count2D, NULL); 


  offset2D[0]=my_offset;
  offset2D[1]=0;
  block2D[0]=1;
  block2D[1]=1;
  //count2D[1]=0;
  
  /*printf(ANSI_COLOR_RED"my_rank=%d %d %d\n"ANSI_COLOR_RESET,my_rank,offset2D[0],offset2D[1]);
  printf(ANSI_COLOR_CYAN"my_rank=%d %d %d\n"ANSI_COLOR_RESET,my_rank,total_element_size2D[0],total_element_size2D[1]);
  printf(ANSI_COLOR_MAGENTA"my_rank=%d %d %d\n"ANSI_COLOR_RESET,my_rank,count2D[0],count2D[1]);
 */
 //================================================================================
  //
  //                               write X 
  //
  //================================================================================
  sprintf(str0, "/Hex");
  plist_id = H5Pcreate(H5P_DATASET_CREATE);  
  dset_id = H5Dcreate(file_id, str0, H5T_NATIVE_INT, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  H5Sclose(filespace);
  
  filespace = H5Dget_space(dset_id);
  
     
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset2D, NULL, count2D, block2D);
    
      plist_id = H5Pcreate(H5P_DATASET_XFER);
      //H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
      
      H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
      //status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace,filespace,plist_id, xtemp);
        
      H5Dwrite(dset_id, H5T_NATIVE_INT, memspace,filespace,plist_id,connectivity);
      H5Pclose(plist_id);
  
  
  
  
   H5Sclose(filespace);
   H5Sclose(memspace);
  
  

  delete[] xtemp;
  //delete[] ytemp;
  //delete[] ztemp;
  
  delete[] connectivity;
  
  
  
  XdmfParallelUnstructured(cube,offset,comm,info,ncube_total,nnode_total);
  
 
}

void NodeCoords(double *xyz,double *Xa, double *Xb,double *Ya, double *Yb, double *Za, double *Zb,int index)
{
	
	switch(index)
	{
		case 0:
	    xyz[0]=*Xa;
        xyz[1]=*Ya;
        xyz[2]=*Za;
		break;
		case 1:
	    xyz[0]=*Xb;
        xyz[1]=*Ya;
        xyz[2]=*Za;
		break;
		case 2:
	    xyz[0]=*Xb;
        xyz[1]=*Yb;
        xyz[2]=*Za;
		break;
		case 3:
	    xyz[0]=*Xa;
        xyz[1]=*Yb;
        xyz[2]=*Za;
		break;
		case 4:
	    xyz[0]=*Xa;
        xyz[1]=*Ya;
        xyz[2]=*Zb;
		break;
		case 5:
	    xyz[0]=*Xb;
        xyz[1]=*Ya;
        xyz[2]=*Zb;
		break;
		case 6:
	    xyz[0]=*Xb;
        xyz[1]=*Yb;
        xyz[2]=*Zb;
		break;
		case 7:
	    xyz[0]=*Xa;
        xyz[1]=*Yb;
        xyz[2]=*Zb;
		break;
		
		
	}
	
	
	
	
}


// is in list can be replaced by binary search

void IsInList(const Center_coords &XYZ,const double *xyz, int *bol)
{
	*bol=1;
	
	for(unsigned int i=0;i<XYZ.size();i++)
	{
		if(XYZ.at(i).x==xyz[0] && XYZ.at(i).y==xyz[1] && XYZ.at(i).z==xyz[2])
		{
			*bol=0;
			break;
		}
	}
	
	
}

//======================================================================


void XdmfParallelUnstructured(const Cube &cube,const int offset,const MPI_Comm comm,const MPI_Info info,const int ncube_total,const int nnode_total)
{
  MPI_File fp; 
  int buf[1000], my_rank,np;  
  MPI_Comm_rank(comm, &my_rank); 
  MPI_File_open(comm,XDMF_NAME, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
  MPI_Comm_size(comm,&np); 
  MPI_Request request; 
  unsigned int i;
  int j;

  MPI_Status status;
  const char* names[] = {"X", "Y", "Z"};
  
  char strL[100];
  char strM[100];
  char strN[100];
  char strNcube[1000];
  char stroff[1000];
  int index;

  char str[1000];
 // a counts the offset for header whic is only written by process rank 0 and  and b the hyperslab part for each cube 
 int a=0,b=0;
  if(my_rank==0)
    {      
      sprintf(str,"<?xml version=\"1.0\" ?>\n");
      MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);     
     
      sprintf(str,"<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []> \n");
      MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);
     
      sprintf(str,"<Xdmf Version=\"2.0\">\n");
      MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);
     
      sprintf(str,"<Domain>\n");
      MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);
     
      sprintf(str,"   <Grid Name=\"mesh0\" GridType=\"Uniform\">\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	
	  b=b+strlen(str);  
	  sprintf(str,"        <Topology TopologyType=\"Hexahedron\" NumberOfElements=\"%d\">\n",ncube_total);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
	  sprintf(str,"   	 <DataItem Name=\"Hex\" Dimensions=\"%d %d\" NumberType=\"Int\" Format=\"HDF\">\n",ncube_total,8);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	
	   sprintf(str,"   	 Pxdmf3d.h5:/Hex\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status); 
	   sprintf(str,"         </DataItem>\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	  
	  sprintf(str,"      </Topology>\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);
	  
	  sprintf(str,"          <Geometry GeometryType=\"X_Y_Z\">  \n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);
	  b=b+strlen(str);    
	  sprintf(str,"   	 <DataItem Name=\"X\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nnode_total);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	
	  sprintf(str,"   	 Pxdmf3d.h5:/%s\n",names[0]);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status); 
	  sprintf(str,"         </DataItem>\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	
	  sprintf(str,"   	 <DataItem Name=\"Y\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nnode_total);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	
	  sprintf(str,"   	 Pxdmf3d.h5:/%s\n",names[1]);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);  
	  sprintf(str,"         </DataItem>\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	
	  sprintf(str,"   	 <DataItem Name=\"Z\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nnode_total);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);
	  sprintf(str,"   	 Pxdmf3d.h5:/%s\n",names[2]);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status); 	 
	  sprintf(str,"         </DataItem>\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	
	  sprintf(str,"      </Geometry>   \n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	
	 // sprintf(str,"      </ Topology>\n");
	  //MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	
	   #if(1)
      sprintf(str,"         <Attribute Name=\"Q\" AttributeType=\"Scalar\" Center=\"Cell\"> \n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);
	  sprintf(str,"   	 <DataItem Name=\"Q\" Dimensions=\"%d \" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",ncube_total);
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
	  sprintf(str,"   	 Pxdmf3d.h5:/Q\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	
	  sprintf(str,"         </DataItem>\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
	  sprintf(str," </Attribute>\n");
	  MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	
	  #endif
	  sprintf(str,"  </Grid>\n");  
      MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
	  sprintf(str,"      </Domain>   \n");
      MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	 
      sprintf(str,"  </Xdmf>\n");
      MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);	  	  
      
    }
  
  

      MPI_File_close(&fp);

}


//
