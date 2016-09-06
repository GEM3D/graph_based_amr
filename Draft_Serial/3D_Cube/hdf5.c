#include "typedefs.h"
#include "hdf5.h"

void write_xdmf_serial(Cube& cube, int L, int M, int N);

void write_hdf5_serial(Cube& cube, int L, int M, int N)
{
	
	unsigned int i,j,k,l;
	
	
	hid_t     file_id;
    file_id = H5Fcreate("xdmf3d.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	
	double *xtemp =NULL; 
    double *ytemp =NULL;
    double *ztemp =NULL;
    xtemp=(double *) calloc(L*M*N, sizeof(double));
    ytemp=(double *) calloc(L*M*N, sizeof(double));
    ztemp=(double *) calloc(L*M*N, sizeof(double));
    
    // for now put Q here
    
    double *q =NULL;
    q=(double *) calloc(L*M*N,sizeof(double));

    
    


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
    
    hid_t     dataset_id, dataspace_id;
    hsize_t   dims[3];
    herr_t    status;

	char str0[50];
	char str1[50];
	char str2[50];
	char str3[50];
		
	
	for(i=0;i<cube.size();i++)
	//for(i=0;i<2;i++)
	{
		//i=0;
     
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
	       xtemp[index]=Xa+Xh*j;	      
	       //printf("%lf\n",xtemp[index]);      
	       ytemp[index]=Ya+Yh*k;
	       ztemp[index]=Za+Zh*l;
	       q[index]=(double)i;    
            index++;
        }
      }
	}	
	 		 
    sprintf(str0, "/X%d", i);
    sprintf(str1, "/Y%d", i);
    sprintf(str2, "/Z%d", i);
    sprintf(str3, "/Q%d", i);
    /*
    printf( "%s\n", str0);
    printf( "%s\n", str1);
    printf( "%s\n", str2);
    printf( "%s\n", str3);
    */
   
    /* Write separate coordinate arrays for the x and y and z coordinates. */
    
        dims[0] = L;
        dims[1] = M;
        dims[2] = N;
        
        dataspace_id = H5Screate_simple(3, dims, NULL);
        dataset_id = H5Dcreate(file_id, str0, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT, xtemp);
        status = H5Dclose(dataset_id);
        status = H5Sclose(dataspace_id);
 
        dataspace_id = H5Screate_simple(3, dims, NULL);
        dataset_id = H5Dcreate(file_id, str1, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT, ytemp);
        status = H5Dclose(dataset_id);
        status = H5Sclose(dataspace_id);
	    
	    dataspace_id = H5Screate_simple(3, dims, NULL);
        dataset_id = H5Dcreate(file_id, str2, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT, ztemp);
        status = H5Dclose(dataset_id);
        status = H5Sclose(dataspace_id);
        
        
       dataspace_id = H5Screate_simple(3, dims, NULL);
       dataset_id = H5Dcreate(file_id, str3, H5T_NATIVE_DOUBLE,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
       status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, q);
       status = H5Dclose(dataset_id);
       status = H5Sclose(dataspace_id);

        
        
       // printf("%d\n",i);

	}
	
	write_xdmf_serial(cube,L,M, N);
	
	//printf("%d\n",cube.size());
		
	free(xtemp)	;
	free(ytemp)	;
	free(ztemp)	;
	free(q);
}

void write_xdmf_serial(Cube& cube, int L, int M, int N)
{
	unsigned int i;
	char str0[50];
	char str1[50];
	char str2[50];
	char str3[50];
	char str4[50];
	char str5[50];
	
    FILE *xmf=NULL;
    xmf = fopen("xdmf3d.xmf", "w");
    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
    fprintf(xmf, " <Domain>\n");
    fprintf(xmf, "   <Grid Name=\"AMR0\" GridType=\"Collection\" CollectionType=\"Spatial\">\n");
    
    for(i=0;i<cube.size();i++)
   //   for(i=0;i<2;i++)
    {
	  sprintf(str0, "/X%d", i);
	  sprintf(str1, "/Y%d", i);
      sprintf(str2, "/Z%d", i);
      sprintf(str3, "mesh%d", i);
      sprintf(str4, "/Q%d", i);
   /* printf( "%s\n", str0);
    printf( "%s\n", str1);
    printf( "%s\n", str2);
    printf( "%s\n", str3);  
    printf( "%s\n", str4);  */
	fprintf(xmf, "   <Grid Name=\"%s\" GridType=\"Uniform\">\n",str3);
	fprintf(xmf, "     <Topology TopologyType=\"3DSMesh\" NumberOfElements=\"%d %d %d\"/>\n", L, M, N);
    fprintf(xmf, "     <Geometry GeometryType=\"X_Y_Z\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", L, M, N);
    fprintf(xmf, "        xdmf3d.h5:%s\n",str0);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", L, M, N);
    fprintf(xmf, "        xdmf3d.h5:%s\n",str1);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", L, M, N);
    fprintf(xmf, "        xdmf3d.h5:%s\n",str2);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Geometry>\n");
    // add your variables here
    fprintf(xmf, "     <Attribute Name=\"Q\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", L, M, N);
    fprintf(xmf, "        xdmf3d.h5:%s\n",str4);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");    
    fprintf(xmf, "   </Grid>\n");		
	}
    fprintf(xmf, "   </Grid>\n");
    fprintf(xmf, " </Domain>\n");
    fprintf(xmf, "</Xdmf>\n");
    fclose(xmf);
    
    

}
