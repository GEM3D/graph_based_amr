
//======================================================================
//
//                Compact Version of Adding Kids
//
//======================================================================


void add_kids_compact(Cube& cube, Cube& cube_temp,int id)
{
 unsigned int i,j;
	
 
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

double xmin=cube_temp[0].xyz[0];  	
double ymin=cube_temp[0].xyz[2];
double zmin=cube_temp[0].xyz[4];

double xmax=cube_temp[0].xyz[1];  	
double ymax=cube_temp[0].xyz[3];
double zmax=cube_temp[0].xyz[5];

printf("%lf %lf %lf\n",xmid,ymid,zmid);
// inside the refined element update 
//======================================================================

//                     kid 0 is located at 0,0,0 

//======================================================================

unsigned int map_index[8];
map_index[0]=id;

for(i=1;i<8;i++)
{
	map_index[i]=nquad+1;
}

double Xmin[8]={xmin,xmid,xmid,xmin,xmin,xmid,xmid,xmin};
double Xmax[8]={xmid,xmax,xmax,xmid,xmid,xmax,xmax,xmid};

double Ymin[8]={ymin,ymin,ymid,ymid,ymin,ymin,ymid,ymid};
double Ymax[8]={ymid,ymid,ymax,ymax,ymid,ymid,ymax,ymax};

double Zmin[8]={zmin,zmin,zmin,zmin,zmid,zmid,zmid,zmid};
double Zmax[8]={zmid,zmid,zmid,zmid,zmax,zmax,zmax,zmax};

for(i=0;i<8;i++)
{
  cube[map_index[i]].xyz[0]=Xmin[i];  	
  cube[map_index[i]].xyz[1]=Xmax[i];
  
  cube[map_index[i]].xyz[2]=Ymin[i];
  cube[map_index[i]].xyz[3]=Ymax[i];
  
  cube[map_index[i]].xyz[4]=Zmin[i];
  cube[map_index[i]].xyz[5]=Zmax[i];
 
}
  
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
	     xc_id= cube[id].xyz[0]+cube[id].xyz[1];		 
		 yc_id= cube[id].xyz[2]+cube[id].xyz[3];  
		 
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
	 printf("idx =%d\n",idx);
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
	 
	  cube[index].nbr[2].push_back(cube_temp[0].nbr[2].at(idx)); 
	  
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
	 
	  cube[index].nbr[2].push_back(cube_temp[0].nbr[2].at(idx)); 
	  
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
	 
	  cube[index].nbr[1].push_back(cube_temp[0].nbr[1].at(idx)); 
	  
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
	for(j=0;j<cube_temp[0].nbr[i].size();j++)
	{
	cube_temp[0].nbr[i].pop_back();
    }
}
	
}

//=============================================================================
//
//==============================================================================

/*****************************************************************************
*
* Copyright (c) 2000 - 2006, The Regents of the University of California
* Produced at the Lawrence Livermore National Laboratory
* All rights reserved.
*
* This file is part of VisIt. For details, see http://www.llnl.gov/visit/. The
* full copyright notice is contained in the file COPYRIGHT located at the root
* of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
*
* Redistribution  and  use  in  source  and  binary  forms,  with  or  without
* modification, are permitted provided that the following conditions are met:
*
*  - Redistributions of  source code must  retain the above  copyright notice,
*    this list of conditions and the disclaimer below.
*  - Redistributions in binary form must reproduce the above copyright notice,
*    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
*    documentation and/or materials provided with the distribution.
*  - Neither the name of the UC/LLNL nor  the names of its contributors may be
*    used to  endorse or  promote products derived from  this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED.  IN  NO  EVENT  SHALL  THE  REGENTS  OF  THE  UNIVERSITY OF
* CALIFORNIA, THE U.S.  DEPARTMENT  OF  ENERGY OR CONTRIBUTORS BE  LIABLE  FOR
* ANY  DIRECT,  INDIRECT,  INCIDENTAL,  SPECIAL,  EXEMPLARY,  OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

/*
   WARNING, YOU ALSO NEED THE EXPAT XML PARSER LIBRARY TO COMPILE THIS CODE 
   You can find expat here... http://expat.sourceforge.net/
*/


#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
extern int errno;

#include <assert.h>
#if !defined(NDEBUG) && !defined(DEBUG)
#undef assert
#define assert(a)
#endif

#include <expat.h>

#include <silo.h>


typedef enum _xmlParseState_t
{
    None = 0,
    AMRDecomp,
    Level,
    Patch
} xmlParseState_t;

typedef struct patchInfo_t
{
    int   idx;
    int   vid;
    double rank;
    int level;
    int logSize[3];
    int logExtents[6];
    double spatExtents[6];
    int numChildren;
    int *children;
} patchInfo_t;

typedef struct amrConfig_t
{
    char *meshName;
    int numDims;
    int numLevels;
    int numPatches;
    double dx, dy, dz; /* all -1 ==> curvilinear */

    xmlParseState_t currentState;
    int currentLevel;
    int currentPatchOnLevel;
    int currentChildForPatch;

    int *numPatchesOnLevel;
    int **levels;
    int *ratios;

    patchInfo_t *patches;
    int currentPatch;
} amrConfig_t;

#ifdef XML_LARGE_SIZE
#if defined(XML_USE_MSC_EXTENSIONS) && _MSC_VER < 1400
#define XML_FMT_INT_MOD "I64"
#else
#define XML_FMT_INT_MOD "ll"
#endif
#else
#define XML_FMT_INT_MOD "l"
#endif

#define BUFFSIZE        8192
char Buff[BUFFSIZE];

#define INT     ((int)'i')
#define DBL     ((int)'d')
#define STR     ((int)'s'), 1
#define EOA     "EOA"

static int
getAttrVals(const char **attrs, ...)
{
    int i, j, nvals;
    const char *attrName, *p;
    char *pend;
    int type;
    va_list ap;

    va_start(ap, attrs);
    while (1)
    {
        attrName = va_arg(ap, const char *);
        if (strncmp(attrName, EOA, 3) == 0)
            break;

        /* find this attr's name in the list of xml attrs */
        for (i = 0; attrs[i]; i += 2)
        {
            if (strcmp(attrs[i], attrName) == 0)
                break;
        }
        if (!attrs[i])
            break;

        type = va_arg(ap, int);
        nvals = va_arg(ap, int);
        p = attrs[i+1];
        switch (type)
        {
            case 'i':
            {
                int *ptr = va_arg(ap, int *);
                for (j = 0; j < nvals; j++)
                {
                    ptr[j] = strtol(p,&pend,10);
                    p = pend;
                }
                break;
            }
            case 'd':
            {
                double *ptr = va_arg(ap, double *);
                for (j = 0; j < nvals; j++)
                {
                    ptr[j] = strtod(p,&pend);
                    p = pend;
                }
                break;
            }
            case 's':
            {
                char **ptr = va_arg(ap, char **);
                *ptr = strdup(attrs[i+1]); 
                break;
            }
            default: return -1;
        }
    }

    va_end(ap);
    return 0;
}

static void XMLCALL
start(void *data, const char *el, const char **attr)
{
    int i;
    amrConfig_t *ud = (amrConfig_t *) data;
    
    if      (strcmp(el, "AMRDecomp") == 0)
    {
        getAttrVals(attr, "meshName", STR, &ud->meshName,
                          "numDims", INT, 1, &ud->numDims,
                          "numLevels", INT, 1, &ud->numLevels,
                          "numPatches", INT, 1, &ud->numPatches,
                    EOA);
        ud->currentLevel = -1;
        ud->levels = (int **) malloc(ud->numLevels * sizeof(int*));
        ud->ratios = (int *) malloc(ud->numLevels * ud->numDims * sizeof(int));
        ud->numPatchesOnLevel = (int *) malloc(ud->numLevels * sizeof(int));
        ud->patches = (patchInfo_t *) malloc(ud->numPatches * sizeof(patchInfo_t));
        ud->currentState = AMRDecomp;
    }
    else if (strcmp(el, "Level") == 0)
    {
        int cl;
        ud->currentLevel++;
        cl = ud->currentLevel;
        getAttrVals(attr, "numPatches", INT, 1, &(ud->numPatchesOnLevel[cl]),
                          "ratios", INT, ud->numDims, &(ud->ratios[ud->numDims*cl]),
                    EOA);
        ud->levels[cl] = (int *) malloc(ud->numPatchesOnLevel[cl] * sizeof(int));
        ud->currentPatchOnLevel = 0;
        ud->currentPatch = -1;
        ud->currentState = Level;
    }
    else if (strcmp(el, "Patch") == 0)
    {
        ud->currentPatch++;
        patchInfo_t *p = &(ud->patches[ud->currentPatch]);
        getAttrVals(attr, "iDx", INT, 1, &(p->idx),
                          "vId", INT, 1, &(p->vid),
                          "level", INT, 1, &(p->level),
                          "rank", DBL, 1, &(p->rank),
                          "numChildren", INT, 1, &(p->numChildren),
                          "logSize", INT, ud->numDims, &(p->logSize),
                          "logExtents", INT, 2*(ud->numDims), &(p->logExtents),
                          "spatExtents", DBL, 2*(ud->numDims), &(p->spatExtents),
                    EOA);
        if (p->numChildren == 0)
            p->children = 0;
        else
            p->children = (int *) malloc(p->numChildren * sizeof(int));
        ud->currentChildForPatch = 0;
        ud->currentState = Patch;
    } 
}

static void XMLCALL
end(void *data, const char *el)
{
}

static void XMLCALL
cdstart(void *userData)
{
}

static void XMLCALL
cdend(void *userData)
{
}

static int
readVals(const char *s, int type, void *pvals)
{
    const char *p1;
    char *p2;
    int n, done;

    p1 = s;
    n = 0;
    done = 0;
    errno = 0;
    while (!done)
    {
        switch (type)
        {
            case 'i':
            {
                int *iptr = (int *) pvals;
                int ival = strtol(p1, &p2, 10);
                if (ival == 0 && p1 == p2 || errno != 0)
                {
                    done = 1;
                    break;
                }
                if (*p2 == '\0')
                    done = 1;
                if (iptr) iptr[n] = ival;
                break;
            }
            case 'd':
            {
                double *dptr = (double *) pvals;
                double dval = strtod(p1, &p2);
                if (dval == 0 && p1 == p2 || errno != 0)
                {
                    done = 1;
                    break;
                }
                if (*p2 == '\0')
                    done = 1;
                if (dptr) dptr[n] = dval;
                break;
            }
        }
        p1 = p2;
        if (!done) n++;
    }
    return n;
}

static void XMLCALL
chardata(void *userData, const XML_Char *s, int len)
{
    amrConfig_t *ud = (amrConfig_t *) userData;
    int i, val, done;
    char *s1, *p1, *p2;

    if (len == 0)
        return;

    s1 = malloc(len+1);
    strncpy(s1, s, len);
    s1[len] = '\0';

    /* if its all whitespace, ignore it */
    if (strspn(s1, " \\\f\\\n\\\r\\\t\\\v") == len)
    {
        free(s1);
        return;
    }

    switch (ud->currentState)
    {
        case Level:
        {
            int cl = ud->currentLevel;
            int cp = ud->currentPatchOnLevel;
            int *levels = &(ud->levels[cl][cp]);
            int npatches = readVals(s1, INT, levels);
            int i;
            for (i = 0; i < npatches; i++)
                printf("ud->levels[%d][%d]=%d\n", cl, cp+i, ud->levels[cl][cp+i]);
            ud->currentPatchOnLevel += npatches;
            break;
        }
        case Patch:
        {
            int cp = ud->currentPatch;
            int cc = ud->currentChildForPatch;
            int *children = &(ud->patches[cp].children[cc]);
            int nchildren = readVals(s1, INT, children);
            int i;
            for (i = 0; i < nchildren; i++)
                printf("ud->patches[%d].children[%d]=%d\n", cp, cc+i, ud->patches[cp].children[cc+i]);
            ud->currentChildForPatch += nchildren;
            break;
        }
    }

    free(s1);
}

int ProcessXMLAMRConfigFile(const char *xmlFileName, amrConfig_t *amrconfig)
{
    FILE *acf = fopen(xmlFileName, "r");
    XML_Parser p = XML_ParserCreate(NULL);
    if (!p)
        return -1;

    XML_SetUserData(p, amrconfig);
    XML_SetElementHandler(p, start, end);
    XML_SetCdataSectionHandler(p, cdstart, cdend);
    XML_SetCharacterDataHandler(p, chardata);

    for (;;) {
        int done;
        int len;

      len = (int)fread(Buff, 1, BUFFSIZE, acf);
      if (ferror(acf))
          return -1;
      done = feof(acf);

      if (XML_Parse(p, Buff, len, done) == XML_STATUS_ERROR) {
        fprintf(stderr, "Parse error at line %" XML_FMT_INT_MOD "u:\n%s\n",
              XML_GetCurrentLineNumber(p),
              XML_ErrorString(XML_GetErrorCode(p)));
          return -1;
      }

      if (done)
          break;
    }
    XML_ParserFree(p);
    return 0;
}

/*-------------------------------------------------------------------------
 * Function: main
 *
 * Purpose: Add an mrgtree object to an existing silo file containing a
 * multi-mesh object whose individual pieces form the patches of some AMR
 * mesh.
 *
 * Return:	0
 *
 * Programmer:	Mark C. Miller, Wed Aug 27 09:17:38 PDT 2008
 *
 *-------------------------------------------------------------------------
 */
int
main(int argc, char *argv[])
{
    int i;
    int copy = 1;
    DBfile *dbfile = 0;
    DBmrgtree *mrgTree;
    DBmultimesh *mm;
    DBoptlist *optList;
    char *siloFileName = 0;
    char *amrconfigFileName = "amr_config.xml";
    amrConfig_t amrconf;
    char tmpName[256];
    char lvlMapsName[256];
    char chldMapsName[256];

    assert(copy==0);

    /* Parse command-line */
    for (i=1; i<argc; i++)
    {
        if (strstr(argv[i], "help"))
        {
            print("This application will take an existing Silo file containing a\n
                   multi-block quadmesh and an xml file indicating how the multi-block\n
                   pieces fit together in an AMR hierarchy and add the mrgree to the\n
                   silo file.\n");
            print("Usage: add_amr_mrgtree <silo-filename> <amr-config-xml-filename>\n");
            exit(1);
        }
	else if (!strcmp(argv[i], "-nocp"))
        {
            copy = 0;
        }
        else if (siloFileName == 0)
        {
            siloFileName = strdup(argv[i]);
	}
        else
        {
            amrconfigFileName = strdup(argv[i]);
        }
    }

    DBShowErrors(DB_ABORT, NULL);


    /* by default, we make a copy of the specified file */
    if (copy)
    {
        char syscmd[256];

        /* copy the file in the filesystem */
        snprintf(syscmd, sizeof(syscmd), "cp -f %s %s_wmrgtree", siloFileName, siloFileName);
        system(syscmd);

        /* use syscmd as buffer for constructing new filename */
        snprintf(syscmd, sizeof(syscmd), "%s_wmrgtree", siloFileName);
        free(siloFileName);
        siloFileName = strdup(syscmd);
    }

    /* open and process the amr configuration file */
    ProcessXMLAMRConfigFile(amrconfigFileName, &amrconf);

    /* open the silo file */
    dbfile = DBOpen(siloFileName, DB_UNKNOWN, DB_APPEND);

    /* get the specific multi-mesh object we want to add an mrg tree too */
    mm = DBGetMultimesh(dbfile, amrconf.meshName);
    sprintf(tmpName, "%s_wmrgtree", amrconf.meshName);
    optList = DBMakeOptlist(10);
    DBAddOption(optList, DBOPT_MRGTREE_NAME, "mrgTree");
    DBPutMultimesh(dbfile, tmpName, mm->nblocks, mm->meshnames, mm->meshtypes, optList);
    DBClearOptlist(optList);

#warning HACK FOR SINGLE VARIABLE
    {
        DBmultivar *mv = DBGetMultivar(dbfile, "foo");
        DBoptlist *optList2 = DBMakeOptlist(10);
        DBAddOption(optList2, DBOPT_MMESH_NAME, tmpName);
        DBPutMultivar(dbfile, "foo_wmrgtree", mv->nvars, mv->varnames, mv->vartypes, optList2);
	DBFreeOptlist(optList2);
    }
    
    /* write this multi-mesh object back to the file, with a different name
       and attache the mrg tree name to it */

    /* Write the groupel maps that specify which blocks of the multi-block mesh
       belong to which levels */
    {
        int *segTypes = (int *) malloc(amrconf.numLevels * sizeof(int));
        for (i = 0; i < amrconf.numLevels; i++)
            segTypes[i] = DB_BLOCKCENT;
        sprintf(lvlMapsName, "%s_wmrgtree_lvlMaps", amrconf.meshName);
        DBPutGroupelmap(dbfile, lvlMapsName, amrconf.numLevels, segTypes, amrconf.numPatchesOnLevel,
            0, amrconf.levels, 0, 0, 0); 
        free(segTypes);
    }

    /* Write the groupel maps that specify which blocks are children (refinements)
       of other blocks */
    {
        int *segTypes = (int *) malloc(amrconf.numPatches * sizeof(int));
        int *segLens = (int *) malloc(amrconf.numPatches * sizeof(int));
        int **segData = (int **) malloc(amrconf.numPatches * sizeof(int*));
        for (i = 0; i < amrconf.numPatches; i++)
        {
            segTypes[i] = DB_BLOCKCENT;
            segLens[i] = amrconf.patches[i].numChildren;
            segData[i] = amrconf.patches[i].children;
        }
        sprintf(chldMapsName, "%s_wmrgtree_chldMaps", amrconf.meshName);
        DBPutGroupelmap(dbfile, chldMapsName, amrconf.numPatches, segTypes, segLens,
            0, segData, 0, 0, 0); 
        free(segTypes);
        free(segLens);
        free(segData);
    }

    /* Create an mrg tree for inclusion in the file */
    mrgTree = DBMakeMrgtree(DB_MULTIMESH, 0, 2, 0);

    /* Add a region for AMR configuration */
    DBAddRegion(mrgTree, "amr_decomp", 0, 2, 0, 0, 0, 0, 0, 0); 
    DBSetCwr(mrgTree, "amr_decomp");
    DBAddRegion(mrgTree, "levels", 0, amrconf.numLevels, 0, 0, 0, 0, 0, 0); 
    DBSetCwr(mrgTree, "levels");

    /* Handle the regions representing each level in the AMR mesh */
    {
        char *levelRegnNames[1];
        int *segTypes = (int *) malloc(amrconf.numLevels * sizeof(int));
        int *segIds = (int *) malloc(amrconf.numLevels * sizeof(int));
        for (i = 0; i < amrconf.numLevels; i++)
        {
            segIds[i] = i;
            segTypes[i] = DB_BLOCKCENT;
        }
        levelRegnNames[0] = "@level%d@n";
        DBAddRegionArray(mrgTree, amrconf.numLevels, levelRegnNames, 0, lvlMapsName, 1,
            segIds, amrconf.numPatchesOnLevel, segTypes, 0);
    }
    DBSetCwr(mrgTree, "..");

    DBAddRegion(mrgTree, "patches", 0, amrconf.numPatches, 0, 0, 0, 0, 0, 0); 
    DBSetCwr(mrgTree, "patches");

    /* Handle the regions representing each patch */
    {
        char *patchRegnNames[1];
        int *segTypes = (int *) malloc(amrconf.numPatches * sizeof(int));
        int *segIds = (int *) malloc(amrconf.numPatches * sizeof(int));
        int *segLens = (int *) malloc(amrconf.numPatches * sizeof(int));
        for (i = 0; i < amrconf.numPatches; i++)
        {
            segIds[i] = i;
            segTypes[i] = DB_BLOCKCENT;
            segLens[i] = amrconf.patches[i].numChildren;
        }
        patchRegnNames[0] = "@patch%d@n";
        DBAddRegionArray(mrgTree, amrconf.numPatches, patchRegnNames, 0, chldMapsName, 1,
            segIds, segLens, segTypes, 0);
    }

    {
        char *mrgv_onames[5];
	sprintf(tmpName, "%s_wmrgtree_lvlRatios", amrconf.meshName);
        mrgv_onames[0] = strdup(tmpName);
	sprintf(tmpName, "%s_wmrgtree_ijkExts", amrconf.meshName);
        mrgv_onames[1] = strdup(tmpName);
	sprintf(tmpName, "%s_wmrgtree_xyzExts", amrconf.meshName);
        mrgv_onames[2] = strdup(tmpName);
        mrgv_onames[3] = strdup("rank");
        mrgv_onames[4] = 0;

        DBAddOption(optList, DBOPT_MRGV_ONAMES, mrgv_onames);
        DBPutMrgtree(dbfile, "mrgTree", "amr_mesh", mrgTree, optList);
        DBFreeMrgtree(mrgTree);
    }

    /* Output level refinement ratios as an mrg variable on the array of regions
       representing the levels */
    {
        char *compnames[3] = {"iRatio","jRatio","kRatio"};
        char *levelRegnNames[1];
        int *data[3];
        for (i = 0; i < amrconf.numDims; i++)
            data[i] = (int *) malloc(amrconf.numLevels * sizeof(int));
        for (i = 0; i < amrconf.numLevels; i++)
        {
            data[0][i] = amrconf.ratios[i*amrconf.numDims+0];
            data[1][i] = amrconf.ratios[i*amrconf.numDims+1];
            if (amrconf.numDims == 3)
                data[2][i] = amrconf.ratios[i*amrconf.numDims+2];
        }
        levelRegnNames[0] = "@level%d@n";
	sprintf(tmpName, "%s_wmrgtree_lvlRatios", amrconf.meshName);
        DBPutMrgvar(dbfile, tmpName, "mrgTree", amrconf.numDims,
            compnames, amrconf.numLevels, levelRegnNames, DB_INT, (void**)data, 0);
        for (i = 0; i < amrconf.numDims; i++)
            free(data[i]);
    }

    /* Output logical extents of the patches as an mrg variable on the
       array of regions representing the patches */
    {
        char *compnames[6] = {"iMin","iMax","jMin","jMax","kMin","kMax"};
        char *scompnames[6] = {"xMin","xMax","yMin","yMax","zMin","zMax"};
        char *patchRegnNames[1];
        int *data[6];
        float *rdata[1];
        float *sdata[6];
        patchRegnNames[0] = "@patch%d@n";
        for (i = 0; i < 2*amrconf.numDims; i++)
        {
            data[i] = (int *) malloc(amrconf.numPatches * sizeof(int));
            sdata[i] = (float *) malloc(amrconf.numPatches * sizeof(float));
        }
        rdata[0] = (float *) malloc(amrconf.numPatches * sizeof(float));
        for (i = 0; i < amrconf.numPatches; i++)
        {
            data[0][i] = amrconf.patches[i].logExtents[0];
            data[1][i] = amrconf.patches[i].logExtents[1];
            data[2][i] = amrconf.patches[i].logExtents[2];
            data[3][i] = amrconf.patches[i].logExtents[3];
            sdata[0][i] = amrconf.patches[i].spatExtents[0];
            sdata[1][i] = amrconf.patches[i].spatExtents[1];
            sdata[2][i] = amrconf.patches[i].spatExtents[2];
            sdata[3][i] = amrconf.patches[i].spatExtents[3];

            if (amrconf.numDims == 3)
            {
                data[4][i] = amrconf.patches[i].logExtents[4];
                data[5][i] = amrconf.patches[i].logExtents[5];
                sdata[4][i] = amrconf.patches[i].spatExtents[4];
                sdata[5][i] = amrconf.patches[i].spatExtents[5];
            }
            rdata[0][i] = amrconf.patches[i].rank;
        }
	sprintf(tmpName, "%s_wmrgtree_ijkExts", amrconf.meshName);
        DBPutMrgvar(dbfile, tmpName, "mrgTree", 2*amrconf.numDims,
            compnames, amrconf.numPatches, patchRegnNames, DB_INT, (void**)data, 0);
	sprintf(tmpName, "%s_wmrgtree_xyzExts", amrconf.meshName);
        DBPutMrgvar(dbfile, tmpName, "mrgTree", 2*amrconf.numDims,
            scompnames, amrconf.numPatches, patchRegnNames, DB_FLOAT, (void**)sdata, 0);
        for (i = 0; i < 2*amrconf.numDims; i++)
        {
            free(data[i]);
            free(sdata[i]);
        }
        DBPutMrgvar(dbfile, "rank", "mrgTree", 1, 0,
            amrconf.numPatches, patchRegnNames, DB_FLOAT, (void**)rdata, 0);
        free(rdata[0]);
    }

    DBClose(dbfile);

    return (0);
}







