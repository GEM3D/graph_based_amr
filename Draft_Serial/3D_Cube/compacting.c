#include "typedef.h"
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
	while(cube_temp[0].nbr[i].size()!=0)
	{
	cube_temp[0].nbr[i].pop_back();
    }
}
	
}

