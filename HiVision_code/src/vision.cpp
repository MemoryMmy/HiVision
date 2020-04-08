/*
* 
*    data: 2019-5-20 
*    auther : Ma Mengyu@National University of Defense Technology
*    e-mail : mamengyu10@nudt.edu.cn
*    description : 
*    Hybird-Parallel Visualization Engine
*    run :
* 	 
*/

#include "ogrsf_frmts.h"
#include "ogr_p.h"
#include "cpl_conv.h"
#include "cpl_string.h"
#include "Redis.h"
#include <stdio.h>
#include <cstring>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <dirent.h>
#include <omp.h>
#include <sys/time.h>
#include <regex.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <queue>
#include <utility>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/foreach.hpp>
#include <boost/interprocess/managed_mapped_file.hpp>
#include <boost/lambda/lambda.hpp>
#include <png.h> 
#include <stdlib.h>

#define MAX_NODE_SIZE 8
#define MAX_TILE_PARAMS 256
#define TILE_SIZE 256
#define L 20037508.34

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
namespace bgm = boost::geometry::model;
namespace bi = boost::interprocess;

using namespace std;


typedef bgm::d2::point_xy<double> point;
typedef bgm::box<point> box;
typedef bgm::segment<point> segment;

typedef boost::tuple<segment,unsigned long,bool> polygon_segment;
typedef std::pair<box,unsigned long> polygon_box;  


typedef bgi::quadratic<MAX_NODE_SIZE> params;
typedef bgi::indexable<point> indexable_point;
typedef bgi::equal_to<point> equal_to_point;
typedef bi::allocator<point, bi::managed_mapped_file::segment_manager> allocator_point;
typedef bgi::rtree<point, params, indexable_point, equal_to_point, allocator_point> rtree_point;

typedef bgi::indexable<segment> indexable_segment;
typedef bgi::equal_to<segment> equal_to_segment;
typedef bi::allocator<segment, bi::managed_mapped_file::segment_manager> allocator_segment;
typedef bgi::rtree<segment, params, indexable_segment, equal_to_segment, allocator_segment> rtree_segment;


typedef bgi::indexable<polygon_segment> indexable_psegment;
typedef bgi::equal_to<polygon_segment> equal_to_psegment;
typedef bi::allocator<polygon_segment, bi::managed_mapped_file::segment_manager> allocator_psegment;
typedef bgi::rtree<polygon_segment, params, indexable_psegment, equal_to_psegment, allocator_psegment> rtree_psegment;

typedef bgi::indexable<polygon_box> indexable_box;
typedef bgi::equal_to<polygon_box> equal_to_box;
typedef bi::allocator<polygon_box, bi::managed_mapped_file::segment_manager> allocator_box;
typedef bgi::rtree<polygon_box, params, indexable_box, equal_to_box, allocator_box> rtree_box;

//Polygon with smaller x span has higher priority
bool SortPolygon(polygon_box a, polygon_box b){
	return (bg::get<1,0>(a.first)-bg::get<0,0>(a.first))<(bg::get<1,0>(b.first)-bg::get<0,0>(b.first));
}

void GetList(char* argv, char* result[], char* flag, int& count)
{
	char* string = strdup(argv); 
	char* p;
	int i = 0;
	while((p = strsep(&string, flag)) != NULL)
	{
		result[i] = p;
		i++;
	}
	result[i] = string;
	count = i;
}

//Visualize Linestring Objects
void  LineVision(int z,int x,int y,char* indexPath, char *dataId,char *buffer_area)
{
	char* tmp=new char[256];
	double r,rp,rp2,rp4,rbox;
	rtree_segment * rtree_segment_ptr;
	rp=L/(128<<z);
	rp2=0.5303*rp;
	rp4=0.25*rp;
	r=2*rp;
	rbox=r+rp4;
	sprintf(tmp,"%s%s",indexPath, dataId);
	bi::managed_mapped_file file(bi::open_only, tmp);
	rtree_segment_ptr = file.find<rtree_segment>("rtree").first;
	#pragma omp parallel for num_threads(2) 
	for(int i = 0; i < 256; i++)
	{
		for(int j = 0; j < 256; j ++)
		{
			int buffer_index=i*256+j;
			buffer_area[buffer_index]= 0;
			double web_mercartor_x = (256*x+j+0.5)*rp - L;
			double web_mercartor_y = L-(256*y+i-0.5)*rp;
			box pointbuffer_small(point(web_mercartor_x-rp2,web_mercartor_y-rp2),point(web_mercartor_x+rp2,web_mercartor_y+rp2));
			rtree_segment::const_query_iterator it=rtree_segment_ptr->qbegin(bgi::intersects(pointbuffer_small));
			if(it!=rtree_segment_ptr->qend())
			{
				buffer_area[buffer_index] += 4;
			}
			else
			{
				box pointbuffer(point(web_mercartor_x-rbox,web_mercartor_y-rbox),point(web_mercartor_x+rbox,web_mercartor_y+rbox));
				std::vector<segment> segment_result;			
				rtree_segment_ptr->query(bgi::intersects(pointbuffer)&&bgi::nearest(point(web_mercartor_x,web_mercartor_y),1),
				std::back_inserter(segment_result));
				int fcount=segment_result.size();
				if(fcount>0)
				{
					if(bg::distance(point(web_mercartor_x+rp4,web_mercartor_y+rp4),segment_result.front())<r)
						buffer_area[buffer_index] += 1;
					if(bg::distance(point(web_mercartor_x-rp4,web_mercartor_y-rp4),segment_result.front())<r)
						buffer_area[buffer_index] += 1;
					if(bg::distance(point(web_mercartor_x+rp4,web_mercartor_y-rp4),segment_result.front())<r)
						buffer_area[buffer_index] += 1;
					if(bg::distance(point(web_mercartor_x-rp4,web_mercartor_y+rp4),segment_result.front())<r)
						buffer_area[buffer_index] += 1;
				}
			}
		}
	}				
}

//Visualize Point Objects
void  PointVision(int z,int x,int y,char* indexPath, char *dataId,char *buffer_area)
{	
	char* tmp=new char[256];
	double r,rp,rp2,rp4,rbox;
	rtree_point * rtree_point_ptr;	
	rp=L/(128<<z);
	rp2=1.2374*rp;
	rp4=0.25*rp;
	r=2*rp;
	rbox=r+rp4;
	sprintf(tmp,"%s%s",indexPath,dataId);
	bi::managed_mapped_file file(bi::open_only, tmp);
	rtree_point_ptr = file.find<rtree_point>("rtree").first;
	#pragma omp parallel for num_threads(2) 
	for(int i = 0; i < 256; i++)
	{
		for(int j = 0; j < 256; j ++)
		{
			int buffer_index=i*256+j;
			buffer_area[buffer_index]= 0;
			double web_mercartor_x = (256*x+j+0.5)*rp - L;
			double web_mercartor_y = L-(256*y+i-0.5)*rp;
			box pointBuffer_small(point(web_mercartor_x-rp2,web_mercartor_y-rp2),point(web_mercartor_x+rp2,web_mercartor_y+rp2));	
			rtree_point::const_query_iterator it=rtree_point_ptr->qbegin(bgi::intersects(pointBuffer_small));
						
			if(it!=rtree_point_ptr->qend())
			{
				buffer_area[buffer_index] += 4;
			}
			else
			{
				box pointBuffer(point(web_mercartor_x-rbox,web_mercartor_y-rbox),point(web_mercartor_x+rbox,web_mercartor_y+rbox));		
				std::vector<point> point_result;
				rtree_point_ptr->query(bgi::intersects(pointBuffer)&&bgi::nearest(point(web_mercartor_x,web_mercartor_y),1),
				std::back_inserter(point_result));
				int fcount=point_result.size();
				if(fcount>0)
				{
					if(bg::distance(point(web_mercartor_x+rp4,web_mercartor_y+rp4),point_result.front())<r)
					{
						buffer_area[buffer_index] += 1;
					}
					if(bg::distance(point(web_mercartor_x-rp4,web_mercartor_y-rp4),point_result.front())<r)
					{
						buffer_area[buffer_index] += 1;
					}
					if(bg::distance(point(web_mercartor_x+rp4,web_mercartor_y-rp4),point_result.front())<r)
					{
						buffer_area[buffer_index] += 1;
					}
					if(bg::distance(point(web_mercartor_x-rp4,web_mercartor_y+rp4),point_result.front())<r)
					{
						buffer_area[buffer_index] += 1;
					}
				}
								
			}
										
		}
	}	
}

//Visualize Pologon Objects
void  PologonVision(int z,int x,int y,char* indexPath, char *dataId,char *buffer_area)
{
	char* tmp=new char[256];
	double r,rp,rp2,rp4,rbox;
	rtree_psegment * rtree_psegment_ptr;
	rtree_box * rtree_box_ptr;
	rp=L/(128<<z);
	rp2=0.5303*rp;
	rp4=0.25*rp;
	r=rp;
	rbox=r+rp4;
	sprintf(tmp,"%s%s",indexPath,dataId);
	bi::managed_mapped_file file(bi::open_only, tmp);
	rtree_psegment_ptr=file.find<rtree_psegment>("rtree").first;
	sprintf(tmp,"%s%s_mbr",indexPath,dataId);
	bi::managed_mapped_file file_mbr(bi::open_only, tmp);
	rtree_box_ptr=file_mbr.find<rtree_box>("rtree").first;
	#pragma omp parallel for num_threads(2) 
	for(int i = 0; i < 256; i++)
	{
		for(int j = 0; j < 256; j ++)
		{
			int buffer_index=i*256+j;
			buffer_area[buffer_index]= 0;
			double web_mercartor_x = (256*x+j+0.5)*rp - L;
			double web_mercartor_y = L-(256*y+i-0.5)*rp;
			std::vector<polygon_box> box_result;
			rtree_box_ptr->query(bgi::intersects(point(web_mercartor_x,web_mercartor_y)),std::back_inserter(box_result));
			sort(box_result.begin(),box_result.end(),SortPolygon);				
			BOOST_FOREACH(polygon_box const &v, box_result)
			{
				unsigned long pID= v.second;
				box pBox= v.first;
				double minx= bg::get<0,0>(pBox);
				double maxx= bg::get<1,0>(pBox);
				int lcount=0;		
				std::vector<polygon_segment> psegment_result;
				if(web_mercartor_x-minx<maxx-web_mercartor_x)									
					rtree_psegment_ptr->query(bgi::intersects(segment(point(minx,web_mercartor_y),point(web_mercartor_x,web_mercartor_y))),
					std::back_inserter(psegment_result));
				else
					rtree_psegment_ptr->query(bgi::intersects(segment(point(web_mercartor_x,web_mercartor_y),point(maxx,web_mercartor_y))),
				std::back_inserter(psegment_result));
				BOOST_FOREACH(polygon_segment const &sv, psegment_result)
				{
					bool level=boost::get<2>(sv);
					unsigned long pIDtmp= boost::get<1>(sv);
					if (level&&pID==pIDtmp)
						lcount++;
				}
				if(lcount%2>0)
				{	
					buffer_area[buffer_index]=4;
					break;
				}
			}										
			box pointBuffer_small(point(web_mercartor_x-rp2,web_mercartor_y-rp2),point(web_mercartor_x+rp2,web_mercartor_y+rp2));				
			rtree_psegment::const_query_iterator it=rtree_psegment_ptr->qbegin(bgi::intersects(pointBuffer_small));
						
			if(it!=rtree_psegment_ptr->qend())
			{
				if(buffer_area[buffer_index]>0)
					buffer_area[buffer_index] += 4;
				else
					buffer_area[buffer_index] -= 4;
			}
			else
			{							
				box pointBuffer(point(web_mercartor_x-rbox,web_mercartor_y-rbox),point(web_mercartor_x+rbox,web_mercartor_y+rbox));
                std::vector<polygon_segment> segment_result;			
				rtree_psegment_ptr->query(bgi::intersects(pointBuffer)&&bgi::nearest(point(web_mercartor_x,web_mercartor_y),1),
				std::back_inserter(segment_result));
                int fcount=segment_result.size();
				if(fcount>0)
				{
					segment pSegment= boost::get<0>(segment_result.front());
					if(buffer_area[buffer_index]>0)
					{
						if(bg::distance(point(web_mercartor_x+rp4,web_mercartor_y+rp4),pSegment)<r)
						{
							buffer_area[buffer_index] += 1;
						}
						if(bg::distance(point(web_mercartor_x-rp4,web_mercartor_y-rp4),pSegment)<r)
						{
							buffer_area[buffer_index] += 1;
						}
						if(bg::distance(point(web_mercartor_x+rp4,web_mercartor_y-rp4),pSegment)<r)
						{
							buffer_area[buffer_index] += 1;
						}
						if(bg::distance(point(web_mercartor_x-rp4,web_mercartor_y+rp4),pSegment)<r)
						{
							buffer_area[buffer_index] += 1;
						}	
					}else
					{
						if(bg::distance(point(web_mercartor_x+rp4,web_mercartor_y+rp4),pSegment)<r)
						{
							buffer_area[buffer_index] -= 1;
						}
						if(bg::distance(point(web_mercartor_x-rp4,web_mercartor_y-rp4),pSegment)<r)
						{
							buffer_area[buffer_index] -= 1;
						}
						if(bg::distance(point(web_mercartor_x+rp4,web_mercartor_y-rp4),pSegment)<r)
						{
							buffer_area[buffer_index] -= 1;
						}
						if(bg::distance(point(web_mercartor_x-rp4,web_mercartor_y+rp4),pSegment)<r)
						{
							buffer_area[buffer_index] -= 1;
						}
								
					}
							
				}
			}			
												
		}
	}					
}


int main( int nArgc, char ** papszArgv )
{ 
	int myId, numProcs;
	MPI_Init(&nArgc,&papszArgv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myId);
	MPI_Comm_size(MPI_COMM_WORLD,&numProcs);
	double t1,t2;
	
	char* indexPath = papszArgv[1];
	char* redisHost=papszArgv[2];
	int redisPort=atoi(papszArgv[3]);
	
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
	CPLSetConfigOption("SHAPE_ENCODING", "UTF-8");
	
	Redis *redis = new Redis();
	if(!redis->connect(redisHost, redisPort))
	{
		printf("connect redis error!\n");
		exit(0);
		MPI_Finalize();
	}
	
	char *buffer_area = (char*)malloc(TILE_SIZE*TILE_SIZE);
	char *buffer_tmp = (char*)malloc(TILE_SIZE*TILE_SIZE);
	char* task= new char[256];
	char* tile_params[MAX_TILE_PARAMS];
	int count=0;
	int x,y,z;
	if(myId==0)
	{
		printf("indexPath:%s\n",indexPath);
    	printf("Service Start. cores:%d\n",numProcs);
    }
    while(1){
		sprintf(task,"%s",redis->brpop("HiVisiontasklist").c_str());//Get tasks from Task Pool
		try{
			if (strlen(task)>0)
			{
				t1=MPI_Wtime();
				GetList(task, tile_params, (char*)"/", count);
				z=atoi(tile_params[1]);
				x=atoi(tile_params[2]);
				y=atoi(tile_params[3]); 
			    if (tile_params[0][0]=='l')
					LineVision(z,x,y,indexPath,tile_params[0],buffer_area);					
				else if(tile_params[0][0]=='p')
					PointVision(z,x,y,indexPath,tile_params[0],buffer_area);
				else if(tile_params[0][0]=='a')
					PologonVision(z,x,y,indexPath,tile_params[0],buffer_area);
				redis->zset(task,buffer_area,TILE_SIZE*TILE_SIZE);
				while (! redis->zget(task,buffer_tmp))
					redis->zset(task,buffer_area,TILE_SIZE*TILE_SIZE);
				redis->expire(task,"1000");//Expire time: 1000s
			    redis->pub("HiVisiontiles",task);
				t2 = MPI_Wtime();
				printf("tile-%s-%d-%f\n",task, myId, t2-t1);
			}
		}
		catch(...)
		{
			printf("Error task %s \n",task);		
		}
	}
	MPI_Finalize();
}
