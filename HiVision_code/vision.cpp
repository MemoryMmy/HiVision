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
#include <hiredis/hiredis.h>
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

bool SortPolygon(polygon_box a, polygon_box b){
	return (bg::get<1,0>(a.first)-bg::get<0,0>(a.first))<(bg::get<1,0>(b.first)-bg::get<0,0>(b.first));
}

class Redis
{
	public:

    	Redis(){}
	
	~Redis()
	{
		this->_connect = NULL;
		this->_reply = NULL;	    	    
	}

	bool connect(string host, int port)
	{
		this->_connect = redisConnect(host.c_str(), port);
		if (this->_connect != NULL && this->_connect->err)
		{
			printf("connect error: %s\n", this->_connect->errstr);
			return false;
		}
		return true;
	}

    string get(string key)
	{
		this->_reply = (redisReply*)redisCommand(this->_connect, "GET %s", key.c_str());
		
		if(this->_reply->type == REDIS_REPLY_NIL)
		{
			printf("Get Data %s failed\n",key.c_str());
			return "";
		}
		string str = this->_reply->str;
		freeReplyObject(this->_reply);
		return str;
	}

	void set(string key, string value)
	{
		redisCommand(this->_connect, "SET %s %s", key.c_str(), value.c_str());
	}

	void zset(string key, char *value, long unsigned int size)
	{ 
		const char *v[4];
		size_t vlen[4];
		v[0] = (char *)"zadd";
		vlen[0] = strlen("zadd");
		v[1] = key.c_str();
		vlen[1] = strlen(key.c_str());
		std::stringstream ss;
		ss << time(0);
		v[2] = ss.str().c_str();
		vlen[2] = ss.str().size();
		v[3] = (const char *)value;
		vlen[3] = size;
	    //~ printf("sizeof(v) :%ld sizeof(v[0]) %d\n ",sizeof(v),sizeof(v[0]));
		redisCommandArgv(this->_connect, sizeof(v) / sizeof(v[0]), v, vlen);
	}
	
	void expire(string key, string time)
	{
		printf("expire %s %s\n", key.c_str(), time.c_str());
		redisCommand(this->_connect, "expire %s %s", key.c_str(), time.c_str());
	}
	
	void pub(string channel, string message)
	{
		redisCommand(this->_connect, "publish %s %s", channel.c_str(), message.c_str());
	}
	
	void del(string key)
	{
		redisCommand(this->_connect, "del %s", key.c_str());
	}
	
	int getllen(string key)
	{
		this->_reply = (redisReply*)redisCommand(this->_connect, "llen %s", key.c_str());
		int str = this->_reply->integer;
		freeReplyObject(this->_reply);
		return str;
	}
	
	string getlindex(string key, int index )
	{
		this->_reply = (redisReply*)redisCommand(this->_connect, "lindex %s %d", key.c_str(),index);
		string str = this->_reply->str;
		freeReplyObject(this->_reply);
		return str;
	}
	void lpush(string key,string value)
	{
		redisCommand(this->_connect, "lpush %s %s", key.c_str(), value.c_str());
	}
	
	string rpop(string key)
	{
		this->_reply = (redisReply*)redisCommand(this->_connect, "rpop %s", key.c_str());
		if(this->_reply->type == REDIS_REPLY_NIL)
			return "";
		string str = this->_reply->str;
		freeReplyObject(this->_reply);
		return str;
	}
	
	string brpop(string key)
	{
		this->_reply = (redisReply*)redisCommand(this->_connect, "brpop %s 0", key.c_str());
		string str = this->_reply->element[1]->str;
		freeReplyObject(this->_reply);
		return str;
	}
	
	bool zget(char *key, char *value)
	{
		this->_reply = (redisReply *)redisCommand(this->_connect, "zrange %s 0 -1", key);
		if (this->_reply->elements<1) 
			return false;
		else
			memcpy(value, this->_reply->element[0]->str, this->_reply->element[0]->len);
		return true;
	}
	
	
	
	private:

    	redisContext* _connect;
		redisReply* _reply;
				
};

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
	char* tmp=new char[256];
	char* tile_params[MAX_TILE_PARAMS];
	int count=0;
	
	rtree_segment * rtree_segment_ptr;
	rtree_point * rtree_point_ptr;
	rtree_psegment * rtree_psegment_ptr;
	rtree_box * rtree_box_ptr;
	
	int x,y,z,fcount;
	double r,rp,rp2,rp4,rbox;
	

	if(myId==0)
	{
		printf("indexPath:%s\n",indexPath);
    	printf("Service Start. cores:%d\n",numProcs);
    }
    while(1){
		
		sprintf(task,"%s",redis->brpop("HiVisiontasklist").c_str());
		try{
			if (strlen(task)>0)
			{
				t1=MPI_Wtime();
				GetList(task, tile_params, (char*)"/", count);
			    //~ r=atof(tile_params[1]);
				z=atoi(tile_params[1]);
				x=atoi(tile_params[2]);
				y=atoi(tile_params[3]); 
                
			    if (tile_params[0][0]=='l')
			    {
					rp=L/(128<<z);
					rp2=0.5303*rp;
					rp4=0.25*rp;
					r=2*rp;
					rbox=r+rp4;
					sprintf(tmp,"%s%s",indexPath,tile_params[0]);
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
							box pointBuffer_small(point(web_mercartor_x-rp2,web_mercartor_y-rp2),point(web_mercartor_x+rp2,web_mercartor_y+rp2));
							rtree_segment::const_query_iterator it=rtree_segment_ptr->qbegin(bgi::intersects(pointBuffer_small));
							if(it!=rtree_segment_ptr->qend())
							{
								buffer_area[buffer_index] += 4;
							}
							else
							{
								box pointBuffer(point(web_mercartor_x-rbox,web_mercartor_y-rbox),point(web_mercartor_x+rbox,web_mercartor_y+rbox));
								std::vector<segment> segment_result;			
								rtree_segment_ptr->query(bgi::intersects(pointBuffer)&&bgi::nearest(point(web_mercartor_x,web_mercartor_y),1),
								 std::back_inserter(segment_result));
	                            fcount=segment_result.size();
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
				else if(tile_params[0][0]=='p')
				{	
					
					rp=L/(128<<z);
					rp2=1.2374*rp;
					rp4=0.25*rp;
					r=2*rp;
					rbox=r+rp4;
					sprintf(tmp,"%s%s",indexPath,tile_params[0]);
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
								fcount=point_result.size();
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
				else if(tile_params[0][0]=='a')
				{
					rp=L/(128<<z);
					rp2=0.5303*rp;
					rp4=0.25*rp;
					r=rp;
					rbox=r+rp4;
					sprintf(tmp,"%s%s",indexPath,tile_params[0]);
				    bi::managed_mapped_file file(bi::open_only, tmp);
					rtree_psegment_ptr=file.find<rtree_psegment>("rtree").first;
					sprintf(tmp,"%s%s_mbr",indexPath,tile_params[0]);
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
								//~ buffer_area[buffer_index] = 1;
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
	                            fcount=segment_result.size();
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
				redis->zset(task,buffer_area,TILE_SIZE*TILE_SIZE);
				while (! redis->zget(task,buffer_tmp))
					redis->zset(task,buffer_area,TILE_SIZE*TILE_SIZE);
				redis->expire(task,"1000");
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
