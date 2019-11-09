/*
* 
*    data: 2018-11-20 
*    auther : Ma Mengyu@DBRG
*    description :
*    create indexes(Polygon objects)
*    dependency :
*    gdal2.10  boost1.64 
*    run :
* 	 ./Polygonindex --shp ./alaska_OSM/gis_osm_buildings_a_free_1.shp --output ../indexes/abuilding --rediskey abuilding --redishost 127.0.0.1 --redisport 6379
*/


#include "ogrsf_frmts.h"
#include "ogr_p.h"
#include "cpl_conv.h"
#include "cpl_string.h"
#include <stdio.h>
#include <string>
#include <iostream>
#include <cstdio>
#include <sstream>
#include <vector>
#include <dirent.h>
#include <omp.h>
#include <sys/time.h>
#include <regex.h>
#include <math.h>
#include <hiredis/hiredis.h>
#include <queue>
#include <utility>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/foreach.hpp>
#include <boost/interprocess/managed_mapped_file.hpp>
#define MAX_NODE_SIZE 8
#define MAX_DOUBLE 100000000

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
namespace bgm = boost::geometry::model;
namespace bi = boost::interprocess;

using namespace std;
class Redis
{
	public:

    	Redis(){}
	
	~Redis()
	{
		this->_connect = NULL;
		this->_reply = NULL;	    	    
	}

	bool connect(char * host, int port)
	{
		this->_connect = redisConnect(host, port);
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

	void set(char * key, char * value)
	{
		redisCommand(this->_connect, "SET %s %s", key, value);
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
	
	void freeredis()
	{
		redisFree(this->_connect);
		//~ freeReplyObject(this->_reply);
	}
	
	private:

    	redisContext* _connect;
		redisReply* _reply;
				
};


void GetFilelist(char* argv, char* result[], char* flag, int& count)
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

void Usage()
{
    printf( "Usage:           [--shp:       input shapefile     ]\n"
	    "                     [--output:    output indexes      ]\n"
	    "                     [--rediskey:    rediskey      ]\n"
	    "                     [--redishost:    redishost        ]\n"
	    "                     [--redisport:    redisport        ]\n"
	    "                     [--maxsize:   max index file size(MB)  (default:100000MB) ]\n"		    	    
	    );
}

int main( int nArgc, char ** papszArgv )
{
	typedef bgm::d2::point_xy<double> point;
	typedef bgm::box<point> box;
	typedef bgm::segment<point> segment;
	typedef boost::tuple<segment,unsigned long,bool> polygon_segment;
	typedef std::pair<box,unsigned long> polygon_box;  
	 
    struct timeval t1,t2;
    double timeuse;
    gettimeofday(&t1,NULL);
	
	typedef bgi::quadratic<MAX_NODE_SIZE> params;
    typedef bgi::indexable<polygon_segment> indexable_segment;
    typedef bgi::equal_to<polygon_segment> equal_to_segment;
    typedef bi::allocator<polygon_segment, bi::managed_mapped_file::segment_manager> allocator_segment;
    typedef bgi::rtree<polygon_segment, params, indexable_segment, equal_to_segment, allocator_segment> rtree_segment;
	

	//~ typedef bgi::quadratic<MAX_NODE_SIZE> params;
    typedef bgi::indexable<polygon_box> indexable_box;
    typedef bgi::equal_to<polygon_box> equal_to_box;
    typedef bi::allocator<polygon_box, bi::managed_mapped_file::segment_manager> allocator_box;
    typedef bgi::rtree<polygon_box, params, indexable_box, equal_to_box, allocator_box> rtree_box;
	
	
	const char* shpFile = NULL;     
	const char* outIndex=NULL;
	char* redishost=NULL;
	char* rediskey=NULL;
	int redisport=0;
	char* outMBRIndex=new char[256];
	long size =100000;
	double tolerance=0.000001;//大于24级瓦片精度
	double minXOut=MAX_DOUBLE,minYOut=MAX_DOUBLE,maxXout=-1*MAX_DOUBLE,maxYOut=-1*MAX_DOUBLE; 
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
	CPLSetConfigOption("SHAPE_ENCODING", "UTF-8");
	for( int iArg = 1; iArg < nArgc; iArg++ )
	{
		if( EQUAL(papszArgv[iArg], "--shp") )
	    {
			shpFile = papszArgv[iArg+1];		
		}else if( EQUAL(papszArgv[iArg], "--output") )
	    {
			outIndex = papszArgv[iArg+1];
		}else if( EQUAL(papszArgv[iArg], "--redishost") )
	    {
			redishost = papszArgv[iArg+1];
		}else if( EQUAL(papszArgv[iArg], "--rediskey") )
	    {
			rediskey = papszArgv[iArg+1];
		}
		else if( EQUAL(papszArgv[iArg], "--redisport") )
	    {
			redisport = atoi(papszArgv[iArg+1]);
		}
		else if( EQUAL(papszArgv[iArg], "--maxsize") )
	    {
			size = atol(papszArgv[iArg+1]);
		}
	}
	
	if (shpFile == NULL || outIndex==NULL ||size==0||redishost==NULL||redisport==0)
	{
		Usage();
		return 0;
	}
	sprintf(outMBRIndex,"%s_mbr",outIndex);
	printf("--shp: %s\n", shpFile);
	printf("--output: %s  %s \n", outIndex,outMBRIndex);
	printf("--rediskey: %s\n", rediskey);
	printf("--redishost: %s\n", redishost);
	printf("--redisport: %d\n", redisport);
	printf("--maxsize: %ld MB\n", size);
	
	Redis *redis = new Redis();
	if(!redis->connect(redishost, redisport))
	{
		printf("connect redis error!\n");
		return 0;
	}

	remove(outIndex);
	remove(outMBRIndex);
	
	bi::managed_mapped_file file(bi::create_only, outIndex, size*1014*1024);
    allocator_segment alloc(file.get_segment_manager());
    rtree_segment * rtree_ptr = file.construct<rtree_segment>("rtree")(params(), indexable_segment(), equal_to_segment(), alloc);
	
	
	bi::managed_mapped_file filembr(bi::create_only, outMBRIndex, size*1014*1024);
    allocator_segment alloc_mbr(filembr.get_segment_manager());
    rtree_box * rtree_ptr_mbr = filembr.construct<rtree_box>("rtree")(params(), indexable_box(), equal_to_box(), alloc_mbr);
	

	std::cout << rtree_ptr->size()<<" "<<rtree_ptr_mbr->size() << std::endl;
	
	OGRRegisterAll();
	OGRLayer* shpLayer;
	OGRFeature *shpFeature;
	GDALDataset* shpDS =(GDALDataset*)GDALOpenEx(shpFile, GDAL_OF_VECTOR,NULL, NULL, NULL );
	if(shpDS==NULL)
	{
		printf("[ERROR] Open shpFile failed.\n");
		return 0;
	}
	shpLayer=shpDS->GetLayer(0);
	
	OGREnvelope env;
	shpLayer->GetExtent(&env);
	minXOut=minXOut<env.MinX? minXOut:env.MinX;
	minYOut=minYOut<env.MinY? minYOut:env.MinY;
	maxXout=maxXout>env.MaxX? maxXout:env.MaxX;
	maxYOut=maxYOut>env.MaxY? maxYOut:env.MaxY;
	shpLayer->ResetReading();
	
	
	OGRSpatialReference* fRef;
	OGRSpatialReference tRef;
	fRef = shpLayer->GetSpatialRef();
	tRef.importFromEPSG(3857);
	OGRCoordinateTransformation *coordTrans;
	coordTrans = OGRCreateCoordinateTransformation(fRef, &tRef);
	

	coordTrans ->Transform(1,&minXOut,&minYOut);
    coordTrans ->Transform(1,&maxXout,&maxYOut);
     char * tmp=new char[256];
    sprintf(tmp,"%lf,%lf,%lf,%lf", minXOut,minYOut,maxXout,maxYOut);
	redis->set(rediskey,tmp);
		printf("%s\n",tmp);	
	unsigned long polygonID=0;
	double px0,py0,px1,py1;
	double change;
	double minx,miny,maxx,maxy;
	while((shpFeature=shpLayer-> GetNextFeature())!= NULL)
	{	
		
		
		OGRGeometry *poGeometry=shpFeature->GetGeometryRef();
		int eType = wkbFlatten(poGeometry->getGeometryType());
		if(eType == wkbPolygon)
		{	
			
			OGRPolygon* pOGRPolygon=(OGRPolygon*) poGeometry;
			OGRLinearRing *pLinearRing = pOGRPolygon->getExteriorRing();
			pLinearRing->closeRings();
			int pointCount=pLinearRing->getNumPoints();
			px0 =pLinearRing->getX(0);
			py0 =pLinearRing->getY(0);
			change=py0-pLinearRing->getY(pointCount-2);
			coordTrans ->Transform(1,&px0,&py0);
			minx=px0;
			maxx=px0;
			miny=py0;
			maxy=py0;
			for(int i=1;i<pointCount;i++)
			{
				px1=pLinearRing->getX(i);
				py1=pLinearRing->getY(i);
				coordTrans ->Transform(1,&px1,&py1);
				minx=minx<px1? minx:px1;
				miny=miny<py1? miny:py1;
				maxx=maxx>px1? maxx:px1;
				maxy=maxy>py1? maxy:py1;
				if (py1==py0)
					rtree_ptr->insert(boost::make_tuple(segment(point(px0, py0),point(px1,py1)),polygonID,false));
				else if(change*(py0-py1)<0)
				{
					if(py0>py1)
						rtree_ptr->insert(boost::make_tuple(segment(point(px0-tolerance*(px1-px0)/(py1-py0), py0-tolerance),point(px1,py1)),polygonID,true));
					else
						rtree_ptr->insert(boost::make_tuple(segment(point(px0+tolerance*(px1-px0)/(py1-py0), py0+tolerance),point(px1,py1)),polygonID,true));
				}else
					rtree_ptr->insert(boost::make_tuple(segment(point(px0, py0),point(px1,py1)),polygonID,true));
				//~ printf(" %ld: %lf  %lf  %lf  %lf  \n",polygonID,px0,py0,px1,py1);
				px0=px1;
				py0=py1;
				change=py1-py0;
				
			}
	
			int innerCount=pOGRPolygon->getNumInteriorRings();
			for (int j=0;j<innerCount;j++)
			{
				pLinearRing = pOGRPolygon->getInteriorRing(j);
				pLinearRing->closeRings();
			    pointCount=pLinearRing->getNumPoints();
				px0 =pLinearRing->getX(0);
				py0 =pLinearRing->getY(0);
				change=py0-pLinearRing->getY(pointCount-2);
				coordTrans ->Transform(1,&px0,&py0);
				minx=minx<px0? minx:px0;
				miny=miny<py0? miny:py0;
				maxx=maxx>px0? maxx:px0;
				maxy=maxy>py0? maxy:py0;
				for(int i=1;i<pointCount;i++)
				{
					px1=pLinearRing->getX(i);
					py1=pLinearRing->getY(i);
					coordTrans ->Transform(1,&px1,&py1);
					minx=minx<px1? minx:px1;
					miny=miny<py1? miny:py1;
					maxx=maxx>px1? maxx:px1;
					maxy=maxy>py1? maxy:py1;
					if (py1==py0)
						rtree_ptr->insert(boost::make_tuple(segment(point(px0, py0),point(px1,py1)),polygonID,false));
					else if(change*(py0-py1)<0)
					{
						if(py0>py1)
							rtree_ptr->insert(boost::make_tuple(segment(point(px0-tolerance*(px1-px0)/(py1-py0), py0-tolerance),point(px1,py1)),polygonID,true));
						else
							rtree_ptr->insert(boost::make_tuple(segment(point(px0+tolerance*(px1-px0)/(py1-py0), py0+tolerance),point(px1,py1)),polygonID,true));
					}else
						rtree_ptr->insert(boost::make_tuple(segment(point(px0, py0),point(px1,py1)),polygonID,true));
					px0=px1;
					py0=py1;
					change=py1-py0;
				}
			}
			rtree_ptr_mbr->insert(std::make_pair(box(point(minx, miny),point(maxx,maxy)),polygonID));
			polygonID++;		
		}
		else if( eType == wkbMultiPolygon)
		{
			OGRMultiPolygon* pOGRMultiPolygon=(OGRMultiPolygon*) poGeometry;
			int polygonCount =pOGRMultiPolygon->getNumGeometries();
			for(int k=0;k<polygonCount;k++)
			{
				OGRPolygon* pOGRPolygon=(OGRPolygon*)pOGRMultiPolygon->getGeometryRef(k);
				OGRLinearRing *pLinearRing = pOGRPolygon->getExteriorRing();
				pLinearRing->closeRings();
				int pointCount=pLinearRing->getNumPoints();
				px0 =pLinearRing->getX(0);
				py0 =pLinearRing->getY(0);
				change=py0-pLinearRing->getY(pointCount-2);
				coordTrans ->Transform(1,&px0,&py0);
				minx=px0;
				maxx=px0;
				miny=py0;
				maxy=py0;
				for(int i=1;i<pointCount;i++)
				{
					px1=pLinearRing->getX(i);
					py1=pLinearRing->getY(i);
					coordTrans ->Transform(1,&px1,&py1);
					minx=minx<px1? minx:px1;
					miny=miny<py1? miny:py1;
					maxx=maxx>px1? maxx:px1;
					maxy=maxy>py1? maxy:py1;
					if (py1==py0)
						rtree_ptr->insert(boost::make_tuple(segment(point(px0, py0),point(px1,py1)),polygonID,false));
					else if(change*(py0-py1)<0)
					{
						if(py0>py1)
							rtree_ptr->insert(boost::make_tuple(segment(point(px0-tolerance*(px1-px0)/(py1-py0), py0-tolerance),point(px1,py1)),polygonID,true));
						else
							rtree_ptr->insert(boost::make_tuple(segment(point(px0+tolerance*(px1-px0)/(py1-py0), py0+tolerance),point(px1,py1)),polygonID,true));
					}else
						rtree_ptr->insert(boost::make_tuple(segment(point(px0, py0),point(px1,py1)),polygonID,true));
					px0=px1;
					py0=py1;
					change=py1-py0;
				}
		
				int innerCount=pOGRPolygon->getNumInteriorRings();
				for (int j=0;j<innerCount;j++)
				{
					pLinearRing = pOGRPolygon->getInteriorRing(j);
					pLinearRing->closeRings();
				    pointCount=pLinearRing->getNumPoints();
					px0 =pLinearRing->getX(0);
					py0 =pLinearRing->getY(0);
					change=py0-pLinearRing->getY(pointCount-2);
					coordTrans ->Transform(1,&px0,&py0);
					minx=minx<px0? minx:px0;
					miny=miny<py0? miny:py0;
					maxx=maxx>px0? maxx:px0;
					maxy=maxy>py0? maxy:py0;
					for(int i=1;i<pointCount;i++)
					{
						px1=pLinearRing->getX(i);
						py1=pLinearRing->getY(i);
						coordTrans ->Transform(1,&px1,&py1);
						minx=minx<px1? minx:px1;
						miny=miny<py1? miny:py1;
						maxx=maxx>px1? maxx:px1;
						maxy=maxy>py1? maxy:py1;
						if (py1==py0)
							rtree_ptr->insert(boost::make_tuple(segment(point(px0, py0),point(px1,py1)),polygonID,false));
						else if(change*(py0-py1)<0)
						{
							if(py0>py1)
								rtree_ptr->insert(boost::make_tuple(segment(point(px0-tolerance*(px1-px0)/(py1-py0), py0-tolerance),point(px1,py1)),polygonID,true));
							else
								rtree_ptr->insert(boost::make_tuple(segment(point(px0+tolerance*(px1-px0)/(py1-py0), py0+tolerance),point(px1,py1)),polygonID,true));
						}else
							rtree_ptr->insert(boost::make_tuple(segment(point(px0, py0),point(px1,py1)),polygonID,true));
						px0=px1;
						py0=py1;
						change=py1-py0;
					}
				}
				rtree_ptr_mbr->insert(std::make_pair(box(point(minx, miny),point(maxx,maxy)),polygonID));
				polygonID++;
			}		
		}
	}
	std::cout << rtree_ptr->size()<<" "<<rtree_ptr_mbr->size() << std::endl;
	bi::managed_mapped_file::shrink_to_fit(outIndex);
	bi::managed_mapped_file::shrink_to_fit(outMBRIndex);
	gettimeofday(&t2,NULL);
    timeuse=t2.tv_sec-t1.tv_sec+(t2.tv_usec-t1.tv_usec)/1000000.0;
    printf("Use Time:%f\n",timeuse);
}
