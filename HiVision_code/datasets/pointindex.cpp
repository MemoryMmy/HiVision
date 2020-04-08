/*
* 
*    data: 2018-08-10 
*    auther : Ma Mengyu@DBRG
*    description :
*    create indexes(Point objects)
*    dependency :
*    gdal2.10  boost1.64 
*    run :
* 	 ./Pointindex --shp ./alaska_OSM/gis_osm_places_free_1.shp --output ../indexes/pplaces --rediskey pplaces --redishost 127.0.0.1 --redisport 6379
*/


#include "ogrsf_frmts.h"
#include "ogr_p.h"
#include "cpl_conv.h"
#include "cpl_string.h"
#include <stdio.h>
#include <string>
#include <iostream>
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
#define MAX_DOUBLE 1000000000

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
namespace bgm = boost::geometry::model;
namespace bi = boost::interprocess;


using namespace std;

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
    
    struct timeval t1,t2;
    double timeuse;
    gettimeofday(&t1,NULL);

    typedef bgi::quadratic<MAX_NODE_SIZE> params_t;
    typedef bgi::indexable<point> indexable_t;
    typedef bgi::equal_to<point> equal_to_t;
    typedef bi::allocator<point, bi::managed_mapped_file::segment_manager> allocator_t;
    typedef bgi::rtree<point, params_t, indexable_t, equal_to_t, allocator_t> rtree_t;
    
	char* shpFile = NULL;     
	char* outIndex=NULL;
	char* redishost=NULL;
	char* rediskey=NULL;
	int redisport=0;
	long size =100000;
		
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
	printf("--shp: %s\n", shpFile);
	printf("--output: %s\n", outIndex);
	printf("--rediskey: %s\n", rediskey);
	printf("--redishost: %s\n", redishost);
	printf("--redisport: %d\n", redisport);
	printf("--maxsize: %ld MB\n", size);
	
	remove(outIndex);
	
	Redis *redis = new Redis();
	if(!redis->connect(redishost, redisport))
	{
		printf("connect redis error!\n");
		return 0;
	}
	
	bi::managed_mapped_file file(bi::create_only, outIndex, size*1024*1024);
	allocator_t alloc(file.get_segment_manager());
	rtree_t * rtree_ptr = file.find_or_construct<rtree_t>("rtree")(params_t(), indexable_t(), equal_to_t(), alloc);
	
	std::cout << rtree_ptr->size() << std::endl;

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
	

	double px0,py0;
	while((shpFeature=shpLayer-> GetNextFeature())!= NULL)
	{	
		OGRGeometry *poGeometry=shpFeature->GetGeometryRef();
		int eType = wkbFlatten(poGeometry->getGeometryType());
		if( eType == wkbPoint)
		{
			OGRPoint* pOGRPoint=(OGRPoint*) poGeometry;
			px0 = pOGRPoint->getX();
			py0 = pOGRPoint->getY();
			coordTrans ->Transform(1,&px0,&py0);
			rtree_ptr->insert(point(px0,py0));
		}else if( eType == wkbMultiPoint)
		{
			OGRMultiPoint* pOGRMultiPoint=(OGRMultiPoint*) poGeometry;
			int pointCount =pOGRMultiPoint->getNumGeometries();
			for(int j=0;j<pointCount;j++)
			{
				OGRPoint* pOGRPoint=(OGRPoint*)pOGRMultiPoint->getGeometryRef(j);
				px0 = pOGRPoint->getX();
				py0 = pOGRPoint->getY();
				coordTrans ->Transform(1,&px0,&py0);
				rtree_ptr->insert(point(px0,py0));
			}		
		}
	}
	std::cout << rtree_ptr->size() << std::endl;
	bi::managed_mapped_file::shrink_to_fit(outIndex);
	gettimeofday(&t2,NULL);
    timeuse=t2.tv_sec-t1.tv_sec+(t2.tv_usec-t1.tv_usec)/1000000.0;
    printf("Use Time:%f\n",timeuse);
}
