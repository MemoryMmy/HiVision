/*
* 
*    data: 2019-5-20 
*    auther : Ma Mengyu@National University of Defense Technology
*    e-mail : mamengyu10@nudt.edu.cn
*    description : 
*    Multi-Thread Visualization Server
*    run :
* 	 
*/

#include "ogrsf_frmts.h"
#include "ogr_p.h"
#include "crow.h"
#include "cpl_conv.h"
#include "cpl_string.h"
#include <stdio.h>
#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <dirent.h>
#include <omp.h>
#include <sys/time.h>
#include <regex.h>
#include <math.h>
#include <stdlib.h>
#include <png.h> 
#include <stdlib.h>
#include <stack>
#include <signal.h>


#include <hiredis/hiredis.h>
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

#define TILE_SIZE 256
#define L 20037508.34
#define MAX_LAYER_NUMBER 10
#define MAX_OPERATION_NUMBER 32


class HiVisionPara
{
	public:
	char *shpPath=NULL;
	char *indexPath=NULL;
	char *patternPath=NULL;
	char *pgPath=NULL;
	char *redisHost=NULL;
	int redisPort=6379;
	int servicePort=10080;
	
};
HiVisionPara hvpara;

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
class ServerLogHandler : public crow::ILogHandler {
    public:
        void log(std::string /*message*/, crow::LogLevel /*level*/) override {
//            cerr << "ServerLogHandler -> " << message;
        }
};
struct ServerMiddleware 
{
    std::string message;

    ServerMiddleware() 
    {
        message = "foo";
    }

    void setMessage(std::string newMsg)
    {
        message = newMsg;
    }

    struct context
    {
		
    };

    void before_handle(crow::request& /*req*/, crow::response& /*res*/, context& /*ctx*/)
    {
        CROW_LOG_DEBUG << " - MESSAGE: " << message;
    }

    void after_handle(crow::request& /*req*/, crow::response& /*res*/, context& /*ctx*/)
    {
        // no-op
    }
};
bool Linestring( char* shpFile, char* outIndex,char* redishost,char* rediskey,int redisport ,bool ispg=false)
{
	typedef bgm::d2::point_xy<double> point;
	//~ typedef bgm::box<point> box;
	typedef bgm::segment<point> segment;

    struct timeval t1,t2;
    double timeuse;
    gettimeofday(&t1,NULL);
	

    typedef bgi::quadratic<MAX_NODE_SIZE> params;
    typedef bgi::indexable<segment> indexable_segment;
    typedef bgi::equal_to<segment> equal_to_segment;
    typedef bi::allocator<segment, bi::managed_mapped_file::segment_manager> allocator_segment;
    typedef bgi::rtree<segment, params, indexable_segment, equal_to_segment, allocator_segment> rtree_segment;
	long size =300000;
		
	double minXOut=MAX_DOUBLE,minYOut=MAX_DOUBLE,maxXout=-1*MAX_DOUBLE,maxYOut=-1*MAX_DOUBLE; 
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
	CPLSetConfigOption("SHAPE_ENCODING", "UTF-8");
	
	fstream f;
	f.open(outIndex,ios::in);
	if (f)
	{
		printf("[ERROR] Dataset already registered! shpFile:%s outIndex:%s\n",shpFile,outIndex);
		return false;
	}
	f.close();
	
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
		printf("[ERROR] connect redis error! shpFile:%s outIndex:%s\n",shpFile,outIndex);
		return false;
	}
	
	GDALAllRegister();
	OGRRegisterAll();
	OGRLayer* shpLayer;
	OGRFeature *shpFeature;
	GDALDataset* shpDS;
	
	if(ispg)
	{
		shpDS =(GDALDataset*)GDALOpenEx(hvpara.pgPath, GDAL_OF_VECTOR,NULL, NULL, NULL );
		if(shpDS==NULL)
		{
			printf("[ERROR] Connect postgre failed. shpFile:%s outIndex:%s\n",shpFile,outIndex);
			return false;
		}
		shpLayer=shpDS->GetLayerByName(shpFile);
	}else
	{
		shpDS =(GDALDataset*)GDALOpenEx(shpFile, GDAL_OF_VECTOR,NULL, NULL, NULL );
		if(shpDS==NULL)
		{
			printf("[ERROR] Open shpFile failed. shpFile:%s outIndex:%s\n",shpFile,outIndex);
			return false;
		}
		shpLayer=shpDS->GetLayer(0);
	}
	if(shpLayer==NULL)
	{
		printf("[ERROR] Open shpFile failed. shpFile:%s outIndex:%s\n",shpFile,outIndex);
		return false;
	}
	
	OGREnvelope env;
	if(shpLayer->GetExtent(&env)!=0)
	{
		printf("[ERROR] Get extent failed. shpFile:%s outIndex:%s\n",shpFile,outIndex);
		return false;
	}
	minXOut=minXOut<env.MinX? minXOut:env.MinX;
	minYOut=minYOut<env.MinY? minYOut:env.MinY;
	maxXout=maxXout>env.MaxX? maxXout:env.MaxX;
	maxYOut=maxYOut>env.MaxY? maxYOut:env.MaxY;
	shpLayer->ResetReading();
	
	OGRSpatialReference* fRef;
	OGRSpatialReference tRef;
	fRef = shpLayer->GetSpatialRef();
	if(fRef==NULL)
	{
		printf("[ERROR] No SpatialRef. shpFile:%s outIndex:%s\n",shpFile,outIndex);
		return false;
	}
	tRef.importFromEPSG(3857);
	OGRCoordinateTransformation *coordTrans;
	coordTrans = OGRCreateCoordinateTransformation(fRef, &tRef);
	
	coordTrans ->Transform(1,&minXOut,&minYOut);
    coordTrans ->Transform(1,&maxXout,&maxYOut); 
    
    char * tmp=new char[256];
    sprintf(tmp,"%lf,%lf,%lf,%lf", minXOut,minYOut,maxXout,maxYOut);
	redis->set(rediskey,tmp);
	
	bi::managed_mapped_file file(bi::create_only, outIndex, size*1024*1024);
    allocator_segment alloc(file.get_segment_manager());
    rtree_segment * rtree_ptr = file.construct<rtree_segment>("rtree")(params(), indexable_segment(), equal_to_segment(), alloc);
	
	std::cout << rtree_ptr->size() << std::endl;
	double px0,py0,px1,py1;
	while((shpFeature=shpLayer-> GetNextFeature())!= NULL)
	{	
		if(shpFeature==NULL)
		continue;
		OGRGeometry *poGeometry=shpFeature->GetGeometryRef();
		if(poGeometry==NULL)
		continue;
		int eType = wkbFlatten(poGeometry->getGeometryType());
		if(eType == wkbLineString)
		{	
			OGRLineString* pOGRLineString=(OGRLineString*) poGeometry;
			int pointCount = pOGRLineString->getNumPoints();
			px0=pOGRLineString->getX(0);
			py0=pOGRLineString->getY(0);
			coordTrans ->Transform(1,&px0,&py0);
			for(int i=1;i<pointCount;i++)
			{
				px1=pOGRLineString->getX(i);
				py1=pOGRLineString->getY(i);
				coordTrans ->Transform(1,&px1,&py1);
				rtree_ptr->insert(segment(point(px0, py0),point(px1,py1)));
				px0=px1;
				py0=py1;
			}	
		}
		else if( eType == wkbMultiLineString)
		{
			OGRMultiLineString* pOGRMultiLineString=(OGRMultiLineString*) poGeometry;
			int lineCount =pOGRMultiLineString->getNumGeometries();
			for(int j=0;j<lineCount;j++)
			{
				OGRLineString* pOGRLineString=(OGRLineString*)pOGRMultiLineString->getGeometryRef(j);
				int pointCount = pOGRLineString->getNumPoints();
				px0=pOGRLineString->getX(0);
				py0=pOGRLineString->getY(0);
				coordTrans ->Transform(1,&px0,&py0);
				for(int i=1;i<pointCount;i++)
				{
					px1=pOGRLineString->getX(i);
					py1=pOGRLineString->getY(i);
					coordTrans ->Transform(1,&px1,&py1);
					rtree_ptr->insert(segment(point(px0, py0),point(px1,py1)));
					px0=px1;
					py0=py1;
				}
			}		
		}
	}
	if(rtree_ptr->size()<1)
	{
		printf("[ERROR] No feature found. shpFile:%s outIndex:%s\n",shpFile,outIndex);
		remove(outIndex);
		return false;
	}
	std::cout << rtree_ptr->size() << std::endl;
	bi::managed_mapped_file::shrink_to_fit(outIndex);
	gettimeofday(&t2,NULL);
    timeuse=t2.tv_sec-t1.tv_sec+(t2.tv_usec-t1.tv_usec)/1000000.0;
     printf("[DONE] %s Use Time:%f\n", outIndex,timeuse);
	
    return true;
}
bool Point( char* shpFile, char* outIndex,char* redishost,char* rediskey,int redisport  ,bool ispg=false)
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
    
	long size =300000;
		
	double minXOut=MAX_DOUBLE,minYOut=MAX_DOUBLE,maxXout=-1*MAX_DOUBLE,maxYOut=-1*MAX_DOUBLE; 
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
	CPLSetConfigOption("SHAPE_ENCODING", "UTF-8");
	
	fstream f;
	f.open(outIndex,ios::in);
	if (f)
	{
		printf("[ERROR] Dataset already registered! shpFile:%s outIndex:%s\n",shpFile,outIndex);
		return false;
	}
	f.close();
	
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
		printf("[ERROR] connect redis error!\n");
		return false;
	}
	
	OGRRegisterAll();
	OGRLayer* shpLayer;
	OGRFeature *shpFeature;
	GDALDataset* shpDS;
	if(ispg)
	{
		shpDS =(GDALDataset*)GDALOpenEx(hvpara.pgPath, GDAL_OF_VECTOR,NULL, NULL, NULL );
		if(shpDS==NULL)
		{
			printf("[ERROR] Connect postgre failed. shpFile:%s outIndex:%s\n",shpFile,outIndex);
			return false;
		}
		shpLayer=shpDS->GetLayerByName(shpFile);
	}else
	{
		shpDS =(GDALDataset*)GDALOpenEx(shpFile, GDAL_OF_VECTOR,NULL, NULL, NULL );
		if(shpDS==NULL)
		{
			printf("[ERROR] Open shpFile failed. shpFile:%s outIndex:%s\n",shpFile,outIndex);
			return false;
		}
		shpLayer=shpDS->GetLayer(0);
	}
	if(shpLayer==NULL)
	{
		printf("[ERROR] Open shpFile failed. shpFile:%s outIndex:%s\n",shpFile,outIndex);
		return false;
	}

	
	
	OGREnvelope env;
	if(shpLayer->GetExtent(&env)!=0)
	{
		printf("[ERROR] Get extent failed. shpFile:%s outIndex:%s\n",shpFile,outIndex);
		return false;
	}
	minXOut=minXOut<env.MinX? minXOut:env.MinX;
	minYOut=minYOut<env.MinY? minYOut:env.MinY;
	maxXout=maxXout>env.MaxX? maxXout:env.MaxX;
	maxYOut=maxYOut>env.MaxY? maxYOut:env.MaxY;
	shpLayer->ResetReading();

	OGRSpatialReference* fRef;
	OGRSpatialReference tRef;
	fRef = shpLayer->GetSpatialRef();
	if(fRef==NULL)
	{
		printf("[ERROR] No SpatialRef. shpFile:%s outIndex:%s\n",shpFile,outIndex);
		return false;
	}
	tRef.importFromEPSG(3857);
	OGRCoordinateTransformation *coordTrans;
	coordTrans = OGRCreateCoordinateTransformation(fRef, &tRef);
	
    coordTrans ->Transform(1,&minXOut,&minYOut);
    coordTrans ->Transform(1,&maxXout,&maxYOut);
    char * tmp=new char[256];
    sprintf(tmp,"%lf,%lf,%lf,%lf", minXOut,minYOut,maxXout,maxYOut);
	redis->set(rediskey,tmp);
	bi::managed_mapped_file file(bi::create_only, outIndex, size*1024*1024);
	allocator_t alloc(file.get_segment_manager());
	rtree_t * rtree_ptr = file.find_or_construct<rtree_t>("rtree")(params_t(), indexable_t(), equal_to_t(), alloc);
	
	std::cout << rtree_ptr->size() << std::endl;

	double px0,py0;
	while((shpFeature=shpLayer-> GetNextFeature())!= NULL)
	{	
		if(shpFeature==NULL)
		continue;
		OGRGeometry *poGeometry=shpFeature->GetGeometryRef();
		if(poGeometry==NULL)
		continue;
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
	if(rtree_ptr->size()<1)
	{
		printf("[ERROR] No feature found. shpFile:%s outIndex:%s\n",shpFile,outIndex);
		remove(outIndex);
		return false;
	}
	std::cout << rtree_ptr->size() << std::endl;
	bi::managed_mapped_file::shrink_to_fit(outIndex);
	gettimeofday(&t2,NULL);
    timeuse=t2.tv_sec-t1.tv_sec+(t2.tv_usec-t1.tv_usec)/1000000.0;
    printf("[DONE] %s Use Time:%f\n", outIndex,timeuse);
    return true;
}
bool Polygon( char* shpFile, char* outIndex,char* redishost,char* rediskey,int redisport,bool ispg=false)
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
	
	char* outMBRIndex=new char[256];
	long size =300000;
	double tolerance=0.00001;//大于24级瓦片精度
	double minXOut=MAX_DOUBLE,minYOut=MAX_DOUBLE,maxXout=-1*MAX_DOUBLE,maxYOut=-1*MAX_DOUBLE; 
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
	CPLSetConfigOption("SHAPE_ENCODING", "UTF-8");
	
	
	fstream f;
	f.open(outIndex,ios::in);
	if (f)
	{
		printf("[ERROR] Dataset already registered! shpFile:%s outIndex:%s\n",shpFile,outIndex);
		return false;
	}
	f.close();

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
		printf("[ERROR] connect redis error! shpFile:%s outIndex:%s\n",shpFile,outIndex);
		return false;
	}

	remove(outIndex);
	remove(outMBRIndex);
	
	OGRRegisterAll();
	OGRLayer* shpLayer;
	OGRFeature *shpFeature;
	GDALDataset* shpDS;
	if(ispg)
	{
		shpDS =(GDALDataset*)GDALOpenEx(hvpara.pgPath, GDAL_OF_VECTOR,NULL, NULL, NULL );
		if(shpDS==NULL)
		{
			printf("[ERROR] Connect postgre failed. shpFile:%s outIndex:%s\n",shpFile,outIndex);
			return false;
		}
		shpLayer=shpDS->GetLayerByName(shpFile);
	}else
	{
		shpDS =(GDALDataset*)GDALOpenEx(shpFile, GDAL_OF_VECTOR,NULL, NULL, NULL );
		if(shpDS==NULL)
		{
			printf("[ERROR] Open shpFile failed. shpFile:%s outIndex:%s\n",shpFile,outIndex);
			return false;
		}
		shpLayer=shpDS->GetLayer(0);
	}
	if(shpLayer==NULL)
	{
		printf("[ERROR] No such layer. shpFile:%s outIndex:%s\n",shpFile,outIndex);
		return false;
	}
	OGREnvelope env;
	if(shpLayer->GetExtent(&env)!=0)
	{
		printf("[ERROR] Get extent failed. shpFile:%s outIndex:%s\n",shpFile,outIndex);
		return false;
	}
	minXOut=minXOut<env.MinX? minXOut:env.MinX;
	minYOut=minYOut<env.MinY? minYOut:env.MinY;
	maxXout=maxXout>env.MaxX? maxXout:env.MaxX;
	maxYOut=maxYOut>env.MaxY? maxYOut:env.MaxY;
	shpLayer->ResetReading();
	OGRSpatialReference* fRef;
	OGRSpatialReference tRef;
	fRef = shpLayer->GetSpatialRef();
	if(fRef==NULL)
	{
		printf("[ERROR] No SpatialRef. shpFile:%s outIndex:%s\n",shpFile,outIndex);
		return false;
	}

	tRef.importFromEPSG(3857);
	OGRCoordinateTransformation *coordTrans;
	coordTrans = OGRCreateCoordinateTransformation(fRef, &tRef);
	
	coordTrans ->Transform(1,&minXOut,&minYOut);
    coordTrans ->Transform(1,&maxXout,&maxYOut);
    char * tmp=new char[256];
    sprintf(tmp,"%lf,%lf,%lf,%lf", minXOut,minYOut,maxXout,maxYOut);
	redis->set(rediskey,tmp);
	bi::managed_mapped_file file(bi::create_only, outIndex, size*1024*1024);
    allocator_segment alloc(file.get_segment_manager());
    rtree_segment * rtree_ptr = file.construct<rtree_segment>("rtree")(params(), indexable_segment(), equal_to_segment(), alloc);
	
	bi::managed_mapped_file filembr(bi::create_only, outMBRIndex, size*128*1024);
    allocator_segment alloc_mbr(filembr.get_segment_manager());
    rtree_box * rtree_ptr_mbr = filembr.construct<rtree_box>("rtree")(params(), indexable_box(), equal_to_box(), alloc_mbr);
	std::cout << rtree_ptr->size()<<" "<<rtree_ptr_mbr->size() << std::endl;
	unsigned long polygonID=0;
	double px0,py0,px1,py1;
	double change;
	double minx,miny,maxx,maxy;
	while((shpFeature=shpLayer-> GetNextFeature())!= NULL)
	{	
		if(shpFeature==NULL)
		continue;
		OGRGeometry *poGeometry=shpFeature->GetGeometryRef();
		if(poGeometry==NULL)
		continue;
		
		//~ printf("polygonID:%ld\n",polygonID);
		int eType = wkbFlatten(poGeometry->getGeometryType());
		if(eType == wkbPolygon)
		{	
			//~ printf("wkbPolygon: polygonID %ld \n",polygonID);
			OGRPolygon* pOGRPolygon=(OGRPolygon*) poGeometry;
			OGRLinearRing *pLinearRing = pOGRPolygon->getExteriorRing();
			if(pLinearRing==NULL)
			{
				continue;
			}
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
				if (fabs(py1-py0)<tolerance)
				{
					//~ printf("py1:%lf py0:%lf abs(py1-py0):%20.18lf tolerance:%20.18lf\n",py1,py0,abs(py1-py0),tolerance);
					rtree_ptr->insert(boost::make_tuple(segment(point(px0, py0),point(px1,py1)),polygonID,false));
				}
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
					if (fabs(py1-py0)<tolerance)
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
			}
			rtree_ptr_mbr->insert(std::make_pair(box(point(minx, miny),point(maxx,maxy)),polygonID));
			polygonID++;		
		}
		else if( eType == wkbMultiPolygon)
		{
			//~ printf("wkbMultiPolygon: polygonID %ld \n",polygonID);
			OGRMultiPolygon* pOGRMultiPolygon=(OGRMultiPolygon*) poGeometry;
			int polygonCount =pOGRMultiPolygon->getNumGeometries();
			for(int k=0;k<polygonCount;k++)
			{
				OGRPolygon* pOGRPolygon=(OGRPolygon*)pOGRMultiPolygon->getGeometryRef(k);
				OGRLinearRing *pLinearRing = pOGRPolygon->getExteriorRing();
				if(pLinearRing==NULL)
				{
					continue;
				}
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
					if (fabs(py1-py0)<tolerance)
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
					//~ pLinearRing->closeRings();
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
						if (fabs(py1-py0)<tolerance)
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
				}
				rtree_ptr_mbr->insert(std::make_pair(box(point(minx, miny),point(maxx,maxy)),polygonID));
				polygonID++;
			}		
		}
	}
	if(rtree_ptr->size()<1)
	{
		printf("[ERROR] No feature found. shpFile:%s outIndex:%s\n",shpFile,outIndex);
		remove(outIndex);
		remove(outMBRIndex);
		return false;
	}
	std::cout << rtree_ptr->size()<<" "<<rtree_ptr_mbr->size() << std::endl;
	bi::managed_mapped_file::shrink_to_fit(outIndex);
	bi::managed_mapped_file::shrink_to_fit(outMBRIndex);
	gettimeofday(&t2,NULL);
    timeuse=t2.tv_sec-t1.tv_sec+(t2.tv_usec-t1.tv_usec)/1000000.0;
    printf("[DONE] %s Use Time:%f\n", outIndex,timeuse);
    return true;
}


int main(int nArgc, char ** papszArgv)
{
	hvpara.shpPath=papszArgv[1];
	hvpara.indexPath=papszArgv[2];
	hvpara.patternPath=papszArgv[3];
	hvpara.redisHost=papszArgv[4];
	hvpara.redisPort=atoi(papszArgv[5]);
	hvpara.servicePort=atoi(papszArgv[6]);
    
    crow::App<ServerMiddleware> app;
    CROW_ROUTE(app,"/HiVision/add/<string>/<string>/<int>").name("HiVisionAdd")
    ([](const crow::request& req, crow::response& res, string shp, string appid, int type){
        char* redisHost=hvpara.redisHost;
        int redisPort=hvpara.redisPort;
        char* shpPath=hvpara.shpPath;
        char* indexPath=hvpara.indexPath;    
        try{
			char * tmp1=new char[256];
			char * tmp2=new char[256];
			char * tmp3=new char[256];
			if(type<3)
			{
				sprintf(tmp1,"%s/%s.shp",shpPath,shp.c_str());
				fstream f;
				f.open(tmp1,ios::in);
				if (!f)
				{
					res.write("[ERROR] Dataset not found!");
					throw "Dataset not found!";
				}
				f.close();
			}else
			{
				sprintf(tmp1,"%s",shp.c_str());
			}
			if(type==0)
			{
				sprintf(tmp2,"%sp%s",indexPath,appid.c_str());
				sprintf(tmp3,"p%s",appid.c_str());
				if (Point(tmp1,tmp2,redisHost,tmp3,redisPort))
					res.write("[DONE] Register done!");
				else 
					throw "Failed";
			}else if(type==1)
			{
				sprintf(tmp2,"%sl%s",indexPath,appid.c_str());
				sprintf(tmp3,"l%s",appid.c_str());
				if (Linestring(tmp1,tmp2,redisHost,tmp3,redisPort))
					res.write("[DONE] Register done!");
				else 
					throw "Failed";
			}else if(type==2)
			{
				sprintf(tmp2,"%sa%s",indexPath,appid.c_str());
				sprintf(tmp3,"a%s",appid.c_str());
				if (Polygon(tmp1,tmp2,redisHost,tmp3,redisPort))
					res.write("[DONE] Register done!");
				else 
					throw "Failed";
			}else
			{
				res.write("[ERROR] Wrong type!");
			}			
			res.end();
		}
        catch(const char* msg)
		{	
			res.write("[ERROR] Register failed");
			res.end();
		}
    });

    CROW_ROUTE(app, "/HiVision/<string>/<int>/<int>/<int>/<double>/<int>/<int>/<int>.png").name("HiVision")
    ([](const crow::request& req, crow::response& res,string appid, int R, int G, int B, double AD,int z,int x, int y){
        char* redisHost=hvpara.redisHost;
        int redisPort=hvpara.redisPort;
        std::ostringstream os;
        Redis *redis = new Redis();
		if(!redis->connect(redisHost, redisPort))
		{
			printf("connect redis error!\n");
		}
		int A=AD*64;
		int AH=AD*256;
		int Ax=(1-AD)*64;
		int Ay=-64;
		
		char* key=new char[32];
		char* newkey= new char[32];
		
		vector<char> pos;
		long size;
		//~ int layercount;
		png_bytep * row_pointers=(png_bytep*)malloc(256*sizeof(png_bytep));
		for(int i = 0; i < TILE_SIZE; i++)
			row_pointers[i] = (png_bytep)malloc(1024);	
        try{
			char * tmp=new char[256];
			int count;
			char* bbox[8];
			sprintf(tmp,"%s",redis->get(appid).c_str());
			GetList(tmp, bbox, (char*)",", count);
			//~ double r=atof(radius[0]);
			double minx=atof(bbox[0]);
			double miny=atof(bbox[1]);
			double maxx=atof(bbox[2]);
			double maxy=atof(bbox[3]);
			double tile_minx = ((256*x+0.5)/(128<<z)-1)* L;
			double tile_miny = (1-(256*y+255.5)/(128<<z))* L;
			double tile_maxx = ((256*x+255.5)/(128<<z)-1)* L;
			double tile_maxy = (1-(256*y-0.5)/(128<<z))* L;
			delete tmp;
			if (tile_minx<maxx && tile_miny<maxy && tile_maxx>minx && tile_maxy> miny)
			{
				sprintf(key,"%s/%d/%d/%d",appid.c_str(),z,x,y);
				char (*buffer_area)[TILE_SIZE] = (char(*)[TILE_SIZE])malloc(TILE_SIZE*TILE_SIZE);
				if (! redis->zget(key,(char *)buffer_area))
				{
					memset(newkey, 0, 32);
					redisContext* rc = redisConnect(redisHost, redisPort);
					redisReply* reply = (redisReply*) redisCommand(rc, "subscribe HiVisiontiles");
					redis->lpush("HiVisiontasklist",key);
					while (strcmp(newkey,key)!=0)
					{
						if(redisGetReply(rc, (void **)&reply) == REDIS_OK&& reply->type == REDIS_REPLY_ARRAY)
							sprintf(newkey,"%s", reply->element[2]->str);
					}
					redis->zget(key,(char *)buffer_area);
					freeReplyObject(reply);
					redisFree(rc);	
				}
				png_structp png_ptr=png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
				png_infop info_ptr = png_create_info_struct(png_ptr);
				FILE *temp_png=tmpfile();
				png_init_io(png_ptr, temp_png);
				png_set_IHDR(png_ptr, info_ptr, TILE_SIZE, TILE_SIZE, 8, PNG_COLOR_TYPE_RGB_ALPHA, PNG_INTERLACE_NONE,PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
				png_write_info(png_ptr, info_ptr);	
				for(int i = 0; i < TILE_SIZE; i++)
					for(int j = 0; j < TILE_SIZE; j++)
					{
						if (buffer_area[i][j]>4)
						{
							row_pointers[i][4*j]  = R;// red
							row_pointers[i][4*j+1] = G;// green
							row_pointers[i][4*j+2] = B;// blue
							row_pointers[i][4*j+3] = AH+Ax*(buffer_area[i][j]-4)-1;
						}
						else if (buffer_area[i][j]>0)
						{
							
							row_pointers[i][4*j]  = R;// red
							row_pointers[i][4*j+1] = G;// green
							row_pointers[i][4*j+2] = B;// blue
							row_pointers[i][4*j+3] = buffer_area[i][j]*A-1;
						}
						else if (buffer_area[i][j]<0)
						{
							row_pointers[i][4*j]  = R;// red
							row_pointers[i][4*j+1] = G;// green
							row_pointers[i][4*j+2] = B;// blue
							row_pointers[i][4*j+3] = buffer_area[i][j]*Ay-1;
							
						}
						else
						{
							row_pointers[i][4*j]  = 0;// red
							row_pointers[i][4*j+1] = 0;// green
							row_pointers[i][4*j+2] = 0;// blue
							row_pointers[i][4*j+3] = 0;
						}
					}
				free(buffer_area);
				png_write_image(png_ptr, row_pointers);
				png_write_end(png_ptr, NULL);
				fseek(temp_png ,0,SEEK_END);
				size=ftell(temp_png);
				rewind(temp_png);
				pos.resize(size);
				fread(&pos[0], 1, size,temp_png);
				fclose(temp_png);
				string pos_tostr=string(pos.begin(),pos.end());
				res.write(pos_tostr);
				res.set_header("Content-Type", "image/png");	
				res.set_header("Access-Control-Allow-Origin", "*");				
			}
			else
				throw "Tile out range!";	
			redis->freeredis();
			delete key;
			delete newkey;
			for (int j=0; j<256; j++)
				free(row_pointers[j]);
			free(row_pointers);
			res.end();				
        }
        catch(const char* msg)
		{	
			redis->freeredis();	
			delete key;
			delete newkey;

			for (int j=0; j<256; j++)
				free(row_pointers[j]);
			free(row_pointers);
			res.code=500;
			os<<"ERROR: " << msg<< "\n";
			res.write(os.str());
			res.set_header("Content-Type", "text/html");
			res.end();
		}
    });
    

    CROW_ROUTE(app, "/HiVision/<string>/<int>/<int>/<int>/<string>/<double>/<int>/<int>/<int>.png").name("HiVision")
    ([](const crow::request& req, crow::response& res,string appid, int R, int G, int B, string pattern, double AD,int z,int x, int y){
        char* redisHost=hvpara.redisHost;
        int redisPort=hvpara.redisPort;
        std::ostringstream os;
        Redis *redis = new Redis();
		if(!redis->connect(redisHost, redisPort))
		{
			printf("connect redis error!\n");
		}
		int A=AD*64;
		int AH=AD*256;
		int Ax=(1-AD)*64;
		int Ay=-64;
		
		char* key=new char[32];
		char* newkey= new char[32];
		
		vector<char> pos;
		long size;
		
		
		FILE *pic_fp;
		png_structp png_ptr,pattern_png_ptr;
		png_infop  info_ptr,pattern_info_ptr;
		
		pattern_png_ptr  = png_create_read_struct(PNG_LIBPNG_VER_STRING, 0, 0, 0);
		pattern_info_ptr = png_create_info_struct(pattern_png_ptr);
		setjmp(png_jmpbuf(pattern_png_ptr)); 
			
		png_ptr  = png_create_read_struct(PNG_LIBPNG_VER_STRING, 0, 0, 0);
		info_ptr = png_create_info_struct(png_ptr);
		setjmp(png_jmpbuf(png_ptr)); 
		
		png_bytep * row_pointers=(png_bytep*)malloc(256*sizeof(png_bytep));
		for(int i = 0; i < TILE_SIZE; i++)
			row_pointers[i] = (png_bytep)malloc(1024);	
        try{
			
			char * pngfile= new char[256];
			sprintf(pngfile,"%s%s.png",hvpara.patternPath,pattern.c_str());
			pic_fp = fopen(pngfile, "rb");
			if(pic_fp == NULL)
				throw "Open pattern file failed!";
			
			rewind(pic_fp);
			png_init_io(pattern_png_ptr, pic_fp);
			png_read_png(pattern_png_ptr, pattern_info_ptr, PNG_TRANSFORM_EXPAND, 0);

			int pattern_width, pattern_height;

			png_bytep* pattern_row_pointers;
			pattern_row_pointers = png_get_rows(pattern_png_ptr, pattern_info_ptr);
			pattern_width = png_get_image_width(pattern_png_ptr, pattern_info_ptr);
			pattern_height = png_get_image_height(pattern_png_ptr, pattern_info_ptr);		
			
			char * tmp=new char[256];
			int count;
			char* bbox[8];
			sprintf(tmp,"%s",redis->get(appid).c_str());
			GetList(tmp, bbox, (char*)",", count);

			double minx=atof(bbox[0]);
			double miny=atof(bbox[1]);
			double maxx=atof(bbox[2]);
			double maxy=atof(bbox[3]);
			double tile_minx = ((256*x+0.5)/(128<<z)-1)* L;
			double tile_miny = (1-(256*y+255.5)/(128<<z))* L;
			double tile_maxx = ((256*x+255.5)/(128<<z)-1)* L;
			double tile_maxy = (1-(256*y-0.5)/(128<<z))* L;
			delete tmp;
			if (tile_minx<maxx && tile_miny<maxy && tile_maxx>minx && tile_maxy> miny)
			{
				sprintf(key,"%s/%d/%d/%d",appid.c_str(),z,x,y);
				char (*buffer_area)[TILE_SIZE] = (char(*)[TILE_SIZE])malloc(TILE_SIZE*TILE_SIZE);
				if (! redis->zget(key,(char *)buffer_area))
				{
					memset(newkey, 0, 32);
					redisContext* rc = redisConnect(redisHost, redisPort);
					redisReply* reply = (redisReply*) redisCommand(rc, "subscribe HiVisiontiles");
					redis->lpush("HiVisiontasklist",key);
					while (strcmp(newkey,key)!=0)
					{
						if(redisGetReply(rc, (void **)&reply) == REDIS_OK&& reply->type == REDIS_REPLY_ARRAY)
							sprintf(newkey,"%s", reply->element[2]->str);
					}
					redis->zget(key,(char *)buffer_area);
					freeReplyObject(reply);
					redisFree(rc);	
				}
				png_structp png_ptr=png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
				png_infop info_ptr = png_create_info_struct(png_ptr);
				FILE *temp_png=tmpfile();
				png_init_io(png_ptr, temp_png);
				png_set_IHDR(png_ptr, info_ptr, TILE_SIZE, TILE_SIZE, 8, PNG_COLOR_TYPE_RGB_ALPHA, PNG_INTERLACE_NONE,PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
				png_write_info(png_ptr, info_ptr);	
				for(int i = 0; i < TILE_SIZE; i++)
					for(int j = 0; j < TILE_SIZE; j++)
					{
						if (buffer_area[i][j]>4)
						{
							
							row_pointers[i][4*j]  = R;
							row_pointers[i][4*j+1] = G;
							row_pointers[i][4*j+2] = B;
							row_pointers[i][4*j+3] = AH+Ax*(buffer_area[i][j]-4)-1;
						}
						else if (buffer_area[i][j]==4)	
						{
							row_pointers[i][4*j]  = pattern_row_pointers[(i+256*y)%pattern_height][4*((j+256*x)%pattern_width)];
							row_pointers[i][4*j+1] = pattern_row_pointers[(i+256*y)%pattern_height][4*((j+256*x)%pattern_width)+1];
							row_pointers[i][4*j+2] = pattern_row_pointers[(i+256*y)%pattern_height][4*((j+256*x)%pattern_width)+2];
							row_pointers[i][4*j+3] = pattern_row_pointers[(i+256*y)%pattern_height][4*((j+256*x)%pattern_width)+3];;
						}
						else if (buffer_area[i][j]>0)
						{
							
							row_pointers[i][4*j]  = R;// red
							row_pointers[i][4*j+1] = G;// green
							row_pointers[i][4*j+2] = B;// blue
							row_pointers[i][4*j+3] = buffer_area[i][j]*A-1;
						}
						else if (buffer_area[i][j]<0)
						{
							row_pointers[i][4*j]  = R;// red
							row_pointers[i][4*j+1] = G;// green
							row_pointers[i][4*j+2] = B;// blue
							row_pointers[i][4*j+3] = buffer_area[i][j]*Ay-1;
							
						}
						else
						{
							row_pointers[i][4*j]  = 0;// red
							row_pointers[i][4*j+1] = 0;// green
							row_pointers[i][4*j+2] = 0;// blue
							row_pointers[i][4*j+3] = 0;
						}
					}
				free(buffer_area);
				png_write_image(png_ptr, row_pointers);
				png_write_end(png_ptr, NULL);
				fseek(temp_png ,0,SEEK_END);
				size=ftell(temp_png);
				rewind(temp_png);
				pos.resize(size);
				fread(&pos[0], 1, size,temp_png);
				fclose(temp_png);
				string pos_tostr=string(pos.begin(),pos.end());
				res.write(pos_tostr);
				res.set_header("Content-Type", "image/png");	
				res.set_header("Access-Control-Allow-Origin", "*");				
			}
			else
				throw "Tile out range!";	
			redis->freeredis();
			delete key;
			delete newkey;
			png_destroy_read_struct(&png_ptr, &info_ptr, 0); 
			png_destroy_read_struct(&pattern_png_ptr, &pattern_info_ptr, 0); 
			fclose(pic_fp);
			for (int j=0; j<256; j++)
				free(row_pointers[j]);
			free(row_pointers);
			res.end();				
        }
        catch(const char* msg)
		{	
			redis->freeredis();	
			delete key;
			delete newkey;
			png_destroy_read_struct(&png_ptr, &info_ptr, 0); 
			fclose(pic_fp);
			for (int j=0; j<256; j++)
				free(row_pointers[j]);
			free(row_pointers);
			res.code=500;
			os<<"ERROR: " << msg<< "\n";
			res.write(os.str());
			res.set_header("Content-Type", "text/html");
			res.end();
		}
    });
    
    // enables all log
    //~ app.loglevel(crow::LogLevel::DEBUG);
    //crow::logger::setHandler(std::make_shared<ServerLogHandler>());
    app.port(hvpara.servicePort)
        .multithreaded()
        .run();
}
