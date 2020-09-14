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

#include <ogrsf_frmts.h>
#include <ogr_p.h>
#include <cpl_conv.h>
#include <cpl_string.h>
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

#define MAX_NODE_SIZE	8
#define MAX_TILE_PARAMS 256
#define TILE_SIZE	256
#define L		20037508.34

namespace bg	= boost::geometry;
namespace bgi	= boost::geometry::index;
namespace bgm	= boost::geometry::model;
namespace bi	= boost::interprocess;

using namespace std;

typedef bgm::d2::point_xy<double>									point;
typedef bgm::box<point>											box;
typedef bgm::segment<point>										segment;
typedef boost::tuple<segment, unsigned long, bool>							polygon_segment;
typedef std::pair<box, unsigned long>									polygon_box;
typedef bgi::quadratic<MAX_NODE_SIZE>									params;
typedef bgi::indexable<point>										indexable_point;
typedef bgi::equal_to<point>										equal_to_point;
typedef bi::allocator<point, bi::managed_mapped_file::segment_manager>					allocator_point;
typedef bgi::rtree<point, params, indexable_point, equal_to_point, allocator_point>			rtree_point;
typedef bgi::indexable<segment>										indexable_segment;
typedef bgi::equal_to<segment>										equal_to_segment;
typedef bi::allocator<segment, bi::managed_mapped_file::segment_manager>				allocator_segment;
typedef bgi::rtree<segment, params, indexable_segment, equal_to_segment, allocator_segment>		rtree_segment;
typedef bgi::indexable<polygon_segment>									indexable_psegment;
typedef bgi::equal_to<polygon_segment>									equal_to_psegment;
typedef bi::allocator<polygon_segment, bi::managed_mapped_file::segment_manager>			allocator_psegment;
typedef bgi::rtree<polygon_segment, params, indexable_psegment, equal_to_psegment, allocator_psegment>	rtree_psegment;
typedef bgi::indexable<polygon_box>									indexable_box;
typedef bgi::equal_to<polygon_box>									equal_to_box;
typedef bi::allocator<polygon_box, bi::managed_mapped_file::segment_manager>				allocator_box;
typedef bgi::rtree<polygon_box, params, indexable_box, equal_to_box, allocator_box>			rtree_box;


/**
 * @brief Compare x span of two polygons
 * @param a     Input polygon
 * @param b     Input polygon
 *
 * @return              If the x span of polygon a is smaller than than the x span of polygon b, return true; otherwise, return true.
 */
bool SortPolygon( polygon_box a, polygon_box b )
{
	return( (bg::get<1, 0>( a.first ) - bg::get<0, 0>( a.first ) ) < (bg::get<1, 0>( b.first ) - bg::get<0, 0>( b.first ) ) );
}


/**
 * @brief Split a string into a list where each word is a list item
 * @param argv          Input string
 * @param result[]  Output Word list
 * @param flag          Separator used to split the string
 * @param count         Output Word list length
 *
 */
void GetList( char* argv, char* result[], char* flag, int & count )
{
	char	* string = strdup( argv );
	char	* p;
	int	i = 0;
	while ( (p = strsep( &string, flag ) ) != NULL )
	{
		result[i] = p;
		i++;
	}
	result[i]	= string;
	count		= i;
}


/**
 * @brief Visualize Linestring Objects
 * @param z                     Zoom level of the tile
 * @param x                     x coordinate of the tile
 * @param y                     y coordinate of the tile
 * @param indexPath     Spatial index path
 * @param dataId                Dataset id
 * @param tile_area     Matrix indentifing different regions(e.g. background region, color transition region, the zones of rasterized spatial objects)
 *
 */
void  LineVision( int z, int x, int y, char* indexPath, char *dataId, char *tile_area )
{
	char		* tmp = new char[256];
	double		Rz, R, R1, R2, r, Rz4;
	rtree_segment	* rtree_segment_ptr;
	Rz	= L / (128 << z);       /* Calculate the pixel resolution of zoom level z */
	R	= 2 * Rz;               /* R=2*Rz ,we set the visualizaiton radius N to 2 (pixels) */
	R1	= 1.6464 * Rz;          /* R1=R-sqrt(2)/4*Rz=(2-sqrt(2)/4)*Rz */
	R2	= 2.3536 * Rz;          /* R2=R+sqrt(2)/4*Rz=(2+sqrt(2)/4)*Rz */
	r	= 0.7071 * R1;          /* r=sqrt(2)/2*R1 */
	Rz4	= 0.25 * Rz;            /* Rz4=1/4*Rz */
	sprintf( tmp, "%s%s", indexPath, dataId );
	bi::managed_mapped_file file( bi::open_only, tmp );
	rtree_segment_ptr = file.find<rtree_segment>( "rtree" ).first;
	#pragma omp parallel for num_threads(2)
	for ( int i = 0; i < 256; i++ )
	{
		for ( int j = 0; j < 256; j++ )
		{
			int tile_index = i * 256 + j;
			tile_area[tile_index] = 0;
			/* Calculate the web mercartor coordinate of the pixel */
			double	web_mercartor_x = (256 * x + j + 0.5) * Rz - L;
			double	web_mercartor_y = L - (256 * y + i - 0.5) * Rz;
			/* Query the spatial objects with InnerBox1(if there are lots of spatial objects within the distance R1 from the pixel, we query the spatial objects intersects the InnerBox1) */
			box					InnerBox1( point( web_mercartor_x - r, web_mercartor_y - r ), point( web_mercartor_x + r, web_mercartor_y + r ) );
			rtree_segment::const_query_iterator	it = rtree_segment_ptr->qbegin( bgi::intersects( InnerBox1 ) );
			if ( it != rtree_segment_ptr->qend() )
			{
				tile_area[tile_index] += 4;
			}else  {
				/* Query the spatial objects with OuterBox1(if there are few spatial objects in the neighbor of the pixel, we use the OuterBox1 to filter out the spatial objects which are far from the pixel) */
				box			OuterBox1( point( web_mercartor_x - R1, web_mercartor_y - R1 ), point( web_mercartor_x + R1, web_mercartor_y + R1 ) );
				std::vector<segment>	segment_result1;
				rtree_segment_ptr->query( bgi::intersects( OuterBox1 ) && bgi::nearest( point( web_mercartor_x, web_mercartor_y ), 1 ),
							  std::back_inserter( segment_result1 ) );
				if ( segment_result1.size() > 0 )
				{
					tile_area[tile_index] += 4;
				}else  {
					/* Determine whether the pixel belongs to the color transition regions (calculate the number of sub-pixels that are in the plotting region) */
					box			OuterBox2( point( web_mercartor_x - R2, web_mercartor_y - R2 ), point( web_mercartor_x + R2, web_mercartor_y + R2 ) );
					std::vector<segment>	segment_result2;
					rtree_segment_ptr->query( bgi::intersects( OuterBox2 ) && bgi::nearest( point( web_mercartor_x, web_mercartor_y ), 1 ),
								  std::back_inserter( segment_result2 ) );
					if ( segment_result2.size() > 0 )
					{
						if ( bg::distance( point( web_mercartor_x + Rz4, web_mercartor_y + Rz4 ), segment_result2.front() ) < R )
							tile_area[tile_index] += 1;
						if ( bg::distance( point( web_mercartor_x - Rz4, web_mercartor_y - Rz4 ), segment_result2.front() ) < R )
							tile_area[tile_index] += 1;
						if ( bg::distance( point( web_mercartor_x + Rz4, web_mercartor_y - Rz4 ), segment_result2.front() ) < R )
							tile_area[tile_index] += 1;
						if ( bg::distance( point( web_mercartor_x - Rz4, web_mercartor_y + Rz4 ), segment_result2.front() ) < R )
							tile_area[tile_index] += 1;
					}
				}
			}
		}
	}
}


/**
 * @brief Visualize Point Objects
 * @param z                     Zoom level of the tile
 * @param x                     x coordinate of the tile
 * @param y                     y coordinate of the tile
 * @param indexPath     Spatial index path
 * @param dataId                Dataset id
 * @param tile_area     Matrix indentifing different regions(e.g. background region, color transition region, the zones of rasterized spatial objects)
 *
 */
void  PointVision( int z, int x, int y, char* indexPath, char *dataId, char *tile_area )
{
	char		* tmp = new char[256];
	double		Rz, R, R1, R2, r, Rz4;
	rtree_point	* rtree_point_ptr;
	Rz	= L / (128 << z);       /* Calculate the pixel resolution of zoom level z */
	R	= 2 * Rz;               /* R=2*Rz ,we set the visualizaiton radius N to 2 (pixels) */
	R1	= 1.6464 * Rz;          /* R1=R-sqrt(2)/4*Rz=(2-sqrt(2)/4)*Rz */
	R2	= 2.3536 * Rz;          /* R2=R+sqrt(2)/4*Rz=(2+sqrt(2)/4)*Rz */
	r	= 0.7071 * R1;          /* r=sqrt(2)/2*R1 */
	Rz4	= 0.25 * Rz;            /* Rz4=1/4*Rz */
	sprintf( tmp, "%s%s", indexPath, dataId );
	bi::managed_mapped_file file( bi::open_only, tmp );
	rtree_point_ptr = file.find<rtree_point>( "rtree" ).first;
	#pragma omp parallel for num_threads(2)
	for ( int i = 0; i < 256; i++ )
	{
		for ( int j = 0; j < 256; j++ )
		{
			int tile_index = i * 256 + j;
			tile_area[tile_index] = 0;
			/* Calculate the web mercartor coordinate of the pixel */
			double	web_mercartor_x = (256 * x + j + 0.5) * Rz - L;
			double	web_mercartor_y = L - (256 * y + i - 0.5) * Rz;
			/* Query the spatial objects with InnerBox1(if there are lots of spatial objects within the distance R1 from the pixel, we query the spatial objects intersects the InnerBox1) */
			box					InnerBox1( point( web_mercartor_x - r, web_mercartor_y - r ), point( web_mercartor_x + r, web_mercartor_y + r ) );
			rtree_point::const_query_iterator	it = rtree_point_ptr->qbegin( bgi::intersects( InnerBox1 ) );
			if ( it != rtree_point_ptr->qend() )
			{
				tile_area[tile_index] += 4;
			}else  {
				/* Query the spatial objects with OuterBox1(if there are few spatial objects in the neighbor of the pixel, we use the OuterBox1 to filter out the spatial objects which are far from the pixel) */
				box			OuterBox1( point( web_mercartor_x - R1, web_mercartor_y - R1 ), point( web_mercartor_x + R1, web_mercartor_y + R1 ) );
				std::vector<point>	point_result1;
				rtree_point_ptr->query( bgi::intersects( OuterBox1 ) && bgi::nearest( point( web_mercartor_x, web_mercartor_y ), 1 ),
							std::back_inserter( point_result1 ) );
				if ( point_result1.size() > 0 )
				{
					tile_area[tile_index] += 4;
				}else  {
					/* Determine whether the pixel belongs to the color transition regions (calculate the number of sub-pixels that are in the plotting region) */
					box			OuterBox2( point( web_mercartor_x - R2, web_mercartor_y - R2 ), point( web_mercartor_x + R2, web_mercartor_y + R2 ) );
					std::vector<point>	point_result2;
					rtree_point_ptr->query( bgi::intersects( OuterBox2 ) && bgi::nearest( point( web_mercartor_x, web_mercartor_y ), 1 ),
								std::back_inserter( point_result2 ) );
					if ( point_result2.size() > 0 )
					{
						if ( bg::distance( point( web_mercartor_x + Rz4, web_mercartor_y + Rz4 ), point_result2.front() ) < R )
							tile_area[tile_index] += 1;
						if ( bg::distance( point( web_mercartor_x - Rz4, web_mercartor_y - Rz4 ), point_result2.front() ) < R )
							tile_area[tile_index] += 1;
						if ( bg::distance( point( web_mercartor_x + Rz4, web_mercartor_y - Rz4 ), point_result2.front() ) < R )
							tile_area[tile_index] += 1;
						if ( bg::distance( point( web_mercartor_x - Rz4, web_mercartor_y + Rz4 ), point_result2.front() ) < R )
							tile_area[tile_index] += 1;
					}
				}
			}
		}
	}
}


/**
 * @brief Visualize Pologon Objects
 * @param z                     Zoom level of the tile
 * @param x                     x coordinate of the tile
 * @param y                     y coordinate of the tile
 * @param indexPath     Spatial index path
 * @param dataId                Dataset id
 * @param tile_area     Matrix indentifing different regions(e.g. background region, color transition region, the zones of rasterized spatial objects)
 *
 */
void  PologonVision( int z, int x, int y, char* indexPath, char *dataId, char *tile_area )
{
	char		* tmp = new char[256];
	double		Rz, R, R1, R2, r, Rz4;
	rtree_psegment	* rtree_psegment_ptr;
	rtree_box	* rtree_box_ptr;
	Rz	= L / (128 << z);       /* Calculate the pixel resolution of zoom level z */
	R	= 2 * Rz;               /* R=2*Rz ,we set the visualizaiton radius N to 2 (pixels) */
	R1	= 1.6464 * Rz;          /* R1=R-sqrt(2)/4*Rz=(2-sqrt(2)/4)*Rz */
	R2	= 2.3536 * Rz;          /* R2=R+sqrt(2)/4*Rz=(2+sqrt(2)/4)*Rz */
	r	= 0.7071 * R1;          /* r=sqrt(2)/2*R1 */
	Rz4	= 0.25 * Rz;            /* Rz4=1/4*Rz */
	sprintf( tmp, "%s%s", indexPath, dataId );
	bi::managed_mapped_file file( bi::open_only, tmp );
	rtree_psegment_ptr = file.find<rtree_psegment>( "rtree" ).first;
	sprintf( tmp, "%s%s_mbr", indexPath, dataId );
	bi::managed_mapped_file file_mbr( bi::open_only, tmp );
	rtree_box_ptr = file_mbr.find<rtree_box>( "rtree" ).first;
	#pragma omp parallel for num_threads(2)
	for ( int i = 0; i < 256; i++ )
	{
		for ( int j = 0; j < 256; j++ )
		{
			int tile_index = i * 256 + j;
			tile_area[tile_index] = 0;
			/* Calculate the web mercartor coordinate of the pixel */
			double	web_mercartor_x = (256 * x + j + 0.5) * Rz - L;
			double	web_mercartor_y = L - (256 * y + i - 0.5) * Rz;
			/* Spatial-Index-Based Filling */
			std::vector<polygon_box> box_result;
			rtree_box_ptr->query( bgi::intersects( point( web_mercartor_x, web_mercartor_y ) ), std::back_inserter( box_result ) ); /* Find the candidate polygons */
			sort( box_result.begin(), box_result.end(), SortPolygon );                                                              /* Sort the polygons(Polygon with smaller x span has higher priority) */
			BOOST_FOREACH( polygon_box const & v, box_result )                                                                      /* measure the spatial relationship between the pixel and each candidate polygon one by one until the polygon which contains the pixel is found */
			{
				unsigned long			pID	= v.second;
				box				pBox	= v.first;
				double				minx	= bg::get<0, 0>( pBox );
				double				maxx	= bg::get<1, 0>( pBox );
				int				lcount	= 0;
				std::vector<polygon_segment>	psegment_result;
				if ( web_mercartor_x - minx < maxx - web_mercartor_x )
					rtree_psegment_ptr->query( bgi::intersects( segment( point( minx, web_mercartor_y ), point( web_mercartor_x, web_mercartor_y ) ) ),
								   std::back_inserter( psegment_result ) );
				else
					rtree_psegment_ptr->query( bgi::intersects( segment( point( web_mercartor_x, web_mercartor_y ), point( maxx, web_mercartor_y ) ) ),
								   std::back_inserter( psegment_result ) );
				BOOST_FOREACH( polygon_segment const & sv, psegment_result )
				{
					bool		level	= boost::get<2>( sv );
					unsigned long	pIDtmp	= boost::get<1>( sv );
					if ( level && pID == pIDtmp )
						lcount++;
				}
				if ( lcount % 2 > 0 ) /* pixel is classified as ’inside the polygon’ if the number of crossings is odd, or ’outside’ if it is an even number */
				{
					tile_area[tile_index] = 4;
					break;
				}
			}

			/* Query the spatial objects with InnerBox1(if there are lots of spatial objects within the distance R1 from the pixel, we query the spatial objects intersects the InnerBox1) */
			box					InnerBox1( point( web_mercartor_x - r, web_mercartor_y - r ), point( web_mercartor_x + r, web_mercartor_y + r ) );
			rtree_psegment::const_query_iterator	it = rtree_psegment_ptr->qbegin( bgi::intersects( InnerBox1 ) );
			if ( it != rtree_psegment_ptr->qend() )
			{
				if ( tile_area[tile_index] > 0 )
					tile_area[tile_index] += 4;
				else
					tile_area[tile_index] -= 4;
			}else  {
				/* Query the spatial objects with OuterBox1(if there are few spatial objects in the neighbor of the pixel, we use the OuterBox1 to filter out the spatial objects which are far from the pixel) */
				box				OuterBox1( point( web_mercartor_x - R1, web_mercartor_y - R1 ), point( web_mercartor_x + R1, web_mercartor_y + R1 ) );
				std::vector<polygon_segment>	segment_result1;
				rtree_psegment_ptr->query( bgi::intersects( OuterBox1 ) && bgi::nearest( point( web_mercartor_x, web_mercartor_y ), 1 ),
							   std::back_inserter( segment_result1 ) );
				if ( segment_result1.size() > 0 )
				{
					if ( tile_area[tile_index] > 0 )
						tile_area[tile_index] += 4;
					else
						tile_area[tile_index] -= 4;
				}else  {
					/* Determine whether the pixel belongs to the color transition regions (calculate the number of sub-pixels that are in the plotting region) */
					box				OuterBox2( point( web_mercartor_x - R2, web_mercartor_y - R2 ), point( web_mercartor_x + R2, web_mercartor_y + R2 ) );
					std::vector<polygon_segment>	segment_result2;
					rtree_psegment_ptr->query( bgi::intersects( OuterBox2 ) && bgi::nearest( point( web_mercartor_x, web_mercartor_y ), 1 ),
								   std::back_inserter( segment_result2 ) );
					if ( segment_result2.size() > 0 )
					{
						segment pSegment = boost::get<0>( segment_result2.front() );
						if ( tile_area[tile_index] > 0 )
						{
							if ( bg::distance( point( web_mercartor_x + Rz4, web_mercartor_y + Rz4 ), pSegment ) < R )
								tile_area[tile_index] += 1;
							if ( bg::distance( point( web_mercartor_x - Rz4, web_mercartor_y - Rz4 ), pSegment ) < R )
								tile_area[tile_index] += 1;
							if ( bg::distance( point( web_mercartor_x + Rz4, web_mercartor_y - Rz4 ), pSegment ) < R )
								tile_area[tile_index] += 1;
							if ( bg::distance( point( web_mercartor_x - Rz4, web_mercartor_y + Rz4 ), pSegment ) < R )
								tile_area[tile_index] += 1;
						}else{
							if ( bg::distance( point( web_mercartor_x + Rz4, web_mercartor_y + Rz4 ), pSegment ) < R )
								tile_area[tile_index] -= 1;
							if ( bg::distance( point( web_mercartor_x - Rz4, web_mercartor_y - Rz4 ), pSegment ) < R )
								tile_area[tile_index] -= 1;
							if ( bg::distance( point( web_mercartor_x + Rz4, web_mercartor_y - Rz4 ), pSegment ) < R )
								tile_area[tile_index] -= 1;
							if ( bg::distance( point( web_mercartor_x - Rz4, web_mercartor_y + Rz4 ), pSegment ) < R )
								tile_area[tile_index] -= 1;
						}
					}
				}
			}
		}
	}
}


/**
 * @brief HiVision visualization engine main function
 * @param nArgc         Input parameters copunt
 * @param papszArgv Input startup parameters
 *
 */
int main( int nArgc, char ** papszArgv )
{
	int myId, numProcs;
	MPI_Init( &nArgc, &papszArgv );
	MPI_Comm_rank( MPI_COMM_WORLD, &myId );
	MPI_Comm_size( MPI_COMM_WORLD, &numProcs );
	double	t1, t2;
	char	* indexPath	= papszArgv[1];         /* Path to store the spatial indexes */
	char	* redisHost	= papszArgv[2];         /* Host IP of Redis */
	int	redisPort	= atoi( papszArgv[3] ); /* Host Port of Redis */
	CPLSetConfigOption( "GDAL_FILENAME_IS_UTF8", "NO" );
	CPLSetConfigOption( "SHAPE_ENCODING", "UTF-8" );
	/* Connect Redis */
	Redis *redis = new Redis();
	if ( !redis->connect( redisHost, redisPort ) )
	{
		printf( "connect redis error!\n" );
		exit( 0 );
		MPI_Finalize();
	}
	char	*tile_area	= (char *) malloc( TILE_SIZE * TILE_SIZE );
	char	*buffer_tmp	= (char *) malloc( TILE_SIZE * TILE_SIZE );
	char	* task		= new char[256];
	char	* tile_params[MAX_TILE_PARAMS];
	int	count = 0;
	int	x, y, z;
	if ( myId == 0 )
	{
		printf( "indexPath:%s\n", indexPath );
		printf( "Service Start. cores:%d\n", numProcs );
	}
	while ( 1 )
	{
		/* Get tasks from Task Pool */
		sprintf( task, "%s", redis->brpop( "HiVisiontasklist" ).c_str() );
		try{
			if ( strlen( task ) > 0 )
			{
				t1 = MPI_Wtime();
				/* Parse the parameters */
				GetList( task, tile_params, (char *) "/", count );
				z	= atoi( tile_params[1] );
				x	= atoi( tile_params[2] );
				y	= atoi( tile_params[3] );
				if ( tile_params[0][0] == 'l' )                         /* Linestring visualizaiton */
					LineVision( z, x, y, indexPath, tile_params[0], tile_area );
				else if ( tile_params[0][0] == 'p' )                    /* Point visualizaiton */
					PointVision( z, x, y, indexPath, tile_params[0], tile_area );
				else if ( tile_params[0][0] == 'a' )                    /* Polygon visulization */
					PologonVision( z, x, y, indexPath, tile_params[0], tile_area );
				redis->zset( task, tile_area, TILE_SIZE * TILE_SIZE );  /* Write the result to the Result Pool */
				while ( !redis->zget( task, buffer_tmp ) )
					redis->zset( task, tile_area, TILE_SIZE * TILE_SIZE );
				redis->expire( task, "1000" );                          /* Set expire time to 1000s(expired results will be cleaned up if memory usage reaches the upper limit) */
				redis->pub( "HiVisiontiles", task );                    /* Send task completion message */
				t2 = MPI_Wtime();
				printf( "tile-%s-%d-%f\n", task, myId, t2 - t1 );
			}
		}
		catch ( ... )
		{
			printf( "Error task %s \n", task );
		}
	}
	MPI_Finalize();
}


