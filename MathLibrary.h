// MathLibrary.h - Contains declarations of math functions
// Author: Daniel
#pragma once

#ifdef MATHLIBRARY_EXPORTS
#define MATHLIBRARY_API __declspec(dllexport)
#else
#define MATHLIBRARY_API __declspec(dllimport)
#endif

MATHLIBRARY_API typedef struct Point {
	int x, y;
} Point;

MATHLIBRARY_API typedef struct Rect {
	int x, y, width, height;  // x, y is topleft corner of rectangle
} Rect;

MATHLIBRARY_API typedef struct Point_And_ID {
	Point point;	// this point characterizes the position of the agent
	int id;			// this id is a unique number that references a rect/agent in the client
} Point_And_ID;

MATHLIBRARY_API typedef struct TileCoord
{
	int x;
	int y;
} TileCoord;

typedef struct Graph2D Graph2D;

typedef int bool;
#define false 0
#define true 1

// Quadtree interface. Used for fast collision detection.
extern MATHLIBRARY_API void quadtree_insert_all_and_delete_previous_quadtree_api(Rect r, Point_And_ID* points, int num_points);
extern MATHLIBRARY_API void quadtree_query_api(Rect r);
extern MATHLIBRARY_API int get_point_count_api();
extern MATHLIBRARY_API Point_And_ID* get_data_pointer();
extern MATHLIBRARY_API int* get_apple_pointer();

// Dijkstra algorithm interface. Used for pathfinding.
extern MATHLIBRARY_API Graph2D* DijkstraCreateGraphAPI(int width, int height, int* tile_data);
extern MATHLIBRARY_API int DijkstraFindPathsAPI(Graph2D* graph, int x_end, int y_end);
extern MATHLIBRARY_API void DijkstraUpdateGraphAPI(Graph2D* graph, unsigned short* xy, double* cost_adjustment, int count);
extern MATHLIBRARY_API TileCoord DijkstraPathAPI(Graph2D* graph, int x_start, int y_start, int x_goal, int y_goal);
