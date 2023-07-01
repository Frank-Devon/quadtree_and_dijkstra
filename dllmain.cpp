#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <malloc.h>
#include <memory.h>
#include <math.h>
#include <time.h>
#include "MathLibrary.h"

#define MAX_POINTS 4

bool point_compare(Point p1, Point p2)
{
	return p1.x == p2.x && p1.y == p2.y;
}

bool point_in_rect(Rect* r, Point p)
{
	bool inside_vert = r->x <= p.x && r->x + r->width > p.x;
	bool inside_horz = r->y <= p.y && r->y + r->height > p.y;
	return inside_vert && inside_horz;
}

bool collide_rect(Rect* r1, Rect* r2)
{
	//r1 = 1,1,2,2
	//r2 = 3,1,2,2
	bool b1 = r1->x + r1->width < r2->x;  //right edge > other left edge
	bool b2 = r1->x > r2->x + r2->width;
	bool b3 = r1->y + r1->height < r2->y;
	bool b4 = r1->y > r2->y + r2->height;
	return !(b1 || b2 || b3 || b4);
}

typedef struct quadtree quadtree;

struct quadtree { // This structure is named "myDataType"
	Point_And_ID* data_ptr[MAX_POINTS];  //point points[MAX_POINTS];
	int iPointCount;
	Rect boundary;
	quadtree* subTrees[4];
	bool bIsDivided;
	int iMaxPoints;
	int iDepth;
};

//typedef struct Point_And_ID Point_And_ID;

// globals 
Point_And_ID g_data[1000]; // large enough for virtually all cases
int g_data_count = 0; // how many point_and_id structs are in use in the statically allocated array.
quadtree qt_main;
int num_points = 0;  // amount of points successfully queried
Point_And_ID points_queried[1000];  // CAREFUL, the intended queries are supposed to yield small amount of points... too much and array will overflow

void quadtree_init(quadtree* qt, Rect boundary, int depth)
{
	qt->iPointCount = 0;
	qt->boundary = boundary;
	qt->subTrees[0] = NULL;
	qt->subTrees[1] = NULL;
	qt->subTrees[2] = NULL;
	qt->subTrees[3] = NULL;
	qt->bIsDivided = false;
	qt->iMaxPoints = MAX_POINTS;
	qt->iDepth = depth;
}

bool quadtree_insert(quadtree* qt, Point_And_ID p_id)
{
	if (!point_in_rect(&(qt->boundary), p_id.point))
	{
		return false;  // if point not in boundary, leave
	}
	if (qt->iPointCount < MAX_POINTS)
	{
		// add point to this quadtree node
		g_data[g_data_count] = p_id;
		qt->data_ptr[qt->iPointCount] = &g_data[g_data_count];
		qt->iPointCount = qt->iPointCount + 1;
		g_data_count += 1;
		return true;
	}
	// can't fit point in this node.
	if (qt->bIsDivided == false)
	{
		qt->bIsDivided = true;
		quadtree* m = (quadtree*)malloc(4 * sizeof(struct quadtree));
		qt->subTrees[0] = m;
		qt->subTrees[1] = (m + 1);
		qt->subTrees[2] = (m + 2);
		qt->subTrees[3] = (m + 3);
		// todo rounding issues? especially when height and width are odd numbers?
		Rect r[4];
		int new_width = qt->boundary.width / 2;
		int new_height = qt->boundary.height / 2;
		//NE rect
		r[0].x = qt->boundary.x + new_width;
		r[0].y = qt->boundary.y;
		r[0].width = new_width;
		r[0].height = new_height;
		//NW rect
		r[1].x = qt->boundary.x;
		r[1].y = qt->boundary.y;
		r[1].width = new_width;
		r[1].height = new_height;
		//SW rect
		r[2].x = qt->boundary.x;
		r[2].y = qt->boundary.y + new_height;
		r[2].width = new_width;
		r[2].height = new_height;
		//SE rect
		r[3].x = qt->boundary.x + new_width;
		r[3].y = qt->boundary.y + new_height;
		r[3].width = new_width;
		r[3].height = new_height;

		for (int i = 0; i < 4; i++)
		{
			quadtree_init(qt->subTrees[i], r[i], qt->iDepth + 1);
		}
	}

	return (quadtree_insert(qt->subTrees[0], p_id) ||
		quadtree_insert(qt->subTrees[1], p_id) ||
		quadtree_insert(qt->subTrees[2], p_id) ||
		quadtree_insert(qt->subTrees[3], p_id));
}

//todo fix return  type
void quadtree_query(quadtree* qt, Rect* r)
{
	if (!collide_rect(&(qt->boundary), r))
	{
		return;
	}
	for (int i = 0; i < qt->iPointCount; i++)
	{
		if (point_in_rect(r, qt->data_ptr[i]->point))
		{
			points_queried[num_points] = *(qt->data_ptr[i]);
			num_points++;
		}
	}
	if (qt->bIsDivided)
	{
		for (int i = 0; i < 4; i++)
		{
			quadtree_query(qt->subTrees[i], r);
		}
	}
}

void quadtree_free(quadtree* qt)
{
	if (!qt->bIsDivided)
	{
		return;
	}
	for (int i = 0; i < 4; i++)
	{
		quadtree_free(qt->subTrees[i]);
	}
	for (int i = 0; i < 4; i++)
	{
		if (qt->subTrees[i]->bIsDivided)
		{
			free(qt->subTrees[i]);
		}
	}
	qt->bIsDivided = false;
	return;
}

void quadtree_insert_all_and_delete_previous_quadtree_api(Rect r, Point_And_ID* points_and_ids, int num_points)
{
	//call delete on qt_main... also zero out data points iPointsCount = 0;
	g_data_count = 0;
	quadtree_free(&qt_main);
	qt_main.iPointCount = 0;
	quadtree_init(&qt_main, r, 0);
	for (int i = 0; i < num_points; i++)
	{
		quadtree_insert(&qt_main, points_and_ids[i]);
		//g_data_count++;
	}
}

void quadtree_query_api(Rect r)
{
	num_points = 0;
	quadtree_query(&qt_main, &r);
}

int get_point_count_api()
{
	return num_points;
}

Point_And_ID* get_data_pointer()
{
	return &points_queried[0];
}


////
//// PART 2
//// pathfinding and other helping functions
////

enum TileSet { None = 0, Unvisited = 1, Visited = 2 };
typedef enum TileSet TileSet;
typedef struct Tile Tile;
typedef struct Edge Edge;
typedef struct Heap Heap;

typedef struct Edge
{
	Tile* vertex;  // if this vertex is null edge is non existant
	double cost;  //
};

typedef struct Tile
{
	// cost to enter tile from adjacent tile without considering other factors. the actual cost
	// from one vertex to another is stored in the edges field struct.
	double cost_base;
	//double cost_total;  // cost to enter tile when considering other factors (agents moving/blocking in tile)
	// considers other factors when determining cost to this node. Stopped/slowed agents make this number bigger
	// or infinite.  NORMALLY THIS SHOULD BE 1.0, never 0.0
	double cost_modifier;
	Tile* parent;
	unsigned char passable;  // // get rid of? use cost_base = INFINITY instead? 0 = impassable, 0 != passable
	TileSet set;  // none, unvisited, visited
	double g_cost, h_cost; // g_cost = tentative cost? h_cost not needed?
	unsigned short x_index, y_index; // might not be needed?
	Edge edges[8];
	//incoming edges... Edge* edges_incoming[8];
};

typedef struct Graph2D
{
	Tile* tiles;  // pointer to 2d array of tiles
	Tile** tiles_dirty;  // pointer to tiles that need to be refreshed
	//Heap* heap;  // increases speed when finding minimum tentative cost in unvisited list
	int width;
	int height;
	int tiles_dirty_count;
};

typedef struct Heap
{
	void* (*compare_nodes)(void*, void*);  // returns void* to highest priority node
	int (*compare_nodes_index)(Heap*, int, int);  // gets index of highest priority node
	void** array;  // must be dynamically allocated  
	int size;  // max size of array
	int count;  // currently used elements/nodes in the array
	//int element_size;  // probably not needed !!!!!!!!!!!!!
};

// if tiles have same g_cost, return the tile_0. Could be improved.
void* compare_tiles_get_min_g_cost(void* tile_0, void* tile_1)
{
	if (tile_0 != NULL && tile_1 != NULL)
		return (((Tile*)tile_0)->g_cost <= ((Tile*)tile_1)->g_cost) ? tile_0 : tile_1;
	else if (tile_0 == NULL && tile_1 != NULL)
		return tile_1;
	else if (tile_0 != NULL && tile_1 == NULL)
		return tile_0;
	else
		return NULL;
}

int compare_tiles_get_min_g_cost_index(Heap* heap, int index_0, int index_1)
{
	// both indexes must be valid
	// no error checking needed?
	//return ((Tile*)*(heap->array + index_0))->g_cost <= ((Tile*)*(heap->array + index_1))->g_cost ?
	//	index_0 : index_1;

	if (index_0 != -1 && index_1 != -1)
		return ((Tile*)*(heap->array + index_0))->g_cost <= ((Tile*)*(heap->array + index_1))->g_cost ?
		index_0 : index_1;
	else if (index_0 == -1 && index_1 != -1)
		return index_1;
	else if (index_0 != -1 && index_1 == -1)
		return index_0;
	else
		return -1;
}

// creates a heap data structure. After calling this, actual data must be inserted into data structure.
Heap* HeapCreate(int heap_size, void* (*compare_nodes)(void*, void*), int (*compare_nodes_index)(Heap*, int, int))
{
	Heap* heap = (Heap*)malloc(sizeof(struct Heap));
	if (heap == NULL)
	{
		return NULL;
	}
	//heap->array = (char*)malloc((unsigned long long)heap_size * node_size);
	heap->array = (void**)malloc((unsigned long long)heap_size * sizeof(void*));
	if (heap->array == NULL)
	{
		free(heap);
		return NULL;
	}
	heap->count = 0;
	heap->size = heap_size;
	heap->compare_nodes = compare_nodes;  // determines if min or max heap
	heap->compare_nodes_index = compare_nodes_index;
	//heap->element_size = node_size;  // probably not needed
	return heap;
}

int HeapFree1(Heap* heap)
{
	free(heap->array);
}

// returns 0 for false, non zero for True. 
int HeapIndexInBounds(Heap* heap, int index)
{
	if ( index < 0 || index > heap->count)
	{
		return 0;  // False
	}
	else
	{
		return 1;  // True
	}
}

void* HeapNodeAtIndex(Heap* heap, int index)
{
	if (HeapIndexInBounds(heap, index))
	{
		//return (char*)(heap->array + index * heap->element_size);
		return *(heap->array + index);
	}
	else
	{
		return NULL;
	}
}

// gets parent node at index
void* HeapParent(Heap* heap, int index)
{
	if ( index == 0 || index > heap->count)
	{
		return NULL;
	}
	else
	{
		int index_new = (index - 1) / 2;
		return *(heap->array + index_new);
	}
}

int HeapParentIndex(Heap* heap, int index)
{
	// too much error checking?
	if ( index >= heap->count || index <= 0)
	{
		return -1;  //invalid
	}
	else
	{
		return (index - 1) / 2;
	}
}

// gets left child node
void* HeapLeftChild(Heap* heap, int index) // returns NULL on out of bounds index
{
	if ( index >= heap->count || (2 * index + 1 >= heap->count))
	{
		return NULL;
	}
	else
	{
		return *(heap->array + (2 * index + 1));
	}
}

int HeapLeftChildIndex(Heap* heap, int index)
{
	// too much error checking?
	if ( index >= heap->count || (2 * index + 1 >= heap->count))
	{
		return -1;  //invalid
	}
	else
	{
		return 2 * index + 1;
	}
}

// gets right child node
void* HeapRightChild(Heap* heap, int index) // returns NULL on out of bounds index
{
	if ( index >= heap->count || (2 * index + 2 >= heap->count))
	{
		return NULL;
	}
	else
	{
		return *(heap->array + (2 * index + 2));
	}
}

int HeapRightChildIndex(Heap* heap, int index)
{
	// too much error checking?
	if ( index >= heap->count || (2 * index + 2 >= heap->count))
	{
		return -1;  //invalid
	}
	else
	{
		return 2 * index + 2;
	}
}

int HeapifyUp(Heap* heap, int index)
{
	void* parent = HeapParent(heap, index);
	void* node = HeapNodeAtIndex(heap, index);
	int parent_index = HeapParentIndex(heap, index);

	if (parent == NULL || node == NULL)
	{
		//printf("Heap: parent or node invalid\n");
		return 1;
	}
	void* result = heap->compare_nodes(parent, node); // compare nodes returns higher priority?
	if (result == parent)
	{
		// heap property is still valid
		return 1;
	}
	else
	{

		HeapSwapNodes(heap, parent_index, index);
		HeapifyUp(heap, (index - 1) / 2);
		return 1;
	}
}

int HeapifyDown(Heap* heap, int index)
{
	void* node_current = HeapNodeAtIndex(heap, index);
	void* node_right_child = HeapRightChild(heap, index);
	void* node_left_child = HeapLeftChild(heap, index);
	int right_child_index = HeapRightChildIndex(heap, index);
	int left_child_index = HeapLeftChildIndex(heap, index);

	if (heap == NULL)
	{
		return 0; //heap ptr is bad
	}
	//// get highest priority child
	//void* node_child_highest_priority = heap->compare_nodes(node_left_child, node_right_child);
	//int index_child_highest_priority = HeapIndexOfNode(heap, node_child_highest_priority);  // don't use this func?
	//void* result = heap->compare_nodes(node_current, node_child_highest_priority);

	int index_child_highest_priority = heap->compare_nodes_index(heap, right_child_index, left_child_index);
	// check for invalid -1
	//int index_highest_priority = heap->compare_nodes_index(heap, index, index_child_highest_priority);
	if (index_child_highest_priority == -1)
	{
		// heap is valid (at the very end of heap)
	}
	else
	{
		int index_highest_priority = heap->compare_nodes_index(heap, index, index_child_highest_priority);
		if (index == index_highest_priority)
		{
			// heap is valid
		}
		else
		{
			HeapSwapNodes(heap, index, index_highest_priority);
			HeapifyDown(heap, index_highest_priority);
		}
	}
}

int HeapSwapNodes(Heap* heap, int index_0, int index_1)
{
	void* temp = *(heap->array + index_0);
	*(heap->array + index_0) = *(heap->array + index_1);
	*(heap->array + index_1) = temp;
}

int HeapInsert(Heap* heap, void* data)
{
	// copy from data to internal heap
	if (!(heap->size > heap->count))
	{
		// internal array not big enough. should never/rarely happen. exit and post message.
		printf("Heap: array to small. Cannot add element to array.\n");
		return 0;
	}

	//char* destination = heap->array + heap->count * heap->element_size;
	//memcpy(destination, data, heap->element_size);
	*(heap->array + heap->count) = data;  // is this right?
	heap->count = heap->count + 1;
	HeapifyUp(heap, heap->count - 1);
	return 1;
}

void* HeapPeek(Heap* heap)
{
	if (heap == NULL || heap->count <= 0)
	{
		return NULL;
	}
	return *(heap->array);
}

// bad code... internal array should just hold pointers to actual data,
// instead of actually containing the data itself
void* HeapPoll(Heap* heap)
{
	if (heap == NULL || heap->count <= 0)
	{
		printf("HeapPoll func: heap is not initialized or count = 0.\n");
		return NULL;
	}

	void* result = HeapNodeAtIndex(heap, 0);
	//HeapSwapNodes(heap, *(heap->array), *(heap->array + heap->count - 1));
	HeapSwapNodes(heap, 0, heap->count - 1);
	heap->count -= 1;
	HeapifyDown(heap, 0);
	return result;
}

//
// DIJKSTRA ALGORITHM
//

//Graph2D* DijkstraCreateGraphAPI(int width, int height, int* tile_data);
//int DijkstraFindPathsAPI(Graph2D* graph, int x_end, int y_end);
//void DijkstraUpdateGraphAPI(Graph2D* graph, unsigned short* xy, double* cost_adjustment, int count);
//TileCoord DijkstraPathAPI(Graph2D* graph, int x_start, int y_start, int x_goal, int y_goal);

void TestDijkstraUpdateGraph(Graph2D* graph);
void DijkstraFindEdges(Graph2D* graph, int x, int y);
Graph2D* DijkstraCreateGraph(int width, int height, Tile* tiles_arg);
Tile* GraphGetVertex(Graph2D* graph, int x, int y);
Tile* TestTiles();


Tile* GraphGetVertex(Graph2D* graph, int x, int y)
{
	return (graph->tiles + y * graph->width + x);
}

// this should be used internally only
Graph2D* DijkstraCreateGraph(int width, int height, Tile* tiles_arg)
{
	Graph2D* graph;
	graph = (Graph2D*)malloc(sizeof(Graph2D)); // probably should check for null
	graph->tiles = tiles_arg;
	// does heap need to be bigger than the maximum number of vertexs? no
	graph->width = width;
	graph->height = height;
	graph->tiles_dirty = NULL;
	graph->tiles_dirty_count = 0;

	// Get graph edges
	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			DijkstraFindEdges(graph, x, y);
		}
	}
	return graph;
}


// intended to be called by consumers of this dll
Graph2D* DijkstraCreateGraphAPI(int width, int height, int* tile_data)
{
	// tile_data is a 2d array dimensions width and height, value 1 means passable, 0 impassable
	// this is an easier datatype for the python code to keep track of
	Tile* tiles;
	tiles = (Tile*)malloc(sizeof(Tile) * width * height);
	if (tiles == NULL)
	{
		printf("couldn't allocate memory for tiles\n");
		return NULL;
	}
	for (int i = 0; i < width * height; i++)
	{
		(*(tiles + i)).cost_base = (*(tile_data + i) == 1) ? 1.0 : INFINITY;  // 1 means passable
		(*(tiles + i)).cost_modifier = 1.0;
		(*(tiles + i)).passable = *(tile_data + i);
		(*(tiles + i)).x_index = i % width;
		(*(tiles + i)).y_index = i / width;
		(*(tiles + i)).set = None;
		(*(tiles + i)).parent = NULL;
		(*(tiles + i)).g_cost = INFINITY;
		(*(tiles + i)).h_cost = INFINITY;  // don't use h_cost?
	}
	// debug, check if data correct
	for (int i = 0; i < width*height; i++)
	{
		if (i % width == 0)
		{
			printf("\n");
		}
		printf("%d", (*(tiles + i)).passable);
	}
	return DijkstraCreateGraph(width, height, tiles);
}


// After Graph Vertexes (Tile map structure labeled Graph2D) have been mostly initialized, the edge array of 
// each vertex must be initialized. This should only be called once when setting up the graph edge list.
// To update the graph (due to temporary changes like agents obstructing paths), call DijkstraUpdateGraph 
// (which should be called every frame)
void DijkstraFindEdges(Graph2D* graph, int x, int y)
{
	Tile* tile_start = GraphGetVertex(graph, x, y);
	if (tile_start->passable == 0 || tile_start->cost_base == INFINITY)
	{
		return;
	}
	Tile* tile_end = NULL;
	Tile* tile_diagonal = NULL;  // used as a temp variable to check for tiles that may hinder diagonal movement
	int ydiff;
	int xdiff;
	double modifier;
	// origin is top left
	int x_offsets[8] = { 0, 1, 1, 1, 0, -1, -1, -1 };
	int y_offsets[8] = { -1, -1, 0, 1, 1, 1, 0, -1 };
	// check for boundary conditions
	int edge_count = 0;
	int reachable;  // 0 = false, 1 = true
	for (int i = 0; i < 8; i++)  // start i=0 Up, then clockwise
	{
		reachable = 0;  // false
		// check boundary conditions
		if (x + x_offsets[i] < 0 || x + x_offsets[i] >= graph->width)
		{
			continue;
		}
		if (y + y_offsets[i] < 0 || y + y_offsets[i] >= graph->height)
		{
			continue;
		}
		tile_end = GraphGetVertex(graph, x + x_offsets[i], y + y_offsets[i]);
		if (tile_end->passable == 0 || tile_end->cost_base == INFINITY)
		{
			continue;
		}
		if (abs(y_offsets[i]) + abs(x_offsets[i]) > 1) // diagonal edge
		{
			reachable = 0;  // False
			modifier = sqrt(2.0);
			// don't consider agents, just base terrain
			tile_diagonal = GraphGetVertex(graph, x + x_offsets[i], y);
			if (tile_diagonal->passable == 1 || tile_diagonal->cost_base != INFINITY)
			{
				reachable = 1;  // True
			}
			tile_diagonal = GraphGetVertex(graph, x, y + y_offsets[i]);
			if (tile_diagonal->passable == 1 || tile_diagonal->cost_base != INFINITY)
			{
				reachable = 1;  // True
			}
		}
		else // non diagonal edge
		{
			reachable = 1;  // True
			modifier = 1.0;
		}
		if (reachable == 1)
		{
			tile_start->edges[edge_count].vertex = tile_end;  //GraphGetVertex(graph, x + x_offsets[i], y + y_offsets[i]);
			// error accumulates when using floats?
			tile_start->edges[edge_count].cost = modifier * tile_start->edges[edge_count].vertex->cost_base;
			edge_count++;
		}
	}
	for (int i = edge_count; i < 8; i++)
	{
		tile_start->edges[i].vertex = NULL;
		tile_start->edges[i].cost = INFINITY;
	}
}


// updates the Dijkstra graph. Makes temporary changes to the graph. 
// Agents stoping on tiles can make them impassable. This function accounts for this.
// Any tiles updated this way must be marked dirty (and updated to normal values) for the next frame of the game.
void DijkstraUpdateGraphAPI(Graph2D* graph, unsigned short* xy, double* cost_adjustment, int count)
{
	//// TODO FIX!!!!! EDGES SHOULD BE DIRTY, NOT VERTEXS?
	if (graph == NULL) // || count == 0 || xy == NULL || cost_adjustment == NULL)
	{
		return;  // no work to be done.
	}
	Tile* tile;
	tile = NULL;
	for (int i = 0; i < graph->tiles_dirty_count; i++)  // clean up old dirty tiles
	{
		tile = *(graph->tiles_dirty + i);
		tile->cost_modifier = 1.0;
	}
	if (graph->tiles_dirty != NULL)
	{
		free(graph->tiles_dirty);
	}
	graph->tiles_dirty = (Tile**)malloc(count * sizeof(Tile*));
	graph->tiles_dirty_count = count;
	for (int i = 0; i < count; i++)
	{
		tile = GraphGetVertex(graph, *(xy + 2 * i), *(xy + 2 * i + 1));
		//tile->cost_total = tile->cost_base * *(cost_adjustment + i);
		tile->cost_modifier += *(cost_adjustment + i);
		*((graph->tiles_dirty) + i) = tile;  // mark this tile as dirty (reset cost_total)
	}

	// TODO reset tentative costs and visited status? etc
}


// Important. Finds shortest path from all vertexs (or tiles) to the end vertex (or tile).
// Pathfinding data saved in Tile data structure (from graph->tiles) 
// Tile* unvisited_array[10000];  // dynamically allocate later?
// returns 0 on failure, 1 on success
int DijkstraFindPathsAPI(Graph2D* graph, int x_end, int y_end)
{
	Tile* current;
	Heap* heap = HeapCreate(graph->width * graph->height * 5, compare_tiles_get_min_g_cost, compare_tiles_get_min_g_cost_index);
	// the "end" coordinates actually serve as the starting point of the algorithm
	// unvisited set
	// visited set
	current = GraphGetVertex(graph, x_end, y_end);
	if (current->passable == 0)
	{
		//printf("cannot path to impassable destination\n");
		return 0;
	}

	// make sure g_costs and parent values are cleared
	Tile* tile_clear;
	for (int y = 0; y < graph->height; y++)
	{
		for (int x = 0; x < graph->width; x++)
		{
			tile_clear = GraphGetVertex(graph, x, y);
			tile_clear->g_cost = INFINITY;
			tile_clear->parent = NULL;
			tile_clear->set = None;
		}
	}

	// Now that any previous pathfinding work has been erased, begin pathfinding algorithm.
	current->g_cost = 0;
	HeapInsert(heap, current); // insert initial tile into heap
	double tentative_cost = INFINITY;
	double new_cost;
	//int insert_count = 0;
	while (HeapPeek(heap) != NULL)
	{
		current = (Tile*)HeapPoll(heap);
		//printf("heap count %d; visiting vertex %d, %d; set %d\n", heap->count, current->x_index, current->y_index, current->set);
		if (current->set == Visited)
		{
			continue;
		}
		for (int i = 0; i < 8; i++)  // check neighbors
		{
			if (current->edges[i].vertex == NULL)
			{
				break;  // NULL means end of edge list
			}
			if (current->edges[i].vertex->set != Visited) 
			{
				new_cost = current->g_cost + current->edges[i].cost * current->edges[i].vertex->cost_modifier;
				if (new_cost < current->edges[i].vertex->g_cost)  // found better path to tile
				{
					// update g_cost and parent of neighboring tile (and also parent indexes)
					current->edges[i].vertex->g_cost = new_cost;
					current->edges[i].vertex->parent = current;
				}
				//printf("  %d, %d added to heap. count = %d\n", current->edges[i].vertex->x_index, current->edges[i].vertex->y_index, heap->count);
				current->edges[i].vertex->set = Unvisited;  // not needed?
				HeapInsert(heap, current->edges[i].vertex);
				//insert_count += 1;
			}
		}
		current->set = Visited;
	}
	HeapFree1(heap);

	//printf("insertions = %d\n", insert_count);
	// need to free heap?
	return 1;  // success
}

int DijkstraFindPathsAPI2(Graph2D* graph, int x_end, int y_end)
{
	clock_t start, end;
	start = clock();
	for (int i = 0; i < 200; i++)
	{
		DijkstraFindPathsAPI2(graph, x_end, y_end);
	}
	end = clock();
	printf("time Find Path %f\n", (double)(end - start) / CLOCKS_PER_SEC);

}


void TestDijkstraUpdateGraph(Graph2D* graph)
{
	unsigned short xy[4] = { 2, 2, 3, 5 }; // points 2,2 and 3,4
	double cost_adjustment[2] = { 2.0, 3.0 };
	int count = 2;
	DijkstraUpdateGraphAPI(graph, xy, cost_adjustment, count);
}

// creates a list of Tiles. Used to test DijkstraCreateGraph function
Tile* TestTiles()
{
	Tile* tiles;
	int width = 9; int height = 8;  // todo, programmatically get width and height of map
	tiles = (Tile*)malloc(sizeof(Tile) * width * height); // height is 8, width = 9
	// read from text file?
	FILE* fp; // declaration of file pointer
	char con[10000]; // variable to read the content
	fp = fopen("test_level.txt", "r");// opening of file
	if (!fp)  // checking for error
	{
		printf("cannot open test_level.txt\n");
		return NULL;
	}

	int x = 0;
	int y = 0;
	int count = 0;
	while (fgets(con, 10000, fp) != NULL)// reading file content
	{
		printf("%s", con);
		x = 0;
		while (con[x] == '1' || con[x] == '0')
		{
			if (con[x] == '1')
			{
				(*(tiles + y * width + x)).cost_base = INFINITY;
				(*(tiles + y * width + x)).cost_modifier = 1.0;
				(*(tiles + y * width + x)).passable = 0;
				(*(tiles + y * width + x)).x_index = x;
				(*(tiles + y * width + x)).y_index = y;
				(*(tiles + y * width + x)).set = None;
				(*(tiles + y * width + x)).parent = NULL;
				(*(tiles + y * width + x)).g_cost = INFINITY;
				(*(tiles + y * width + x)).h_cost = INFINITY;  // don't use h_cost?
				count++;
			}
			else if (con[x] == '0')
			{
				(*(tiles + y * width + x)).cost_base = 1.0;
				(*(tiles + y * width + x)).cost_modifier = 1.0;
				(*(tiles + y * width + x)).passable = 1;
				(*(tiles + y * width + x)).x_index = x;
				(*(tiles + y * width + x)).y_index = y;
				(*(tiles + y * width + x)).set = None;
				(*(tiles + y * width + x)).parent = NULL;
				(*(tiles + y * width + x)).g_cost = INFINITY;
				(*(tiles + y * width + x)).h_cost = INFINITY;  // don't use h_cost?
				count++;
			}
			x++;
		}
		y++;
	}
	fclose(fp); // closing file
	return tiles;
}

void* smallest_int(void* int1, void* int2)
{
	if (int1 != NULL && int2 != NULL)
		return (*(int*)int1 <= *(int*)int2) ? int1 : int2;
	else if (int1 != NULL && int2 != NULL)
		return int2;
	else if (int1 != NULL && int2 != NULL)
		return int1;
	else
		return NULL;
}

int smallest_int_index(Heap* heap, int index0, int index1)
{
	if (index0 != -1 && index1 != -1)
		return (*(int*)*(heap->array + index0)) <= (*(int*)*(heap->array + index1)) ? index0 : index1;
	else if (index0 == -1 && index1 != -1)
		return index1;
	else if (index0 != -1 && index1 == -1)
		return index0;
	else
		return -1;
}

// what direction to go if starting from x,y .. and wanting to go to the goal (x_goal, y_goal)
TileCoord DijkstraPathAPI(Graph2D* graph, int x_start, int y_start, int x_goal, int y_goal)
{
	// x_goal and y_goal aren't used for now.. requires saving all possible paths most likely
	// x_goal and y_goal are starting point of dijkstra algorithm
	// find lowest gcost, or parent
	Tile* tile = GraphGetVertex(graph, x_start, y_start);
	TileCoord coord;
	if (tile->parent != NULL)
	{
		coord.x = tile->parent->x_index;
		coord.y = tile->parent->y_index;
	}
	else  // return a value that signifies there is no next tile
	{
		coord.x = 0xFFFF;  // predetermined "invalid" position
		coord.y = 0xFFFF;  // predetermined "invalid" position
	}
	return coord;
}
