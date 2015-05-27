//Author : Jayant Kumar, jayant@umiacs.umd.edu, UMD college park

#ifndef __GRAPH__
#define __GRAPH__
#include <iostream>
#include <vector>
#include "stdafx.h"

using namespace std;

struct node {
int info;
node *next;
};

class Queue {
	public:
	Queue();
	~Queue();
	bool isEmpty();
	void add(int);
	int get();
	private:
	node *first, *last;
};

class Graph {
	public:
		Graph(int size = 2);
		~Graph();
		bool isConnected(int, int);
		// adds the (x, y) pair to the edge set
		void addEdge(int x, int y);
		// performs a Breadth First Search starting with node x
		vector<int> BFS(int x);
		// searches for the minimum length path
		// between the start and target vertices
		//void minPath(int start, int target);
	private :
		int n;
		int **A;
};
#endif 
