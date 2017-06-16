package com.cacheserverdeploy.deploy;

public class Edge {
	public int from;
	public int to;
	public int cap;
	public int cost;
	Edge(int u, int v, int cap, int cost){
		this.from = u;
		this.to = v;
		this.cap = cap;
		this.cost = cost;		
	}
}
