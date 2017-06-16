package com.cacheserverdeploy.deploy;

import com.filetool.util.FileUtil;

import java.util.Vector;
import org.omg.CORBA.PUBLIC_MEMBER;

import java.util.Set;
import java.util.Timer;
import java.util.TimerTask;
import java.awt.List;
import java.util.ArrayDeque;
//import java.io.StringBuffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map;
import java.util.Queue;
import java.util.Random;


public class Mydepeloy {
	public static final int N = 1000;             //网络节点最大数
	public static final int INFINITY = (1<<30);   //定义无穷大
	public static final int MAX_LINE_LEN = 55000;
	public static final int MAX_EDGE_NUM = 2000 * 20;
	public static final int GENERATION_LOOP = 500;
	public static final long BEGINE_TIME = System.currentTimeMillis();
	public int totalDemand;                     //消费节点总需求
	public int maxExpense;                      //最大费用
	public int curMinExpense;                   //当前最优费用
	public int nodesNum,edgesNum,usersNum;      //网络节点数、网络链路数、消费节点数
	public int serverExpense;                   //服务器费用
	Boolean firstTime;
	public int edgeOffset;                      //边偏移量
	UserNode []usernodes;                 //消费节点数组
	Vector<Edge> edges = new Vector<>();                 //边节点
	Vector<Integer> [] adj  = new Vector[N];                  //邻接边
	Vector<Integer> [] newAdj = new Vector[N];             //复制的邻接边
	
	Set<Integer> consumerNodes = new HashSet<Integer>();              //消费节点所在的网络节点
	Vector<Integer> serverPos = new Vector<Integer>();               //服务器位置
	Vector<Integer> bestSerPos = new Vector<Integer>();;              //最优服务器位置
	Map<Integer,Integer> nodeMap = new HashMap<>();                //网络节点->消费节点
	Queue<Vector<Integer> > rsPath = new LinkedList<Vector<Integer>>() ;          //保存最优结果路径
	public int POPULATION_NUM;
	String myFilename;
	String topo_file;
	GeneAlg alg = new GeneAlg();
	Boolean needResetE;
	
	public Mydepeloy() {
		for (int i = 0; i < adj.length; ++i) {
			adj[i] = new Vector<>();
			newAdj[i] = new Vector<>();
		}
	}
	

	//你要完成的功能总入口
	public String[] deploy_server(String topo[])
	{
		// myFilename = filename;
	    readData(topo);

	    POPULATION_NUM = 50-nodesNum/20;
	    if(POPULATION_NUM<10)
	        POPULATION_NUM = 10;

	    curMinExpense = maxExpense = usersNum * serverExpense;  //总费用最大的情况是为每个消费节点配置一台服务器

	    for(int i=0;i<usersNum;i++){                            //最差的情况
	        Vector<Integer> v = new Vector<Integer>();
	        v.addElement(usernodes[i].vex);      //消费节点所在的网络节点编号
	        v.addElement(i);                     //第几个消费节点
	        v.addElement(usernodes[i].demand);   //消费节点所要的带宽
	        rsPath.add(v);
	    }

	    genFirstGeneration();   //产生第一代群体，第一代很重要

	    firstTime = true;
	    needResetE = true;
	    
	    //在里面进行基因变异和交叉等操作，选出好的基因，也就是满足客户要求的，成本较低的服务器放置和流量分配方案
	    alg.epoch();

		//writeData();
		//FileUtil.write(topo_file.data(), filename);

		System.out.println("cost：" + curMinExpense );		
		System.out.print("position：");
	    if(curMinExpense<maxExpense){
	        for(int i=0;i<bestSerPos.size();i++){
	        	System.out.print(bestSerPos.get(i) + "  ");
	           
	        }
	        System.out.println("endl");
	    }
	    else{
	        for(int i=0;i<usersNum;i++){
	        	System.out.println(usernodes[i].vex);
	        }
	    }
	    System.out.println("");

	    return writeData();
	}
	
	String [] writeData() {
		StringBuffer buffer = new StringBuffer();
		buffer.append(rsPath.size() + ","); 
		buffer.append("\n\n" + ","); 
	    while( !rsPath.isEmpty() ){
	        Vector<Integer> v = rsPath.poll();      //返回队列头部元素  //移除头元素，并返回该元素值
	        //rsPath.poll();                //移除头元素，并返回该元素值
	        for(int i=0; i < v.size(); i++){
	            if(i != v.size() - 1) {
	            	buffer.append(v.get(i) + " ");
	            	//buffer.append(" ");
	            }
	            else 
	            	buffer.append(v.get(i));
	        }
	        buffer.append("," +"\n" );
	    }
	    String tmpStr [] = buffer.toString().split(",");
	    for(int i= 0; i<tmpStr.length;i++){
	    	 System.out.println(tmpStr[i]);
	    }
	   
	    return buffer.toString().split(",");
	}
	
	
	public void readData(String [] topo) {
		int line_num = topo.length;
	    String[] tempTopo = topo[0].split(" ");  //第1行，有3个元素
		nodesNum = Integer.parseInt(tempTopo[0]);
		edgesNum = Integer.parseInt(tempTopo[1]);
		usersNum = Integer.parseInt(tempTopo[2]);
		
		serverExpense = Integer.parseInt(topo[2]); //第2行，只有一个元素
		
		//这里可以考虑用arrayList替代，以后再改
		usernodes  = new UserNode[usersNum];
		for(int i = 0 ; i < usersNum ; i++){
			usernodes[i] = new UserNode();
		}
		
		int from,to,bandwidth,expense;
		String tempStr [];
	    for(int i=0; i < edgesNum; i++){
	    	tempStr = topo[4+i].split(" "); 
	    	from = Integer.parseInt(tempStr[0]);
	    	to = Integer.parseInt(tempStr[1]);
	    	bandwidth = Integer.parseInt(tempStr[2]);
	    	expense = Integer.parseInt(tempStr[3]);
	    	
	        Edge e1 = new Edge(from,to,bandwidth,expense);
	        Edge e2 = new Edge(to,from,bandwidth,expense);
	        edges.addElement(e1);
	        edges.addElement(e2);
	        int m = edges.size();
	        adj[from].addElement(m-2);
	        adj[to].addElement(m-1);
	    } 
	    
	    totalDemand = 0;
	    int offset = 5 + edgesNum;   //endgesNum是网络节点的信息，一个有这么多条. 加上 5 ，接下来就是消费节点的信息
	    int x;
	    for(int i=0; i<usersNum; i++){
	    	tempStr = topo[offset+i].split(" "); 
	    	x = Integer.parseInt(tempStr[0]);
	    	usernodes[i].vex = Integer.parseInt(tempStr[1]);
	    	//System.err.println("i = " + i + ",  usernodes[i].vex = " + usernodes[i].vex + ",usernodes[i].demand = " + Integer.parseInt(tempStr[2]));
	    	usernodes[i].demand = Integer.parseInt(tempStr[2]);
	 	    totalDemand +=  usernodes[i].demand;
	        nodeMap.put(usernodes[i].vex, i);          //网络节点->消费节点
	        consumerNodes.add(usernodes[i].vex);
	    }
	    
	    
	    for(int i=0;i<usersNum;i++){                //添加源点到每个消费节点的边
	        Edge e = new Edge(nodesNum,usernodes[i].vex,usernodes[i].demand,0);
	        edges.addElement(e);
	        int m = edges.size();
	        adj[nodesNum].addElement(m-1);
	    }
	    edgeOffset = edges.size();  
	}
	
	public void genFirstGeneration() {   
		/*1、第一代基因的第1个基因是全部的消费节点邻接的网络节点都放置一个服务器。 2、2-五分之一是至少布置消费节点数的1/3的服务器节点，其中这些点都是
		 * 频率高的汇点。3、剩下五分之一到全部的都是随机产生的服务器了，至少产生leastNum个
		 */
		
		Vector<VexFreq> vfreq = new Vector<VexFreq>();
		for(int i=0; i < nodesNum; i++){
			vfreq.addElement(new VexFreq(i));
		}
		
		Boolean [] used = new Boolean[N];
		Arrays.fill(used, false);
		
		for(int i=0 ; i<usersNum;i++){
	        spfa_bfs(vfreq, usernodes[i].vex, used);         //每个消费节点找与其最近的消费节点，更新路径上的节点频率
	    }
		
        Collections.sort(vfreq,Collections.reverseOrder());    //降序排序
        
                
        Gene gene0 = new Gene();     //把每个消费节点布置一个服务器的情况加进第一代群体   //网络节点就是所谓基因，有服务器和没有服务器的区分开来，这里用字符0,1（1是有服务器的节点）
        //初始化
        gene0.geneLen = nodesNum;                   
        for(int i=0;i<usersNum;i++){
            gene0.pos.addElement(usernodes[i].vex);
        }
        gene0.valToList();                           //把服务器位置编码成二进制基因序列
        alg.genes.addElement(gene0);
        
        int leastNum = usersNum/3;                    //至少布置消费节点数的1/3
        int i=1;    //第一代基因是gene0，所以这里是i=1
        
        for(i=1;i<POPULATION_NUM/5;i++){              //生成1/5可能更好的染色体
        	Random rand = new Random();
            Gene betterGene = new Gene();
            betterGene.geneLen = nodesNum;           //初始化

            int serverNum;
            int x = usersNum - leastNum;
            serverNum = rand.nextInt(x) + leastNum;         //产生随机服务器个数，至少产生leastNum个
            for(int j=0; j<serverNum; j++){           //优先布置出现频率高的位置   //前面已经算出并排序了
                betterGene.pos.addElement(vfreq.get(j).vex);
            }
            betterGene.valToList();                  //把服务器位置编码成二进制基因序列
            alg.genes.addElement(betterGene);
        }
        
        for(;i<POPULATION_NUM;i++){                   //随机生成初始代基因种群
        	Random rand = new Random();
            Gene gene = new Gene();
            gene.geneLen = nodesNum;                //初始化
            int serverNum;
            int x = usersNum - leastNum;
            serverNum = rand.nextInt(x) + leastNum;         //产生随机服务器个数，至少产生leastNum个
            //System.out.println("i" + i + " genServerPos(serverNum,gene.pos);");
            genServerPos(serverNum,gene.pos);        //还需要随机生成serverNum个位置
            gene.valToList();                        //把服务器位置编码成二进制基因序列
            alg.genes.addElement(gene);
        }
	}
	
		
	public void spfa_bfs(Vector<VexFreq> freqVce, int s, Boolean used[]){
	    Deque<Integer> que = new ArrayDeque<Integer>(); //日后在看
	    Boolean inq [] = new Boolean[N];             //表示是否在队列
	    int D [] = new int[N];
	    int p [] = new int[N];                //上一条弧

	    for(int i=0; i<nodesNum+2; i++){
	        inq[i] = false;
	        p[i] = -1;
	        D[i] = INFINITY;
	    }
	    que.addLast(s);
	    inq[s] = true;
	    D[s] = 0;
	    while(!que.isEmpty()){
	        int k = que.poll();
//	        que.pop_front();           //对头出队
	        inq[k] = false;            //消除标记

	        for(int i=0;i<adj[k].size();i++){               //遍历k的邻接边
	            Edge e = edges.get(adj[k].get(i));
	            if((D[k] + e.cost) < D[ e.to ]){
	                D[e.to] = D[k] + e.cost;               //更新当前与s最短距离
	                p[e.to] = adj[k].get(i);

	                if(!inq[e.to]){   //如果该节点不在队列里面，则添加进去队列后面，标志设为true
	                    inq[e.to] = true;      //加上标记
	                    if(!que.isEmpty() && (D[e.to]<D[que.peek()]))
	                        que.addFirst(e.to);
	                    else que.addLast(e.to);
	                    //System.out.println(e.to);
	                }
	            }
	        }
	    }

	    for(int i=0;i<usersNum;i++){
	        if(!used[usernodes[i].vex] && (p[usernodes[i].vex]!=-1)){    //源点到这个节点有路径并且未当过源节点
	            int k = usernodes[i].vex;
	            while(edges.get(p[k]).from != s){
	                freqVce.get(k).freq++;    //频率加1
	                k = edges.get(p[k]).from;
	            }
	            freqVce.get(s).freq++;
	        }
	    }
	    used[s] = true;
	}
	
	public int spfa_bfs(Set<Integer> path, int s, int t){   //节点s和t是虚拟出来的两个点,s是连接到所有消费节点，所有服务器存放的点指向节点t
	    Deque<Integer> que = new ArrayDeque<Integer>();
	    Boolean inq [] = new Boolean[N];             //表示是否在队列
	    int D [] = new int[N];
	    int flow [] = new int[N];             //表示节点能通过的最大量，即‘瓶颈’
	    int p [] = new int[N];                //存放s到某个点的上一条弧

	    for(int i=0;i<nodesNum+2;i++){
	        inq[i] = false;
	        p[i] = -1;
	        D[i] = INFINITY;
	        flow[i] = INFINITY;
	    }
	    que.addLast(s);
	    inq[s] = true;
	    D[s] = 0;
	    while(!que.isEmpty()){
	        int k = que.poll();
	        inq[k] = false;            //消除标记

	        for(int i=0;i<adj[k].size();i++){               //遍历k的邻接边
	            Edge e = edges.get(adj[k].get(i));
	            if(e.cap > 0 && ((D[k] + e.cost) < D[e.to])){
	                D[e.to] = D[k] + e.cost;               //更新当前与s最短距离
	                p[e.to] = adj[k].get(i);                //存放的是边的序号
	                if(e.to == t){                           //如果下一个点是汇点，直接复制当前节点的flow
	                    flow[e.to] = flow[k];
	                }
	                else{
	                    if(e.cap < flow[k])                //记录下瓶颈
	                        flow[e.to] = e.cap;
	                    else flow[e.to] = flow[k];
	                }

	                if(!inq[e.to]){
	                    inq[e.to] = true;
	                    if(!que.isEmpty() && (D[e.to]<D[que.peek()]))
	                        que.addFirst(e.to);
	                    else que.addLast(e.to);
	                }
	            }
	        }
	    }
	    
	    if(p[t] != -1){ 
	        int k = edges.get(p[t]).from;                   //1、服务器到汇点那段路径不用改，往前一段  //2、从t点开始，逆序回去(t前一个节点->t)
	        while(edges.get(p[k]).from != s){
	            edges.get(p[k]).cap -= flow[t];         //正向边容量减少
	            edges.get(p[k]^1).cap += flow[t];       //反向边容量增加
	            if (path.contains(p[k]^1)) {
	            	if(edges.get(p[k]).cap == edges.get(p[k]^1).cap){        //正反边抵消
	            		path.remove(p[k]^1);
	            	}
	                   
	                else if(edges.get(p[k]).cap < edges.get(p[k]^1).cap){ //正向边流量大，删除反向边，加入正向边
	                    path.remove(p[k]^1);
	                    path.add(p[k]);
	                    
	                    //System.out.println("edges.get(d).from1  : " + edges.get(p[k]).from);
                    	//System.out.println("nodeMape1 : " + nodeMap.get(edges.get(p[k]).from));
	                }
				}	
	        	else{
	                path.add(p[k]);               //把经过的路径保存起来，不包括源到消费节点和服务器到汇的路径
	            }

	            k = edges.get(p[k]).from;
	        }
	        edges.get(p[k]).cap -= flow[t];             //源到消费节点正向边容量减少
//	        if(path.contains(18) && path.contains(6)){
//	        	System.out.println("path.toString() :" + path.toString());
//	        }
	        return flow[t];
	    }

	    return -1;

	}

	
	
	int solve(int n){		
	    Set<Integer> path = new HashSet<Integer>();                      //保存经过的路径
	    initEdges();                        //初始化工作的图
	    int demand = totalDemand;
	    int sum = n * serverExpense;        //记录该轮产生的总费用，初始化为服务器的开销

	    while(demand>0){
	    	
	        int f = spfa_bfs(path, nodesNum, nodesNum+1);	 
	        if(f!=-1){
	            demand -= f;
	        }
	        else break;
	    }

	    if(demand==0){
	    	Iterator<Integer> it = path.iterator();	
	        int edge_pos;
	        int flow;
	        while(it.hasNext()){   //遍历
	        	edge_pos = it.next();
	            flow = Math.abs(edges.get(edge_pos).cap - edges.get(edge_pos^1).cap) / 2;
	            sum += flow * edges.get(edge_pos).cost;	           
	        }
	        if(sum < curMinExpense){
	            System.out.println(sum);
	        	needResetE = false;                 //不需要重置每条边
	            curMinExpense = sum;                //记下最小花费
	            bestSerPos.clear();
	            bestSerPos = (Vector<Integer>) serverPos.clone();    //记下最优服务器位置
	            //bestSerPos = serverPos;      
	            
	            genPath(path);                      //构建路径
	           
	        }
	        else {
	            needResetE = true;
	        }
//	        System.out.println("int solve()   return sum Done. ");
	        return sum;
	    }

	    needResetE = true;
//	    System.out.println("int solve()  return -1 Done. ");
	    return -1;

	}
	
	void genPath(Set<Integer> path){	
		    while(!rsPath.isEmpty())                      //清空
		        rsPath.poll();
		    
		    
		    int p;
		    int flow,tmpf;		    
		    Iterator<Integer> it;	
		    int [] link = new int[bestSerPos.size()];     //保存与服务器位置有链接的边数

		    for(int i=0;i<bestSerPos.size();i++){
		        if(reachConsumer(bestSerPos.get(i))){
		            Vector<Integer> v = new Vector<Integer>();
		            v.addElement(bestSerPos.get(i));		            
		            v.addElement(nodeMap.get(bestSerPos.get(i)));
		            v.addElement(usernodes[nodeMap.get(bestSerPos.get(i))].demand);
		            rsPath.add(v);		            
		        }
		        link[i] = 0;
		        it = path.iterator();
		        while(it.hasNext()){
		            int d = it.next();
		            if(edges.get(d).to == bestSerPos.get(i)){
		            	//System.err.println("from = " + edges.get(d).from + ", edges.get(d).to " + edges.get(d).to);
		                link[i] ++;
		            }
		        }
		    }

		    for(int i=0;i<bestSerPos.size();i++){
		    	while(link[i]>0){
		    		Vector<Integer> es = new Vector<Integer>();
		            p = bestSerPos.get(i);
		            Vector<Integer> v = new Vector<Integer>();
		            v.addElement(p);
		            flow = INFINITY;
		            Boolean flag = true;

		            while(flag){			            		
		            	flag = false;
		            	it = path.iterator();
		                while(it.hasNext()){
		                    int d = it.next();
		                   
		                    if((edges.get(d).to == p) && ((edges.get(d).cap)  < (edges.get(d^1).cap) )) {              //找到对应的边		                    			                    	
		                    	//System.out.println("edges.get(d).from =" + edges.get(d).from + ", edges.get(d).to =" + edges.get(d).to +  "  ,edges.get(d).cost =" + edges.get(d).cost);
		                    	v.addElement(edges.get(d).from);
		                    	//System.out.println("d = " + d +  ", " + "edges.get(d).from=" + edges.get(d).from);
		                        tmpf = Math.abs(edges.get(d).cap - edges.get(d^1).cap) / 2;    //这条边的流量
		                        es.addElement(d);                            //记录下这条边的位置，用来修改流量
		                        if(tmpf < flow){                            //记下最小的流量
		                            flow = tmpf;
		                        }
		                        p = edges.get(d).from;
		                        flag = true;
		                        /*if(v.firstElement() == 1 && v.contains(2)){
		                        	System.err.println("begin find edges:");
		                        	for(int iii = 0; iii < edges.size() ;iii++){
		                        		if(edges.get(iii).from == 2 || edges.get(iii).to == 2){
		                        			System.err.println("iii = " + iii + "test from =" + edges.get(iii).from + ", to=" + edges.get(iii).to);
		                        			System.out.println("path.toString() :" + path.toString());
		                        		}
		                        	}
		                        }*/
		                        int tempmm = 0;
		    		            while(tempmm <=  v.size() - 1){
		    		            	 System.out.print( v.get(tempmm) + " "); 
		    		            	 tempmm++;
		    		            }		    		            
		                        break;
		                    }
		                }
		                /*
		                if(!it.hasNext()){                                 //看是否已经找到头了，是的话退出循环
		                    flag = false;
		                    //System.out.println("i = " + i + "  , solve()  flge = false  Done. ");
		                }	
		                */	              
		            }
		            
		            //分别添加消费节点的编号 和 该路线所能提供的带宽到v
		            int m = v.size()-1;
		            
		            
		            v.addElement(nodeMap.get(v.get(m)));                             //把该网络节点下所对应的消费节点放进去
		            
		            //System.err.println("m = " + m +", nodeMap.get(v.get(m)) = " + nodeMap.get(v.get(m)) + ", v.get(m) = " + v.get(m));
		            v.addElement(flow);                                      //这条路径提供的的流量
		            System.out.println("消费节点 ： " + v.get(v.size()-2) + "   ****   流量 ： " + v.get(v.size()-1));  
		            
		            int tempmm = 0;
		            while(tempmm <=  v.size() - 2){
		            	 System.out.print( v.get(tempmm) + " "); 
		            	 tempmm++;
		            }
		            System.out.println("*****");
		            
		            rsPath.add(v);
		            for(int j=0;j<es.size();j++){                           //修改路径每条边的流量
		                edges.get(es.get(j)).cap += flow;
		                edges.get(es.get(j)^1).cap -= flow;
		                if(edges.get(es.get(j)).cap==edges.get(es.get(j)^1).cap){        //跟原图的边比较，相等说明这条边已经用完
		                    path.remove(es.get(j));                             //删除这条边
		                    if(edges.get(es.get(j)).to==bestSerPos.get(i)){
		                        link[i]--;
		                    }
		                }
		            }
		        }		      
		    }
		   // delete[] link;//有疑问适合要想处理，日后再说
		}

	
	Boolean reachConsumer(int k){                        //判断服务器是否直接放置在消费节点
		return consumerNodes.contains(k);
	}
	
	public void initEdges(){                                 //初始化原图，在原图的基础上操作
	    if(firstTime){                                //第一次不用重置
	        firstTime = false;
	    }
	    else{
	        int c,i,j;
	        for(i=edgeOffset;i<edges.size();i++){     //1.重置邻接表    //即，删除跟所有服务器虚拟出来的一个汇点的边，每个服务器都只有一条和汇点相连的边。
	            int k = edges.get(i).from;
	            //System.out.println(k + "<-k :"+adj[k].size());
	            adj[k].removeElementAt(adj[k].size() - 1);   
	            //adj[k].removeAllElements();
	        }
	        //日后在看//
	        if(needResetE){                               //需要重置每条边
	            for(i=0;i<edgesNum*2-1;i+=2){             //重置每条边的容量
	                c = (edges.get(i).cap + edges.get(i^1).cap) / 2;
	                edges.get(i).cap = c;
	                edges.get(i^1).cap = c;
	            }
	        }

	        for(j=0,i=edgesNum*2;i<edgeOffset;i++,j++){        //重置源节点到每个消费节点的容量
	        	 edges.get(i).cap = usernodes[j].demand;
	        }
	        /*
	        int tempS = edges.size();
	        for(int temp = edgeOffset; temp < edges.size(); temp++){
	        	edges.removeElementAt(temp);
	        }
	        */
	        for(int temp = edges.size() - 1;   temp >= edgeOffset;temp--){
	        	edges.removeElementAt(temp);
	        }
	        
	        //日后再说
	        //edges.i
//	        Iterator<Edge> it = edges.iterator();
//	        while(it.hasNext())
//	        edges.removeElementAt(it.);         //删除上一版的图中服务器节点到汇节点的边
	    }


	    for(int i=0;i<serverPos.size();i++){                           //给现版本服务器节点加上到汇点的边，*汇点为nodesNum+1*
	        Edge e = new Edge(serverPos.get(i),nodesNum+1,INFINITY,0);
	        edges.addElement(e);
	        int m = edges.size();
	        adj[serverPos.get(i)].addElement(m-1);
	    }
	}
	
	public void genServerPos(int n, Vector<Integer> pos){
	    int i;
	    Random rand = new Random();
	    while(pos.size()<n){
	        int k = rand.nextInt(nodesNum);            //生成不重复的随机节点
	        for(i=0;i<pos.size();i++){           //判断是否存在重复的，直到不重复为止
	            if(pos.get(i)==k)
	                break;
	        }
	        if(i>=pos.size())
	            pos.addElement(k);
	    }
	}
	
	
	class Gene {
		public static final int INFINITY = (1<<30);
		public double fitness = 0.0;                            //适应值
		public int geneLen;                                     //基因序列长度
		Vector<Integer> pos = new Vector<Integer>() ;                                   //服务器位置
		public char [] geneList ;                  //服务器位置对应的基因序列  

		public Gene copyGene(){                                 //复制基因 
			Gene newGene = new Gene();
			newGene.fitness = this.fitness;
			newGene.geneLen = this.geneLen;
			newGene.geneList = null;
			newGene.pos = (Vector<Integer>) this.pos.clone();
			newGene.geneList = new char[geneLen];
			
			for(int i=0;i<this.geneLen;i++){
				newGene.geneList[i] = this.geneList[i];
			}
			
			//newGene.geneList = this.geneList.clone();
			return newGene;
		}       

		public  void computeFitness(){                  //计算基因适应值
			if(this.pos.size() > usersNum){       //如果这个基因中的服务器比我消费节点还多，也就是比第一代基因还差，那么这个基因是很差的，fitness也就给的几乎为0 了
				this.fitness = 1.0/INFINITY;
				return;
			}
			//serverPos.clear();
			serverPos = (Vector<Integer>) this.pos.clone(); //先将该基因的服务器位置赋给全局变量
			int expense = solve(serverPos.size());       
			if(expense == -1){                       //为-1 说明这个服务放置的位置不满足各个用户的流量带宽需求
				this.fitness = 1.0/INFINITY;
			}
			else{
				this.fitness = 1.0/expense;
			}
		}
		
		                                                
		public  void listToVal(){                       //把二进制基因序列解码成服务器位置
			this.pos.clear();
		    for(int i=0;i<geneLen;i++){
		        if(geneList[i] == '1'){
		            this.pos.addElement(i);		          
		        }
		        	
		    }		    
		}	

		public  void valToList(){                       //把服务器位置编码成二进制基因序列
			if(geneList  == null){
		        geneList = new  char [geneLen];
		    }

		    for(int i=0; i<geneLen; i++){
		        geneList[i]= '0';		      
		    }
		    for(int i=0;i<this.pos.size();i++){            //每个服务器位置设为1
		        geneList[this.pos.get(i)] = '1';
		    }
		}
	}
	
	

class GeneAlg {
	public int generation = 1;                     //当前的基因是第几代
	public double totalFiness = 0.0;               //一代基因的总适应值
	public double mutateRate = 0.15;               //基因突变概率
	public double crossRate = 0.7;                 //基因重组概率
	Gene bestGene = null;                          //最优基因
    Vector<Gene> genes = new  Vector<Gene>();                            //当前代基因
    Vector<Integer> hammings = new Vector<Integer>();
	
    public void computeBestGene(){                  //计算一代中的最优基因    	
    	int bestIndex = 0;
    	double bestFitness = 0.0;
        for(int i=0; i<POPULATION_NUM; i++){
            if(genes.get(i).fitness > bestFitness){
                bestIndex = i;
                bestFitness = genes.get(i).fitness;
            }
        }
        bestGene = genes.get(bestIndex);        
    }
    
    public void crossOver(Gene best,Gene gene2){   //基因重组    	
    	Random rand =new Random();    
    	int pos = rand.nextInt(nodesNum);     //随机找一个位置
        for(int i=pos;i<nodesNum;i++){                      //从该位置开始把更优基因的因子复制到另一个基因
            if(best.geneList[i] != gene2.geneList[i]){
            	gene2.geneList[i] = best.geneList[i];                
            }
        }     
    }
    
    public void mutate(Gene gene){                  //基因突变
    	Random rand =new Random(); 
    	int pos = rand.nextInt(nodesNum);                         //随机找一个位置
        if(gene.geneList[pos] == '0')
            gene.geneList[pos] = '1';
        else  gene.geneList[pos] = '0';
    	
    }
    
    public Gene getRoulette(){            //轮盘赌，随机选出一个基因个体
    	Random rand =new Random(); 
    	double rd = rand.nextInt(100000)/100000.0;
        double slice = rd * totalFiness;
        double cursor = 0.0;
        int i;
        for(i=0; i< POPULATION_NUM; i++){
            cursor += genes.get(i).fitness;            
            if(cursor>slice){
                break;
            }
        }
        
        Gene newGene = genes.get(i).copyGene();                //复制一份该基因

        return newGene;
    }
    
//    public int hammingDist(int k){}
    	
    public void epoch(){
    	Vector<Gene> newGenes = new Vector<Gene>();
    	double bestFitness = 0.0;
    	long CURRENT_TIME = 0;  //记录当前的时间，程序到一定时间退出循环
    	 do{
    	        System.out.println("第" + generation + "代");
    	        totalFiness = 0.0;
    	        for(int i=0;i<POPULATION_NUM;i++){               //计算每个个体的适应值
    	            Gene g = genes.get(i);                       //取出基因
    	            g.computeFitness();                         //调用solve()函数
    	            totalFiness += g.fitness;                   //累加适应值
    	        }
    	        computeBestGene();                               //计算最优基因
    	        //cout<<this->bestGene->fitness<<endl;

    	        newGenes.clear();                                //初始化新一代种群
    	        Gene bestG = this.bestGene.copyGene();        //复制一份最优基因直接遗传到下一代
    	        newGenes.addElement(bestG);
    	        Random rand = new Random();
    	        while(newGenes.size()<POPULATION_NUM){           //按照轮盘赌算法选出一个个体与最优基因交叉
    	            Gene gene2 = getRoulette();
    	            
    	            if((rand.nextInt(100000)/100000.0) < crossRate){   //对2个个体按交叉概率执行交叉操作
    	                crossOver(bestG,gene2);
    	            }
    	            if((rand.nextInt(100000)/100000.0) < mutateRate){  //对个体按突变概率执行突变操作
    	                mutate(gene2);
    	            }
    	            
    	            gene2.listToVal();       //把二进制基因序列解码成服务器位置
    	            
    	            if(newGenes.size() < POPULATION_NUM && gene2.pos.size() < usersNum)
    	                newGenes.addElement(gene2);
    	            
    	        }
    	        
    	        genes = (Vector<Gene>) newGenes.clone();                                //新的种群取代旧的种群
    	        generation++;
    	        CURRENT_TIME = System.currentTimeMillis() - BEGINE_TIME;
    	    }while(generation < 20 && CURRENT_TIME <= 100000);//while(generation < GENERATION_LOOP && CURRENT_TIME <= 100000);
    	 System.err.println("we use time is " + CURRENT_TIME);
    }  

}
	
}
