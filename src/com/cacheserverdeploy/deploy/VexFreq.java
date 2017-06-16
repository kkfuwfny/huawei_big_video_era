package com.cacheserverdeploy.deploy;

//网络节点
class VexFreq  implements Comparable<VexFreq>{
	public int vex;                        //网络节点
	public int freq;                       //网络节点出现的频率
    VexFreq(int v){
    	this.vex = v;
    	this.freq = 1;
    }
    @Override
    public int compareTo(VexFreq vf) {
    	if (this.freq > vf.freq)
    		return 1;
    	else 
    		return -1;
    }
}