package com.sanger.validateStructuralVariants;

import java.util.*;


public class svInformation {
	private String line = null;
	private String chromLeft = null;
	private Integer startLeft = null;
	private Integer endLeft = null;
	private String chromRight = null;
	private Integer startRight = null;
	private Integer endRight = null;
	private boolean del = false;
	private int count = 0;
	private HashSet<String> reads = new HashSet<String>(100);
	public svInformation(String line, String chromLeft, Integer startLeft, Integer endLeft, String chromRight, Integer startRight, Integer endRight, Boolean del) {
		this.line = line;
		this.chromLeft = chromLeft;
		this.startLeft = startLeft;
		this.endLeft = endLeft;
		this.chromRight = chromRight;
		this.startRight = startRight;
		this.endRight = endRight;
		this.del = del;
	}
	
	public String getLine() {
		return line;
	}
	
	public void addRead(String name) {
		if(!reads.contains(name))
			reads.add(name);
	}
	public Integer size() {
		return reads.size();
	}
	public void increaseCount() {
		count++;
	}
	public Integer getCount() {
		return count;
	}
	public String getChromLeft() {
		return chromLeft;
	}
	public Integer getStartLeft() {
		return startLeft;
	}
	public Integer getEndLeft() {
		return endLeft;
	}
	public String getChromRight() {
		return chromRight;
	}
	public Integer getStartRight() {
		return startRight;
	}
	public Integer getEndRight() {
		return endRight;
	}
	public boolean isDeletion() {
		return del;
	}
}
