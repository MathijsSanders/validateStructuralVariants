package com.sanger.validateStructuralVariants;

import java.util.function.*;

public class svLinkedContainer {
	public svLinkedContainer left = null;
	public svLinkedContainer right = null;
	private svInformation info = null;
	private Function<svInformation,String> getChrom = null;
	private Function<svInformation,Integer> getStart = null;
	private Function<svInformation,Integer> getEnd = null;
	private svLinkedContainer otherCopy = null;
	
	public svLinkedContainer(svInformation info, Function<svInformation,String> getChrom, Function<svInformation,Integer> getStart, Function<svInformation,Integer> getEnd, BiConsumer<svInformation,svLinkedContainer> setCopy) {
		this.info = info;
		this.getChrom = getChrom;
		this.getStart = getStart;
		this.getEnd = getEnd;
		setCopy.accept(info,this);
	}
	public String getChromosome() {
		return getChrom.apply(info);
	}
	public Integer getStart() {
		return getStart.apply(info);
	}
	public Integer getEnd() {
		return getEnd.apply(info);
	}
	public void setOtherCopy(svLinkedContainer otherCopy) {
		this.otherCopy = otherCopy;
	}
	public svLinkedContainer getOtherCopy() {
		return otherCopy;
	}
}
