package util;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.BEDFileRecord;
import guttmanlab.core.annotation.BlockedAnnotation;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.FeatureCollection;
import guttmanlab.core.coordinatespace.CoordinateSpace;
import guttmanlab.core.sequence.FastaFileIOImpl;
import guttmanlab.core.sequence.Sequence;
import guttmanlab.core.util.CommandLineParser;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.apache.log4j.Logger;

import net.sf.samtools.util.CloseableIterator;

public class SequenceGenerator {
	
	private Map<String, Sequence> chrsByName;
	private Map<String, FeatureCollection<BEDFileRecord>> features;
	private static Logger logger = Logger.getLogger(SequenceGenerator.class.getName());
	
	/**
	 * @param genomeName Name of genome e.g. mm9
	 * @param chrFasta Fasta file of genome
	 * @param featBed Bed file of features e.g. RefSeq annotation
	 */
	private SequenceGenerator(String genomeName, String chrFasta, String featBed) {
		chrsByName = new FastaFileIOImpl().readFromFileByName(chrFasta);
		try {
			features = BEDFileIO.loadFromFileByReferenceName(new File(featBed), CoordinateSpace.forGenome(genomeName));
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	/**
	 * Write regions of DNA excluding any annotated features overlapping the regions
	 * @param dnaRegions Regions to write
	 * @param bothStrands Write DNA on both strands
	 * @param removeEntireSpan Subtract entire span of DNA features
	 * @param outFasta Output fasta file
	 */
	private void writeDnaMinusFeatures(AnnotationCollection<? extends Annotation> dnaRegions, boolean bothStrands, boolean removeEntireSpan, String outFasta) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outFasta));
			CloseableIterator<? extends Annotation> it = dnaRegions.sortedIterator();
			while(it.hasNext()) {
				Annotation dna = it.next();
				Collection<Sequence> toWrite = dnaMinusFeatures(dna, bothStrands, removeEntireSpan);
				for(Sequence seq : toWrite) {
					FastaFileIOImpl.write(seq, bw, 100);
				}
			}
			bw.close();
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	private Collection<Sequence> dnaMinusFeatures(Annotation dnaRegion, boolean bothStrands, boolean removeEntireSpan) {
		
		logger.info("Getting DNA region minus features: " + dnaRegion.getName());
		
		Sequence chr = chrsByName.get(dnaRegion.getReferenceName());
		
		Annotation dnaReg = new BlockedAnnotation(dnaRegion); // Defensive copy
		
		// Subtract each feature from DNA region
		for(BEDFileRecord feature : features.get(dnaRegion.getReferenceName())) {
			if(feature.overlaps(dnaRegion)) {
				Annotation feat = removeEntireSpan ? new SingleInterval(feature.getReferenceName(), feature.getReferenceStartPosition(),
						feature.getReferenceEndPosition(), feature.getOrientation(), feature.getName() + "_span") 
					: new BlockedAnnotation(feature);
				if(bothStrands) feat.setOrientation(Strand.UNKNOWN);
				dnaReg = dnaReg.minus(feat);
			}
		}
		
		// Make sure strand is correct
		if(bothStrands) {
			if(dnaReg.getOrientation().equals(Strand.NEGATIVE) || dnaReg.getOrientation().equals(Strand.POSITIVE)) {
				throw new IllegalArgumentException("If getting DNA on both strands, DNA region strand can't be plus or minus");
			}
		} else {
			if(!dnaReg.getOrientation().equals(Strand.NEGATIVE) && !dnaReg.getOrientation().equals(Strand.POSITIVE)) {
				throw new IllegalArgumentException("If getting DNA on one strand only, DNA region must have defined strand");
			}
		}
		
		Collection<Sequence> rtrn = new ArrayList<Sequence>();
		Iterator<SingleInterval> blocks = dnaReg.getBlocks();
		int blockNum = 0;
		while(blocks.hasNext()) {
			String name = dnaRegion.getName() + "_block" + blockNum;
			SingleInterval block = blocks.next();
			if(bothStrands) {
				SingleInterval blockPlus = new SingleInterval(block.getReferenceName(), block.getReferenceStartPosition(), block.getReferenceEndPosition(), 
						Strand.POSITIVE, name + "_plus");
				SingleInterval blockMinus = new SingleInterval(block.getReferenceName(), block.getReferenceStartPosition(), block.getReferenceEndPosition(), 
						Strand.NEGATIVE, name + "_minus");
				rtrn.add(chr.getSubsequence(blockPlus));
				rtrn.add(chr.getSubsequence(blockMinus));
			} else {
				rtrn.add(chr.getSubsequence(new SingleInterval(block, name)));
			}
			blockNum++;
		}
		
		return rtrn;
		
	}
	
	/**
	 * Write splice junction sequences to fasta file
	 * If the same splice junction appears in multiple transcripts, one is chosen "randomly"
	 * @param flankLen Length of sequence on each side of splice junction
	 * If an exon is shorter than this length, get the whole exon but not part of the next exon outside of that
	 * @param outFile Output file
	 */
	private void writeSpliceJunctions(int flankLen, String outFile) {
		
		logger.info("Writing splice junction sequences to " + outFile);
		
		Map<String, String> seqToName = new HashMap<String, String>();
		int i = 0;
		for(Collection<BEDFileRecord> coll : features.values()) {
			for(BEDFileRecord feat : coll) {
				String featName = feat.getName();
				String chr = feat.getReferenceName();
				Iterator<SingleInterval> blocks = feat.getBlocks();
				SingleInterval prev = blocks.next();
				int jctNum = 0;
				while(blocks.hasNext()) {
					SingleInterval next = blocks.next();
					int start = prev.getReferenceEndPosition() - flankLen;
					int end = next.getReferenceStartPosition() + flankLen;
					Annotation spliceJct = feat.intersect(new SingleInterval(chr, start, end, feat.getOrientation()));
					seqToName.put(chrsByName.get(chr).getSubsequence(spliceJct).getSequenceBases(), "seq_" + i + "_" + featName + "_jct_" + jctNum);
					prev = next;
					jctNum++;
					i++;
				}
			}
		}
		
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));
			for(String seq : seqToName.keySet()) {
				FastaFileIOImpl.write(new Sequence(seqToName.get(seq), seq), bw, 100);
			}
			bw.close();
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
	}
	
	public static void main(String[] args) {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-g", "Genome name", true);
		p.addStringArg("-f", "Genome fasta", true);
		p.addStringArg("-b", "Feature bed", true);
		p.addStringArg("-d", "Bed file of DNA regions to write DNA regions minus features", false, null);
		p.addStringArg("-od", "Output fasta file of DNA regions minus features", false, null);
		p.addBooleanArg("-s", "Use both strands of DNA for writing regions minus features", true);
		p.addBooleanArg("-rs", "Remove entire span of DNA features instead of just blocks", true);
		p.addStringArg("-os", "Output file for splice junction sequences", false, null);
		p.addIntArg("-sl", "Flank length for splice junctions", false, 50);
		p.parse(args);
		
		SequenceGenerator sg = new SequenceGenerator(p.getStringArg("-g"), p.getStringArg("-f"), p.getStringArg("-b"));
		
		String dnaRegions = p.getStringArg("-d");
		String outDnaRegions = p.getStringArg("-od");
		if(dnaRegions != null && outDnaRegions != null) {
			try {
				sg.writeDnaMinusFeatures(BEDFileIO.loadFromFile(dnaRegions, CoordinateSpace.forGenome(p.getStringArg("-g"))), 
						p.getBooleanArg("-s"), p.getBooleanArg("-rs"), outDnaRegions);
			} catch(IOException e) {
				e.printStackTrace();
				System.exit(-1);
			}
		}
		
		String outSpliceJct = p.getStringArg("-os");
		if(outSpliceJct != null) {
			sg.writeSpliceJunctions(p.getIntArg("-sl"), outSpliceJct);
		}
		
	}

}
