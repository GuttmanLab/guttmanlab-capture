package capture;

import java.io.IOException;

import org.apache.log4j.Logger;

import general.CommandLineParser;

public class RAPOligoDesigner {

	private static Logger logger = Logger.getLogger(RAPOligoDesigner.class.getName());
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-tfa", "Fasta file of target transcripts", true);
		p.addStringArg("-o", "Output file prefix", true);
		p.addIntArg("-pl", "Probe length", true);
		p.addIntArg("-ps", "Probe step size", true);
		p.addStringArg("-p3", "Full path to primer3core executable", true);
		p.addIntArg("-prl", "Primer length", true);
		p.addDoubleArg("-tm", "Optimal primer Tm", true);
		p.addBooleanArg("-aplc", "Apply probe low complexity filter: remove probes with low complexity sequence", true);
		p.addBooleanArg("-appb", "Apply probe poly base filter: remove probes with strings of single nucleotides", true);
		p.addBooleanArg("-apr", "Apply probe repeat filter: remove probes with repeat masked bases", true);
		p.addBooleanArg("-aprpb", "Apply primer poly base filter: remove primers with strings of single nucleotides", true);
		p.addStringArg("-ppbb", "Bases for probe poly base filter, e.g. ACGT", false, "ACGT");
		p.addIntArg("-ppbc", "For probe poly base filter, max number of repeats of single base within window", false, 12);
		p.addIntArg("-ppbl", "For probe poly base filter, length of windows to check for too many single nucleotides", false, 15);
		p.addStringArg("-prpbb", "Bases for primer poly base filter, e.g. ACGT", false, "ACGT");
		p.addIntArg("-prpbc", "For primer poly base filter, max number of repeats of single base within window", false, 12);
		p.addIntArg("-prpbl", "For primer poly base filter, length of windows to check for too many single nucleotides", false, 15);
		p.addDoubleArg("-prmp", "For probe repeat filter,  max percentage of repeat masked bases", false, 0.07);
		p.addBooleanArg("-prlo", "For probe repeat filter, lower case bases count as repeats", false, false);
		p.addBooleanArg("-prn", "For probe repeat filter, Ns count as repeats", false, false);
		p.parse(args);
		String transcriptFasta = p.getStringArg("-tfa");
		String outFilePrefix = p.getStringArg("-o");
		int probeLength = p.getIntArg("-pl");
		int probeStepSize = p.getIntArg("-ps");
		String primer3core = p.getStringArg("-p3");
		int primerSize = p.getIntArg("-prl");
		double optimalPrimerTm = p.getDoubleArg("-tm");
		boolean applyProbeLowComplexityFilter = p.getBooleanArg("-aplc");
		boolean applyProbePolyBaseFilter = p.getBooleanArg("-appb");
		boolean applyProbeRepeatFilter = p.getBooleanArg("-apr");
		boolean applyPrimerPolyBaseFilter = p.getBooleanArg("-aprpb");
		String probePolyBaseFilterBases = p.getStringArg("-ppbb");
		int probePolyBaseFilterCutoff = p.getIntArg("-ppbc");
		int probePolyBaseFilterRepeatLength = p.getIntArg("-ppbl");
		String primerPolyBaseFilterBases = p.getStringArg("-prpbb");
		int primerPolyBaseFilterCutoff = p.getIntArg("-prpbc");
		int primerPolyBaseFilterRepeatLength = p.getIntArg("-prpbl");
		double probeRepeatFilterMaxPct = p.getDoubleArg("-prmp");
		boolean probeRepeatFilterLower = p.getBooleanArg("-prlo");
		boolean probeRepeatFilterN = p.getBooleanArg("-prn");
		
		OligoPool designer = new OligoPool(transcriptFasta, probeLength, probeStepSize, primer3core, primerSize, 
				optimalPrimerTm, applyProbeLowComplexityFilter, applyProbePolyBaseFilter, probePolyBaseFilterBases, probePolyBaseFilterCutoff, probePolyBaseFilterRepeatLength,
				applyProbeRepeatFilter, probeRepeatFilterMaxPct, probeRepeatFilterLower, probeRepeatFilterN,
				applyPrimerPolyBaseFilter, primerPolyBaseFilterBases, primerPolyBaseFilterCutoff, primerPolyBaseFilterRepeatLength);

		designer.createOligos(outFilePrefix, null, false);
		
		designer.writeFiles(outFilePrefix);
		
		logger.info("");
		logger.info("All done.");
		
	}

}
