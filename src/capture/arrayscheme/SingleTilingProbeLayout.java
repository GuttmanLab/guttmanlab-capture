package capture.arrayscheme;

import org.apache.log4j.Logger;

import nextgen.core.pipeline.ConfigFileOptionValue;

import capture.OligoPool;
import capture.Probe;
import capture.ProbeSet;

import broad.core.sequence.Sequence;


/**
 * @author prussell
 * Probes of a specified length tiled at a specified step size
 */
public class SingleTilingProbeLayout implements ProbeLayout {
	
	private int size;
	private int step;
	private int start;
	private boolean antisense;
	private static Logger logger = Logger.getLogger(SingleTilingProbeLayout.class.getName());
	
	/**
	 * 
	 */
	public SingleTilingProbeLayout() {}
	
	/**
	 * @param probeSize Probe size
	 * @param stepSize Step size
	 * @param startPos First transcript position to include
	 * @param antisenseToTranscript Probe sequence is antisense to transcript
	 */
	public SingleTilingProbeLayout(int probeSize, int stepSize, int startPos, boolean antisenseToTranscript) {
		size = probeSize;
		step = stepSize;
		start = startPos;
		antisense = antisenseToTranscript;
	}
	
	@Override
	public ProbeSet getProbes(Sequence transcript) {
		ProbeSet rtrn = new ProbeSet(name());
		int startPos = start;
		while(startPos + size <= transcript.getLength()) {
			int endPos = startPos + size - 1;
			String probeID = toString() + "_" + transcript.getId() + "_" + startPos + "_" + endPos;
			Probe probe = new Probe(transcript, this, probeID, startPos, size, antisense);
			rtrn.addProbe(probe);
			startPos += step;
		}
		return rtrn;
	}

	@Override
	public String name() {
		return "single_tiling_probe_layout";
	}

	@Override
	public String configFileLineDescription() {
		return OligoPool.probeLayoutOptionFlag + "\t" + name() + "\tprobe_size\tstep_size\tfirst_position\t<'sense' OR 'antisense'>";
	}

	@Override
	public boolean validConfigFileValue(ConfigFileOptionValue value) {
		if(!value.asString(0).equals(OligoPool.probeLayoutOptionFlag)) return false;
		if(!value.asString(1).equals(name())) return false;
		if(value.getActualNumValues() != 6) {
			logger.error("Correct config file line format: " + configFileLineDescription());
			return false;
		}
		try {
			int s = value.asInt(2);
			if(s < 1) {
				logger.error("Correct config file line format: " + configFileLineDescription());
				return false;
			}
		} catch(NumberFormatException e) {
			logger.error("Correct config file line format: " + configFileLineDescription());
			return false;
		}
		try {
			int s = value.asInt(3);
			if(s < 1) {
				logger.error("Correct config file line format: " + configFileLineDescription());
				return false;
			}
		} catch(NumberFormatException e) {
			logger.error("Correct config file line format: " + configFileLineDescription());
			return false;
		}
		try {
			int f = value.asInt(4);
			if(f < 0) {
				logger.error("Correct config file line format: " + configFileLineDescription());
				return false;
			}
		} catch(NumberFormatException e) {
			logger.error("Correct config file line format: " + configFileLineDescription());
			return false;
		}
		String sa = value.asString(5);
		if(!sa.equals("sense") && !sa.equals("antisense")) {
			logger.error("Correct config file line format: " + configFileLineDescription());
			return false;
		}
		return true;
	}

	@Override
	public void setParametersFromConfigFile(ConfigFileOptionValue value) {
		if(!validConfigFileValue(value)) {
			throw new IllegalArgumentException("Config file line invalid. Line format:\n" + configFileLineDescription());
		}
		size = value.asInt(2);
		step = value.asInt(3);
		start = value.asInt(4);
		antisense = value.asString(5).equals("antisense");
	}
	
	@Override
	public String toString() {
		String dir = antisense ? "antisense" : "sense";
		return name() + "_" + size + "_" + step + "_" + start + "_" + dir;
	}

}
