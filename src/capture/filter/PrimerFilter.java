package capture.filter;

import java.util.List;

import primer.PrimerPair;

import capture.ArrayFeature;
import capture.OligoPool;
import capture.ProbeSet;

/**
 * @author prussell
 *
 */
public interface PrimerFilter extends ArrayFeature {
	
	/**
	 * @param primer Primer pair
	 * @param probes The probes that will potentially be connected to the primer
	 * @return Whether the filter should reject the primer for this probe set
	 */
	public boolean rejectPrimer(PrimerPair primer, List<OligoPool.FullDesignEntry> probes);

}
