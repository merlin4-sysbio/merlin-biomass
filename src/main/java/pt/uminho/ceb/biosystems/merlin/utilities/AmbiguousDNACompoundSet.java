package pt.uminho.ceb.biosystems.merlin.utilities;

import org.biojava.nbio.core.sequence.compound.DNACompoundSet;

public class AmbiguousDNACompoundSet extends DNACompoundSet {

	public AmbiguousDNACompoundSet() {
		
		super();
		
	}
	
	public void addNucleotideCompound(String base, String complement) {
		
		super.addNucleotideCompound(base, complement);
	}
	
}
