package pt.uminho.sysbio;

/**
 * @author Oscar Dias
 *
 */
public class BiomassMetabolite {

	private String singleLetter;
	private String keggId;
	private String name;
	private String group;
	private double molecularWeight;
	private String modelId;
	
	
	/**
	 * @param singleLetter
	 * @param keggId
	 * @param name
	 * @param group
	 */
	public BiomassMetabolite(String singleLetter, String keggId, String name, String group) {
		
		super();
		this.setSingleLetter(singleLetter);
		this.setKeggId(keggId);
		this.setName(name);
		this.setGroup(group);
	}


	/**
	 * @return the singleLetter
	 */
	public String getSingleLetter() {
		return singleLetter;
	}


	/**
	 * @param singleLetter the singleLetter to set
	 */
	public void setSingleLetter(String singleLetter) {
		this.singleLetter = singleLetter;
	}


	/**
	 * @return the keggId
	 */
	public String getKeggId() {
		return keggId;
	}


	/**
	 * @param keggId the keggId to set
	 */
	public void setKeggId(String keggId) {
		this.keggId = keggId;
	}


	/**
	 * @return the name
	 */
	public String getName() {
		return name;
	}


	/**
	 * @param name the name to set
	 */
	public void setName(String name) {
		this.name = name;
	}


	/**
	 * @return the group
	 */
	public String getGroup() {
		return group;
	}


	/**
	 * @param group the group to set
	 */
	public void setGroup(String group) {
		this.group = group;
	}


	/**
	 * @return the molecularWeight
	 */
	public double getMolecularWeight() {
		return molecularWeight;
	}


	/**
	 * @param molecularWeight the molecularWeight to set
	 */
	public void setMolecularWeight(double molecularWeight) {
		this.molecularWeight = molecularWeight;
	}


	/**
	 * @return the modelId
	 */
	public String getModelId() {
		return modelId;
	}


	/**
	 * @param modelId the modelId to set
	 */
	public void setModelId(String modelId) {
		this.modelId = modelId;
	}


	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return "BiomassMetabolite [singleLetter=" + singleLetter + ", keggId=" + keggId + ", name=" + name + ", group="
				+ group + ", molecularWeight=" + molecularWeight + ", modelId=" + modelId + "]";
	}
	
	
}
