package pt.uminho.ceb.biosystems.merlin.biomass;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import pt.uminho.ceb.biosystems.merlin.biomass.Enumerators.MetabolicDataSource;
import pt.uminho.ceb.biosystems.merlin.core.containers.model.MetaboliteContainer;
import pt.uminho.ceb.biosystems.merlin.datatypes.BiomassMetabolite;
import pt.uminho.ceb.biosystems.merlin.services.model.ModelMetabolitesServices;

public class BiomassUtilities {

	private static final String BIOMASS_IDENTIFIERS = "/biomassIdentifiers.txt";

	/**
	 * Merge RNA maps.
	 * 
	 * @param mRNA
	 * @param mRnaContents 
	 * @param tRNA
	 * @param tRnaContents 
	 * @param rRNA
	 * @param rRnaContents 
	 * @return
	 */
	public static Map<BiomassMetabolite, Double> mergeRNAMaps(Map<BiomassMetabolite, Double> mRNA, double mRnaContents, Map<BiomassMetabolite, Double> tRNA, double tRnaContents, Map<BiomassMetabolite, Double> rRNA, double rRnaContents) {

		Map<BiomassMetabolite, Double> ret = new HashMap<>();
		
		double e_rRNA = 0, e_mRNA = 0, e_tRNA = 0;

		for(BiomassMetabolite mKey : mRNA.keySet()) {

			if(mKey.getName().equals("e-RNA"))
				e_mRNA = mKey.getMolecularWeight();
			
			double rRNA_Contents = 0;
			for(BiomassMetabolite rKey : rRNA.keySet()) 
				if(rKey.getName().equals(mKey.getName())) {
					
					rRNA_Contents = rRNA.get(rKey);
					
					if(rKey.getName().equals("e-RNA"))
						e_rRNA = rKey.getMolecularWeight();
				}
			
			double tRNA_Contents = 0;
			for(BiomassMetabolite tKey : tRNA.keySet())
				if(tKey.getName().equals(mKey.getName())) {
					
					tRNA_Contents = tRNA.get(tKey);
					
					if(tKey.getName().equals("e-RNA"))
						e_tRNA = tKey.getMolecularWeight();
				}
			
			double total = mRnaContents*mRNA.get(mKey) + rRnaContents*rRNA_Contents + tRnaContents*tRNA_Contents;
			
			if(mKey.getName().equals("e-RNA"))
				mKey.setMolecularWeight(mRnaContents*e_mRNA + rRnaContents*e_rRNA + tRnaContents*e_tRNA);

			ret.put(mKey, total);
		}
		return ret;
	}
	
	/**
	 * Build initial information for e-biomass.
	 * @param source 
	 * 
	 * @return
	 */
	public static Map<String, BiomassMetabolite> getBiomassMetabolites (MetabolicDataSource source) {

		Map<String, BiomassMetabolite> map = new HashMap<>();

		try {

			BufferedReader buf = new BufferedReader(new InputStreamReader(BiomassUtilities.class.getResourceAsStream(BIOMASS_IDENTIFIERS)));

			String line;
			while((line = buf.readLine()) != null) {

				if(!line.startsWith("#") && !line.trim().isEmpty()) {

					String[] data = line.split("\t");
					
					String metaboliteRefId;
					if(source.equals(MetabolicDataSource.MODEL_SEED))
						metaboliteRefId = data[3];
					else
						metaboliteRefId = data[2];
					
					BiomassMetabolite bm = new BiomassMetabolite(data[1], metaboliteRefId, data[0], data[4]);
					map.put(data[0], bm);
				}
			}
			
			buf.close();

		}
		catch(Exception e) {

			e.printStackTrace();
		}
		return map;
	}

	/**
	 * @param entityComposition
	 * @param entity
	 * @return
	 */
	public static String getReactionEquation(Map<BiomassMetabolite, Double> entityComposition, String macromolecule) {

		String ret = "";

		String reactants = "", products = "";

		for(BiomassMetabolite monomer : entityComposition.keySet()) {

			if(entityComposition.get(monomer)>=0)
				reactants = reactants.concat(entityComposition.get(monomer)+" ").concat(monomer.getName()).concat(" + ");
			else if (monomer.getName().equalsIgnoreCase(macromolecule))
				products = products.concat("1 ").concat(monomer.getName()).concat(" + ");
			else
				products = products.concat((-entityComposition.get(monomer))+" ").concat(monomer.getName()).concat(" + ");
		}
		
		System.out.println(macromolecule +" re " + reactants.length());
		System.out.println(entityComposition);
		if(reactants.length()>0 && products.length()>0) {
			reactants = reactants.substring(0, reactants.length()-3);
			products = products.substring(0, products.length()-3);
			ret = reactants.concat(" => ").concat(products);
			return ret;
		}
		
		return ret;
	}

	/**
	 * Return BiomassMetabolite from map using name;
	 * 
	 * @param name
	 * @param map
	 * @return
	 */
	public static BiomassMetabolite getElementFromMap(String name, Map<BiomassMetabolite, Double> map) {

		for(BiomassMetabolite bm : map.keySet())
			if(bm.getName().equals(name))
				return bm;

		return null;
	}
	
	/**
	 * Get information for e-biomass.
	 * 
	 * @param databaseName
	 * @param data
	 * @return
	 * @throws Exception 
	 */
	public static Map<String, BiomassMetabolite> getModelInformationForBiomass(String databaseName, Map<String, BiomassMetabolite> data) throws Exception {

		List<String> keggs = new ArrayList<>();

		for(String name : data.keySet())
			keggs.add(data.get(name).getKeggId());
		
		for(String name : data.keySet()) {
			
			String kegg = data.get(name).getKeggId();
			
			MetaboliteContainer compound = ModelMetabolitesServices.getCompoundByExternalIdentifier(databaseName, kegg);
			
			if(compound != null){
				data.get(name).setModelId(compound.getMetaboliteID());
				if(compound.getMolecular_weight() != null && !compound.getMolecular_weight().isEmpty())
					data.get(name).setMolecularWeight(Double.valueOf(compound.getMolecular_weight()));
			}
		}
		return data;
	}
}
