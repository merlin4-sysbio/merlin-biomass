package pt.uminho.ceb.biosystems.merlin;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;

public class Utilities {

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
	 * 
	 * @return
	 */
	public static Map<String, BiomassMetabolite> getBiomassMetabolites () {

		Map<String, BiomassMetabolite> map = new HashMap<>();

		try {

			BufferedReader buf = new BufferedReader(new InputStreamReader(Utilities.class.getResourceAsStream("/biomassIdentifiers.txt")));

			String line;
			while((line = buf.readLine()) != null) {

				if(!line.startsWith("#") && !line.trim().isEmpty()) {

					String[] data = line.split("\t");
					BiomassMetabolite bm = new BiomassMetabolite(data[1], data[2], data[0], data[3]);
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
	public static String getReactionEquation(Map<BiomassMetabolite, Double> entityComposition) {

		String ret = "";

		String reactants = "", products = "";

		for(BiomassMetabolite monomer : entityComposition.keySet()) {

			if(entityComposition.get(monomer)>=0)
				reactants = reactants.concat(entityComposition.get(monomer)+" ").concat(monomer.getName()).concat(" + ");
			else 
				products = products.concat((-entityComposition.get(monomer))+" ").concat(monomer.getName()).concat(" + ");
		}
		
		reactants = reactants.substring(0, reactants.length()-3);
		products = products.substring(0, products.length()-3);

		ret = reactants.concat(" => ").concat(products);

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

}
