package pt.uminho.ceb.biosystems.merlin;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.io.IOUtils;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.RNASequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import pt.uminho.ceb.biosystems.merlin.Enumerators.MetaboliteGroups;
import pt.uminho.ceb.biosystems.merlin.Enumerators.ReturnType;
import pt.uminho.ceb.biosystems.merlin.datatypes.BiomassMetabolite;
import pt.uminho.ceb.biosystems.mew.utilities.datastructures.map.MapUtils;


/**
 * @author Oscar Dias / Sophia Santos
 *
 */
public class EstimateBiomassContents {

	private static final Logger logger = LoggerFactory.getLogger(EstimateBiomassContents.class);


	/**
	 * Get nucleotides relative abundance.
	 * 
	 * @param nucleotideSequences
	 * @param nucCellContent
	 * @param useHydratedMW
	 * @param retType
	 * @param exportFilePath
	 * @param biomassMetabolites
	 * @param rna
	 * @return
	 * @throws Exception
	 */
	public static Map<BiomassMetabolite, Double> getNucleotides_RelativeAbundance(Map<String,AbstractSequence<NucleotideCompound>> nucleotideSequences, double nucCellContent, boolean useHydratedMW, ReturnType retType, String exportFilePath,
			Map<String, BiomassMetabolite> biomassMetabolites, boolean rna) throws Exception {

		Map<BiomassMetabolite, Double> nucGeneData = null, nucG_MolMacromolecule = new HashMap<>(), nucGG_Content = new HashMap<>(), nucMmol_gMacromolecule_Content = new HashMap<>(), nucMmol_gDW_Content = new HashMap<>();

		if(nucleotideSequences == null || nucleotideSequences.isEmpty()) {

			String name = "e-DNA";
			if(rna)
				name = "e-RNA";
			BiomassMetabolite eNuc = new BiomassMetabolite("N", "", name, MetaboliteGroups.OTHER+"");
			nucMmol_gMacromolecule_Content.put(eNuc, -1.0);
			return nucMmol_gMacromolecule_Content;
		}
		else {

			nucGeneData = EstimateBiomassContents.processNucleotides(nucleotideSequences, biomassMetabolites, rna);

			double averageNucleotideMW = 0;
			//Remove water from polymerization
			BiomassMetabolite ppi = biomassMetabolites.get("PPI");
			double ppiMassContents = 0;

			for(BiomassMetabolite nuc : nucGeneData.keySet()) {
				
				if(nucGeneData.get(nuc)>0) {
					
					//Remove ppi for polymerization
					double nucMW = nuc.getMolecularWeight() - ppi.getMolecularWeight();

					double nucMassContent = nucGeneData.get(nuc)*nucMW;
					ppiMassContents += nucGeneData.get(nuc)*ppi.getMolecularWeight();
					
					nucG_MolMacromolecule.put(nuc, nucMassContent);
					averageNucleotideMW+=nucMassContent;
				}
			}
			
			//add macromolecule to list
			String name = "e-DNA";
			if(rna)
				name = "e-RNA";
			BiomassMetabolite eNuc = new BiomassMetabolite("N", "", name, MetaboliteGroups.OTHER+"");
			eNuc.setMolecularWeight(averageNucleotideMW);
			nucG_MolMacromolecule.put(eNuc, -averageNucleotideMW);

			for(BiomassMetabolite nuc : nucG_MolMacromolecule.keySet())
				nucGG_Content.put(nuc, nucG_MolMacromolecule.get(nuc)/averageNucleotideMW);
			

			for(BiomassMetabolite nuc : nucGG_Content.keySet())
				nucMmol_gMacromolecule_Content.put(nuc, (nucGG_Content.get(nuc) * 1000)/(nuc.getMolecularWeight()-ppi.getMolecularWeight()));
			

			for(BiomassMetabolite nuc : nucMmol_gMacromolecule_Content.keySet())
				nucMmol_gDW_Content.put(nuc, (nucMmol_gMacromolecule_Content.get(nuc) * nucCellContent));
			
			//add ppi to equation
			nucG_MolMacromolecule.put(ppi, -ppiMassContents);
			nucGG_Content.put(ppi, nucG_MolMacromolecule.get(ppi)/averageNucleotideMW);
			nucMmol_gMacromolecule_Content.put(ppi, nucGG_Content.get(ppi)*1000/ppi.getMolecularWeight());
			nucMmol_gDW_Content.put(ppi, (nucMmol_gMacromolecule_Content.get(ppi) * nucCellContent));

					
			if(exportFilePath != null) {

				String out = "g/mol\n";
				out += MapUtils.prettyToString(nucG_MolMacromolecule);
				out += "g/gMacromolecule\n";
				out += MapUtils.prettyToString(nucGG_Content);
				out += "mmol/gMacromolecule\n";
				out += MapUtils.prettyToString(nucMmol_gMacromolecule_Content);
				out += "\nmmol/gDW\n";
				out += MapUtils.prettyToString(nucMmol_gDW_Content);

				logger.debug("\n"+out);

				PrintWriter writer = new PrintWriter(exportFilePath, "UTF-8");
				writer.println(out);
				writer.close();
			}

			if(retType.equals(ReturnType.MMol_GMacromolecule))
				return nucMmol_gMacromolecule_Content;
			else if(retType.equals(ReturnType.MMol_GDW)) 
				return nucMmol_gDW_Content;
		}

		return null;
	}

	/**
	 * Process Nucletide.
	 * 
	 * @param nucleotideSequences
	 * @param rna 
	 * @return
	 * @throws Exception 
	 */
	public static Map<BiomassMetabolite, Double>  processNucleotides(Map<String,AbstractSequence<NucleotideCompound>> nucleotideSequences, Map<String, BiomassMetabolite> biomassMetabolites, boolean rna) throws Exception {

//		File file = new File (nucleotideSequences);

		Map <BiomassMetabolite, Double> nucBaseFrequencyMap = new HashMap<>();

//		Map <String, DNASequence> sequencesDNA = FastaReaderHelper.readFastaDNASequence(file);
		
		MetaboliteGroups metaboliteGroup = MetaboliteGroups.DNA;
		if(rna)
			metaboliteGroup = MetaboliteGroups.RNA;

		for(String sequence : nucleotideSequences.keySet()) {

			AbstractSequence<NucleotideCompound> nucSequence = nucleotideSequences.get(sequence);
			if(rna)
				nucSequence = new RNASequence(nucSequence.getSequenceAsString().replaceAll("T", "U"));

			List<NucleotideCompound> nucComp = nucSequence.getCompoundSet().getAllCompounds();

			for (NucleotideCompound nuc : nucComp) {

				BiomassMetabolite nuc_id = null;

				for(String name : biomassMetabolites.keySet()) {

					BiomassMetabolite bm = biomassMetabolites.get(name);

					if(bm.getSingleLetter().equals(nuc.getBase()) && bm.getGroup().equalsIgnoreCase(metaboliteGroup.toString()))
						nuc_id = bm;
				}

				if(nuc_id!= null)
					nucBaseFrequencyMap.put(nuc_id,  new Double(nucSequence.countCompounds(nuc)));
			}
		}

		logger.debug("nuc Bases frequency map {}", nucBaseFrequencyMap);

		return EstimateBiomassContents.getRelativeFrequency(nucBaseFrequencyMap);
	}


	/**
	 * Get proteins relative abundance.
	 * 
	 * @param aaSequences
	 * @param proteinCellContent
	 * @param useHydratedMW
	 * @param retType
	 * @param exportFilePath
	 * @param biomassMetabolites
	 * @param geneData
	 * @param separator
	 * @return
	 * @throws Exception
	 */
	public static Map<BiomassMetabolite, Double> getProteinsRelativeAbundance(Map<String, ProteinSequence> aaSequences, double proteinCellContent, boolean useHydratedMW, ReturnType retType, 
			String exportFilePath, Map<String, BiomassMetabolite> biomassMetabolites, String geneData, String separator) throws Exception {

		Map<BiomassMetabolite, Double> aaMolMolContent = null, aaG_MolMacromolecule = new HashMap<>(), aaGG_Content = new HashMap<>(), 
				aaMmol_gMacromolecule_Content = new HashMap<>(), aaMmol_gDW_Content = new HashMap<>();

		if(aaSequences == null || aaSequences.isEmpty()) {

			BiomassMetabolite eProtein = new BiomassMetabolite("P", "", "e-Protein", MetaboliteGroups.OTHER+"");
			aaMmol_gMacromolecule_Content.put(eProtein, -1.0);
			return aaMmol_gMacromolecule_Content;
		}
		else {

			if (geneData==null || geneData.trim().isEmpty())
				aaMolMolContent = processProteins(aaSequences, biomassMetabolites);
			else
				aaMolMolContent = EstimateBiomassContents.processProteinsExpressionData(aaSequences, biomassMetabolites, geneData, separator);

			double averageProteinMW = 0;
			//Remove water from polymerization
			BiomassMetabolite h2o = biomassMetabolites.get("H2O");
			
			double h2oMassContents = 0;

			for(BiomassMetabolite aa : aaMolMolContent.keySet()) {

				if(aaMolMolContent.get(aa)>0) {
					
					//Remove water from polymerization
					double aaMW = aa.getMolecularWeight() - h2o.getMolecularWeight();

					//				logger.debug(aa.getName()+" "+aa.getMolecularWeight());
					//				logger.debug(h2o.getName()+" "+h2o.getMolecularWeight());
					//				logger.info("");
					
					double aaMassContent = aaMolMolContent.get(aa)*aaMW;
					
					h2oMassContents += aaMolMolContent.get(aa) * h2o.getMolecularWeight();

					aaG_MolMacromolecule.put(aa, aaMassContent);
					averageProteinMW+=aaMassContent;
				}
			}
			
			//add macromolecule to list
			BiomassMetabolite eProtein = new BiomassMetabolite("P", "", "e-Protein", MetaboliteGroups.OTHER+"");
			eProtein.setMolecularWeight(averageProteinMW);
			aaG_MolMacromolecule.put(eProtein, -1*averageProteinMW);

			for(BiomassMetabolite aa : aaG_MolMacromolecule.keySet())
				aaGG_Content.put(aa, aaG_MolMacromolecule.get(aa)/averageProteinMW);

			for(BiomassMetabolite aa : aaGG_Content.keySet()) {
				
				boolean notAvailable = true;
				
				double stoichiometry = (aaGG_Content.get(aa) * 1000)/(aa.getMolecularWeight()-h2o.getMolecularWeight());
				//aaMmol_gMacromolecule_Content.put(aa, stoichiometry);
				
				for(BiomassMetabolite biomassMetabolite : biomassMetabolites.values()) {
					
					if(biomassMetabolite.getSingleLetter().equals(aa.getSingleLetter()+"r") && biomassMetabolite.getGroup().equals(aa.getGroup())) {
						
						notAvailable = false;
						aaMmol_gMacromolecule_Content.put(biomassMetabolite, stoichiometry);
					}				
				
					if(biomassMetabolite.getSingleLetter().equals(aa.getSingleLetter()+"p") && biomassMetabolite.getGroup().equals(aa.getGroup())) {
						
						notAvailable = false;
						aaMmol_gMacromolecule_Content.put(biomassMetabolite, -stoichiometry);
					}
				}
				
				if(notAvailable) 
					aaMmol_gMacromolecule_Content.put(aa, stoichiometry);
			}

			for(BiomassMetabolite aa : aaMmol_gMacromolecule_Content.keySet()) {
				
				double stoichiometry = aaMmol_gMacromolecule_Content.get(aa) * proteinCellContent;
				aaMmol_gDW_Content.put(aa, stoichiometry);
			}
			
			
			//add water to equation
			aaG_MolMacromolecule.put(h2o, -h2oMassContents);
			aaGG_Content.put(h2o, aaG_MolMacromolecule.get(h2o)/averageProteinMW);
			aaMmol_gMacromolecule_Content.put(h2o, aaGG_Content.get(h2o)*1000/h2o.getMolecularWeight());
			aaMmol_gDW_Content.put(h2o, (aaMmol_gMacromolecule_Content.get(h2o) * proteinCellContent));
		

			if(exportFilePath != null) {

				String out = "mol/mol\n";
				out += MapUtils.prettyToString(aaMolMolContent);
				out += "g/gMacromolecule\n";
				out += MapUtils.prettyToString(aaGG_Content);
				out += "mmol/gMacromolecule\n";
				out += MapUtils.prettyToString(aaMmol_gMacromolecule_Content);
				out += "\nmmol/gDW\n";
				out += MapUtils.prettyToString(aaMmol_gDW_Content);

				logger.debug("\n"+out);

				PrintWriter writer = new PrintWriter(exportFilePath, "UTF-8");
				writer.println(out);
				writer.close();
			}

			if(retType.equals(ReturnType.MMol_GMacromolecule))
				return aaMmol_gMacromolecule_Content;
			else if(retType.equals(ReturnType.MMol_GDW)) 
				return aaMmol_gDW_Content;
		}
		return null;
	}

	/**
	 * Process proteins.
	 * 
	 * @param aaSequencesFilePath 
	 * @param biomassMetabolites 
	 * 
	 * @return
	 * @throws Exception
	 */
	public static Map<BiomassMetabolite, Double> processProteins (Map<String, ProteinSequence> aaSequences, Map<String, BiomassMetabolite> biomassMetabolites) throws Exception {

		return processProteinsExpressionData(aaSequences, biomassMetabolites, null, null);
	}

	/**
	 * Process proteins with expression data.
	 * 
	 * @param aaSequencesFilePath
	 * @param biomassMetabolites
	 * @param expressionDataFilePath
	 * @param separator
	 * @return
	 * @throws Exception
	 */
	public static Map<BiomassMetabolite, Double> processProteinsExpressionData(Map<String, ProteinSequence> sequences, Map<String, BiomassMetabolite> biomassMetabolites, String expressionDataFilePath, 
			String separator) throws Exception {

//		File aaSequencesFile= new File (aaSequencesFilePath);

//		Map<String, ProteinSequence> sequences = FastaReaderHelper.readFastaProteinSequence(aaSequencesFile);
		
		Map<String, Double> geneExpressionData = new HashMap<>();

		if(expressionDataFilePath!=null)
			geneExpressionData = readExpressionData_File(new File(expressionDataFilePath), separator);

		Map<BiomassMetabolite, Double> aaBaseFrequencyMap = new HashMap<> ();

		for (String seqName : sequences.keySet()) {

			String geneLocus = parseSequenceName(seqName);
			ProteinSequence proteinSequence = sequences.get(seqName);

			List<AminoAcidCompound> aaComp = proteinSequence.getCompoundSet().getAllCompounds();

			for (AminoAcidCompound aa : aaComp) {

				double expression = 1;

				if(geneExpressionData.containsKey(geneLocus))
					expression = geneExpressionData.get(geneLocus);

				int aaFrequency = proteinSequence.countCompounds(aa);

				BiomassMetabolite aa_id = null;

				for(String name : biomassMetabolites.keySet()) {

					BiomassMetabolite bm = biomassMetabolites.get(name);

					if(bm.getSingleLetter().equalsIgnoreCase(aa.getBase()) && bm.getGroup().equalsIgnoreCase(MetaboliteGroups.AA.toString()))
						aa_id = bm;
				}

				if(aa_id!= null) {

					if (!aaBaseFrequencyMap.containsKey(aa_id))
						aaBaseFrequencyMap.put(aa_id, 0.0);

					double aaFrequencyExpression = aaBaseFrequencyMap.get(aa_id) + (aaFrequency * expression);
					aaBaseFrequencyMap.put(aa_id, aaFrequencyExpression);
				}
			}
		}

		logger.debug("aa Bases frequency map {}", aaBaseFrequencyMap);

		return EstimateBiomassContents.getRelativeFrequency(aaBaseFrequencyMap);
	}

	/**
	 * Get biological compound relative frequency.
	 * 
	 * @param baseFrequencyMap
	 * @return
	 */
	public static Map<BiomassMetabolite, Double> getRelativeFrequency(Map<BiomassMetabolite, Double> baseFrequencyMap) {

		double totalFrequency = 0;
		for (double frequency : baseFrequencyMap.values())
			totalFrequency+=frequency;

		Map<BiomassMetabolite, Double> basesRelativeFrequency = new HashMap<> ();

		for (BiomassMetabolite aa : baseFrequencyMap.keySet()) {

			double frequency = baseFrequencyMap.get(aa);
			basesRelativeFrequency.put(aa, (frequency / totalFrequency));
		}

		return basesRelativeFrequency;		
	}

	/**
	 * Parse sequence name.
	 * 
	 * @param sequenceName
	 * @return
	 */
	public static String parseSequenceName(String sequenceName) {

		return sequenceName.split("\\s+")[0];
	}


	/**
	 * Read expression data file.
	 * 
	 * @param path
	 * @return
	 */
	public static Map<String, Double> readExpressionData_File(File path, String separator) {

		Map<String, Double> geneData = new HashMap<> ();
		try {

			InputStream is = new FileInputStream(path);
			List<String> lines = IOUtils.readLines(is);

			for (int i = 1; i < lines.size(); i++) {

				String line = lines.get(i);
				String[] col = line.split(separator);
				geneData.put(col[0], Double.parseDouble(col[1]));
			}

			IOUtils.closeQuietly(is);

		}
		catch (IOException e) {
			e.printStackTrace();
		}
		return geneData;
	}

	/**
	 * Get cofactors abundance.
	 * 
	 * @param biomassMetabolites
	 * @return
	 */
	public static Map<BiomassMetabolite, Double> getCofactoresAbundance(Map<String, BiomassMetabolite> biomassMetabolites) {

		Map<BiomassMetabolite, Double> ret = new HashMap<BiomassMetabolite, Double>();

		for(String key : biomassMetabolites.keySet()) {

			BiomassMetabolite biomassMetabolite = biomassMetabolites.get(key);
			if(biomassMetabolite.getGroup().equalsIgnoreCase(MetaboliteGroups.COFACTOR.toString()))
				ret.put(biomassMetabolite, 0.000001);
		}

		return ret;
	}

}

