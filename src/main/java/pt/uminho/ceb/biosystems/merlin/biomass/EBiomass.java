package pt.uminho.ceb.biosystems.merlin.biomass;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.io.FastaReader;
import org.biojava.nbio.core.sequence.io.FastaSequenceParser;
import org.biojava.nbio.core.sequence.io.FileProxyDNASequenceCreator;
import org.biojava.nbio.core.sequence.io.GenericFastaHeaderParser;
import org.biojava.nbio.core.sequence.template.AbstractSequence;

import es.uvigo.ei.aibench.core.operation.annotation.Direction;
import es.uvigo.ei.aibench.core.operation.annotation.Operation;
import es.uvigo.ei.aibench.core.operation.annotation.Port;
import es.uvigo.ei.aibench.workbench.Workbench;
import pt.uminho.ceb.biosystems.merlin.gui.datatypes.WorkspaceAIB;
import pt.uminho.ceb.biosystems.merlin.gui.utilities.LoadFromConf;
import pt.uminho.ceb.biosystems.merlin.gui.utilities.MerlinUtils;
import pt.uminho.ceb.biosystems.merlin.bioapis.externalAPI.utilities.Enumerators.FileExtensions;
import pt.uminho.ceb.biosystems.merlin.biomass.Enumerators.MetabolicDataSource;
import pt.uminho.ceb.biosystems.merlin.biomass.Enumerators.MetaboliteGroups;
import pt.uminho.ceb.biosystems.merlin.biomass.Enumerators.ReturnType;
import pt.uminho.ceb.biosystems.merlin.core.containers.model.CompartmentContainer;
import pt.uminho.ceb.biosystems.merlin.core.utilities.Enumerators.EBIOMASSTEMPLATE;
import pt.uminho.ceb.biosystems.merlin.core.utilities.Enumerators.Pathways;
import pt.uminho.ceb.biosystems.merlin.core.utilities.Enumerators.SequenceType;
import pt.uminho.ceb.biosystems.merlin.core.utilities.Enumerators.SourceType;
import pt.uminho.ceb.biosystems.merlin.datatypes.BiomassMetabolite;
import pt.uminho.ceb.biosystems.merlin.processes.verifiers.CompartmentsVerifier;
import pt.uminho.ceb.biosystems.merlin.services.ProjectServices;
import pt.uminho.ceb.biosystems.merlin.services.model.ModelMetabolitesServices;
import pt.uminho.ceb.biosystems.merlin.services.model.ModelPathwaysServices;
import pt.uminho.ceb.biosystems.merlin.services.model.ModelReactionsServices;
import pt.uminho.ceb.biosystems.merlin.services.model.ModelSequenceServices;
import pt.uminho.ceb.biosystems.merlin.utilities.AmbiguousDNACompoundSet;
import pt.uminho.ceb.biosystems.merlin.utilities.io.FileUtils;

/**
 * @author Oscar Dias
 *
 */
@Operation(name="add Biomass equation", description="add biomass components estimated from the genome sequence")
public class EBiomass {

	private static final Double _DEFAULT_STOICHIOMETRY = 0.1;
	private WorkspaceAIB project;
	private boolean isProtein;
	private boolean isDNA;
	private double proteinContents;
	private double dnaContents;
	private double rnaContents;
	private double rRNA_Contents;
	private double tRNA_Contents;
	private double mRNA_Contents;
	private boolean isRNA;
	private String separator;
	private File geneExpressionFile;
	private String compartment;
	private boolean isGeneExpression;
	private EBIOMASSTEMPLATE template;
	private MetabolicDataSource source;
	private Map<String,ProteinSequence> proteinSequences;
	private Map<String,AbstractSequence<NucleotideCompound>> tRnaSequences, rRnaSequences, mRnaSequences, dnaSequences;
	private boolean isCompartimentalisedModel;



	@Port(direction=Direction.INPUT, name="select contents template",defaultValue = "GramNegative", order=1)
	public void setTemplate(EBIOMASSTEMPLATE template) {

		this.template = template;
	}

	@Port(direction=Direction.INPUT, name="calculate protein contents",validateMethod="useProtein", description= "will use loaded fasta file (protein.faa file)", defaultValue="true",order=2)
	public void isProteinSequences(boolean isProtein) {

		this.isProtein = isProtein;
	}

	@Port(direction=Direction.INPUT, name="calculate DNA contents",validateMethod="useDna", description= "will use loaded file with whole genome (genomic.fna file)", defaultValue="true",order=3)
	public void isDNA_Contents(boolean isDNA) {

		this.isDNA = isDNA;
	}

	@Port(direction=Direction.INPUT, name="calculate RNA contents",validateMethod="useRna", description= "will use loaded files with tRNA and rRNA sequences (.fna files)", defaultValue="true",order=5)
	public void isRNA_Contents(boolean isRNA) {

		this.isRNA = isRNA;
	}

	@Port(direction=Direction.INPUT, name="use Gene expression data?", advanced=true, defaultValue="false",order=13)
	public void isGeneExpression(boolean isGeneExpression) {

		this.isGeneExpression = isGeneExpression;
	}

	@Port(direction=Direction.INPUT, name="gene expression", description="gene expression data (optional)", advanced=true, defaultValue="file path", order=14)
	public void setGeneExpression(File file) {

		this.geneExpressionFile = file;
	}

	@Port(direction=Direction.INPUT, name="gene expression data separator",description="gene expression data separator character (Optional).", advanced=true, defaultValue=";",order=15)
	public void setGeneExpressionSeparator(String separator) {

		this.separator = separator;
	}

	/**
	 * @param project
	 */
	@Port(direction=Direction.INPUT, name="workspace",description="select workspace",validateMethod="validateProject",order=16)

	public void setProject(WorkspaceAIB project) {

		this.project = project;

	}

	/**
	 * @param project
	 */
	@Port(direction=Direction.INPUT, name="source",description="select metabolic data database source",validateMethod="",order=17)

	public void setProject(MetabolicDataSource source) {

		this.source = source;

	}

	@Port(direction=Direction.INPUT, name="biomass compartment",description="compartment for allocating the biomass equations.",defaultValue="auto",order=18, validateMethod="checkBiomassCompartment")
	public void setBiomassCompartment(String compartment) throws Exception {

		this.isCompartimentalisedModel = ProjectServices.isCompartmentalisedModel(this.project.getName());

		Map<String, String> contents = LoadFromConf.loadEbiomassContents(FileUtils.getConfFolderPath(), this.template);
		//		if(compartment.equalsIgnoreCase("auto"))
		//			this.compartment = _BIOMASS_COMPARTMENT;

		try {

			this.proteinContents = Double.parseDouble(contents.get("proteinContents"));
			this.dnaContents = Double.parseDouble(contents.get("dnaContents"));
			this.rnaContents = Double.parseDouble(contents.get("rnaContents"));
			this.mRNA_Contents = Double.parseDouble(contents.get("mRNA_Contents"));
			this.rRNA_Contents = Double.parseDouble(contents.get("rRNA_Contents"));
			this.tRNA_Contents = Double.parseDouble(contents.get("tRNA_Contents"));


			if(this.isRNA && (this.mRNA_Contents + this.rRNA_Contents + this.tRNA_Contents) != 1)
				throw new IllegalArgumentException("the sum of the RNA contents should be equal to 1.");

			if((this.proteinContents + this.rnaContents+ this.dnaContents) > 1)
				throw new IllegalArgumentException("the sum of the macromolecules contents should lower than 1.");

			if (template.equals(EBIOMASSTEMPLATE.Custom) || contents.get("proteinContents").equals("0")) {
				Workbench.getInstance().warn("set custom content values in the configuration file in merlin directory at /conf/ebiomass_contents.conf");
			}
		}
		catch (Exception e) {
			e.printStackTrace();
			Workbench.getInstance().error("error! check contents configuration file in merlin directory at /conf/ebiomass_contents.conf");
		}

		try {

			//file paths
			String preffix = "/out";
			String exportFilePath = FileUtils.getWorkspaceTaxonomyFolderPath(this.project.getName(), this.project.getTaxonomyID()).concat("biomass");

			//create dir
			File f = new File(exportFilePath);
			if(!f.exists())
				f.mkdir();

			//get e-biomass composition
			Map<String, BiomassMetabolite> biomassMetabolites = BiomassUtilities.getBiomassMetabolites(this.source);
			biomassMetabolites = BiomassUtilities.getModelInformationForBiomass(this.project.getName(), biomassMetabolites);

			String geneExpressionPath = null;
			if(this.isGeneExpression)
				geneExpressionPath = this.geneExpressionFile.getPath();

			Map<BiomassMetabolite, Double> averageProtein = EstimateBiomassContents.getProteinsRelativeAbundance(proteinSequences, this.proteinContents, true, ReturnType.MMol_GMacromolecule, exportFilePath+preffix+"_Prot.txt", biomassMetabolites, geneExpressionPath, this.separator );
			Map<BiomassMetabolite, Double> averageDNA = EstimateBiomassContents.getNucleotides_RelativeAbundance(dnaSequences, this.dnaContents, true, ReturnType.MMol_GMacromolecule, exportFilePath+preffix+"_DNA.txt", biomassMetabolites, false);
			Map<BiomassMetabolite, Double> average_rRNA = EstimateBiomassContents.getNucleotides_RelativeAbundance(rRnaSequences, this.rnaContents, true, ReturnType.MMol_GMacromolecule, exportFilePath+preffix+"_rRNA.txt", biomassMetabolites, true);
			Map<BiomassMetabolite, Double> average_tRNA = EstimateBiomassContents.getNucleotides_RelativeAbundance(tRnaSequences, this.rnaContents, true, ReturnType.MMol_GMacromolecule, exportFilePath+preffix+"_tRNA.txt", biomassMetabolites, true);
			Map<BiomassMetabolite, Double> average_mRNA = EstimateBiomassContents.getNucleotides_RelativeAbundance(mRnaSequences, this.rnaContents, true, ReturnType.MMol_GMacromolecule, exportFilePath+preffix+"_mRNA.txt", biomassMetabolites, true);
			Map<BiomassMetabolite, Double> averageRNA = BiomassUtilities.mergeRNAMaps(average_mRNA, this.mRNA_Contents, average_tRNA, this.tRNA_Contents, average_rRNA, this.rRNA_Contents);

			Map<BiomassMetabolite, Double> averageCofactor = EstimateBiomassContents.getCofactoresAbundance(biomassMetabolites);

			// biomass pathway ID
			Map<String, Set<String>> pathway  = new HashMap<>();
			pathway.put(ModelPathwaysServices.addPathway(this.project.getName(), Pathways.BIOMASS), new HashSet<String>());

			//insert data to model
			double lowerBound = 0;
			double upperBound = 999999;



			//Biomass equation
			{
				Map<BiomassMetabolite, Double> averageBiomass = new HashMap<>();

				BiomassMetabolite eProtein = BiomassUtilities.getElementFromMap("e-Protein", averageProtein); 
				if(eProtein.getModelId() == null)
					eProtein.setModelId(ModelMetabolitesServices.insertCompoundToDatabase(this.project.getName(), eProtein.getName(), eProtein.getName(), eProtein.getMolecularWeight()+""));
				averageBiomass.put(eProtein, this.proteinContents);

				BiomassMetabolite eDNA = BiomassUtilities.getElementFromMap("e-DNA", averageDNA);
				if(eDNA.getModelId() == null)
					eDNA.setModelId(ModelMetabolitesServices.insertCompoundToDatabase(this.project.getName(), eDNA.getName(), eDNA.getName(), eDNA.getMolecularWeight()+""));
				averageBiomass.put(eDNA, this.dnaContents);

				BiomassMetabolite eRNA = BiomassUtilities.getElementFromMap("e-RNA", averageRNA);
				if(eRNA.getModelId() == null)
					eRNA.setModelId(ModelMetabolitesServices.insertCompoundToDatabase(this.project.getName(), eRNA.getName(), eRNA.getName(), eRNA.getMolecularWeight()+""));
				averageBiomass.put(eRNA, this.rnaContents);

				BiomassMetabolite cofactor = new BiomassMetabolite("C","e-Cofactor","e-Cofactor", MetaboliteGroups.OTHER.toString());
				if(cofactor.getModelId() == null)
					cofactor.setModelId(ModelMetabolitesServices.insertCompoundToDatabase(this.project.getName(), cofactor.getName(), cofactor.getName(), cofactor.getMolecularWeight()+""));
				averageBiomass.put(cofactor, _DEFAULT_STOICHIOMETRY);
				averageCofactor.put(cofactor, -1.0);

				BiomassMetabolite lipid = new BiomassMetabolite("L","e-Lipid","e-Lipid", MetaboliteGroups.OTHER.toString());
				if(lipid.getModelId() == null)
					lipid.setModelId(ModelMetabolitesServices.insertCompoundToDatabase(this.project.getName(), lipid.getName(), lipid.getName(), lipid.getMolecularWeight()+""));
				averageBiomass.put(lipid, _DEFAULT_STOICHIOMETRY);

				BiomassMetabolite carbohydrate = new BiomassMetabolite("T","e-Carbohydrate","e-Carbohydrate",MetaboliteGroups.OTHER.toString());
				if(carbohydrate.getModelId() == null)
					carbohydrate.setModelId(ModelMetabolitesServices.insertCompoundToDatabase(this.project.getName(), carbohydrate.getName(), carbohydrate.getName(), carbohydrate.getMolecularWeight()+""));
				averageBiomass.put(carbohydrate, _DEFAULT_STOICHIOMETRY);

				BiomassMetabolite eBiomass = new BiomassMetabolite("B", "e-Biomass", "e-Biomass", MetaboliteGroups.OTHER.toString());
				if(eBiomass.getModelId() == null)
					eBiomass.setModelId(ModelMetabolitesServices.insertCompoundToDatabase(this.project.getName(), eBiomass.getName(), eBiomass.getName(), eBiomass.getMolecularWeight()+""));
				//averageBiomass.put(biomassMetabolites.get("e-Biomass"), -1.0);
				averageBiomass.put(eBiomass, -1.0);


				String entity = "e-Biomass";
				String equation = BiomassUtilities.getReactionEquation(averageBiomass, entity);


				List<Integer> aux = ModelReactionsServices.getModelReactionLabelIdByName(this.project.getName(), entity, isCompartimentalisedModel);

				Integer reactionID = -1;
				if(aux != null && aux.size() > 0)
					reactionID = aux.get(0);

				if(reactionID>0) 
					Workbench.getInstance().warn("R_"+entity+" already available in model, skipping reaction");

				else {

					Map<String, String> compartments  = new HashMap<>(), sthoichiometry = new HashMap<String, String>(), chains= new HashMap<String, String>();

					for(BiomassMetabolite bm : averageBiomass.keySet()) {

						compartments.put(bm.getModelId()+"", this.compartment);
						sthoichiometry.put(bm.getModelId()+"", (-1*averageBiomass.get(bm))+"");
						chains.put(bm.getModelId()+"", "0");
					}

					ModelReactionsServices.insertNewReaction_BiomassRelated(this.project.getName(), entity, equation, false, chains, compartments, sthoichiometry, true, pathway, this.compartment,
							false, false, false, lowerBound, upperBound, SourceType.EBIOMASS, null, isCompartimentalisedModel);
				}
			}

			// Cofactor
			{

				String entity = "e-Cofactor";
				String equation = BiomassUtilities.getReactionEquation(averageCofactor, entity);

				String name = ("R_"+entity);

				List<Integer> aux = ModelReactionsServices.getModelReactionLabelIdByName(this.project.getName(), name, isCompartimentalisedModel);

				Integer reactionID = -1;
				if(aux != null && aux.size() > 0)
					reactionID = aux.get(0);


				if(reactionID>0)  
					Workbench.getInstance().warn("R_"+entity+" already available in model, skipping reaction");

				else {

					Map<String, String> compartments  = new HashMap<>(), sthoichiometry = new HashMap<String, String>(), chains= new HashMap<String, String>();

					for(BiomassMetabolite bm : averageCofactor.keySet()) {

						compartments.put(bm.getModelId()+"", this.compartment);
						sthoichiometry.put(bm.getModelId()+"", -1*averageCofactor.get(bm)+"");
						chains.put(bm.getModelId()+"", "0");
					}

					ModelReactionsServices.insertNewReaction_BiomassRelated(this.project.getName(), entity, equation, false, chains, compartments, sthoichiometry, true, pathway, this.compartment,
							false, false, false, lowerBound, upperBound, SourceType.EBIOMASS, null, isCompartimentalisedModel);
				}
			}

			//Protein
			if(this.isProtein) {

				String entity = "e-Protein";
				String equation = BiomassUtilities.getReactionEquation(averageProtein, entity);

				String name = ("R_"+entity);

				List<Integer> aux = ModelReactionsServices.getModelReactionLabelIdByName(this.project.getName(), name, isCompartimentalisedModel);

				Integer reactionID = -1;
				if(aux != null && aux.size() > 0)
					reactionID = aux.get(0);

				if(reactionID>0) {
					Workbench.getInstance().warn("R_"+entity+" already available in model, skipping reaction");
				}
				else {

					Map<String, String> compartments  = new HashMap<>(), metabolites = new HashMap<String, String>(), chains= new HashMap<String, String>();

					for(BiomassMetabolite bm : averageProtein.keySet()) {

						compartments.put(bm.getModelId()+"", this.compartment);

						if (bm.getName().equalsIgnoreCase(entity))
							metabolites.put(bm.getModelId()+"", "1");
						else
							metabolites.put(bm.getModelId()+"", (-1*averageProtein.get(bm))+"");

						chains.put(bm.getModelId()+"", "0");
					}

					ModelReactionsServices.insertNewReaction_BiomassRelated(this.project.getName(), entity, equation, false, chains, compartments, metabolites, true, pathway, this.compartment,
							false, false, false, lowerBound, upperBound, SourceType.EBIOMASS, null, isCompartimentalisedModel);
				}
			}

			//DNA
			if(this.isDNA) {


				String entity = "e-DNA";
				String equation = BiomassUtilities.getReactionEquation(averageDNA, entity);

				String name = ("R_"+entity);

				List<Integer> aux = ModelReactionsServices.getModelReactionLabelIdByName(this.project.getName(), name, isCompartimentalisedModel);

				Integer reactionID = -1;
				if(aux != null && aux.size() > 0)
					reactionID = aux.get(0);

				if(reactionID>0) 
					Workbench.getInstance().warn("R_"+entity+" already available in model, skipping reaction");
				else {

					Map<String, String> compartments  = new HashMap<>(), metabolites = new HashMap<String, String>(), chains= new HashMap<String, String>();

					for(BiomassMetabolite bm : averageDNA.keySet()) {

						compartments.put(bm.getModelId()+"", this.compartment);

						if (bm.getName().equalsIgnoreCase(entity))
							metabolites.put(bm.getModelId()+"", "1");
						else
							metabolites.put(bm.getModelId()+"", (-1*averageDNA.get(bm))+"");

						chains.put(bm.getModelId()+"", "0");
					}

					ModelReactionsServices.insertNewReaction_BiomassRelated(this.project.getName(), entity, equation, false, chains, compartments, metabolites, true, pathway, this.compartment,
							false, false, false, lowerBound, upperBound, SourceType.EBIOMASS, null, isCompartimentalisedModel);
				}
			}

			//RNA
			if(this.isRNA) {

				String entity = "e-RNA";
				String equation = BiomassUtilities.getReactionEquation(averageRNA, entity);
				String name = ("R_"+entity);

				List<Integer> aux = ModelReactionsServices.getModelReactionLabelIdByName(this.project.getName(), name, isCompartimentalisedModel);

				Integer reactionID = -1;
				if(aux != null && aux.size() > 0)
					reactionID = aux.get(0);

				if(reactionID>0) 
					Workbench.getInstance().warn("R_"+entity+" already available in model, skipping reaction");
				else {

					Map<String, String> compartments  = new HashMap<>(), metabolites = new HashMap<String, String>(), chains= new HashMap<String, String>();

					for(BiomassMetabolite bm : averageRNA.keySet()) {

						compartments.put(bm.getModelId()+"", this.compartment);

						if (bm.getName().equalsIgnoreCase(entity))
							metabolites.put(bm.getModelId()+"", "1");
						else
							metabolites.put(bm.getModelId()+"", (-1*averageRNA.get(bm))+"");
						chains.put(bm.getModelId()+"", "0");
					}

					ModelReactionsServices.insertNewReaction_BiomassRelated(this.project.getName(), "e-RNA", equation, false, chains, compartments, metabolites, true, pathway, this.compartment,
							false, false, false, lowerBound, upperBound, SourceType.EBIOMASS, null, isCompartimentalisedModel);
				}
			}

			MerlinUtils.updateAllViews(project.getName());
			Workbench.getInstance().info("e-Biomass equations added to the model!");
		} 
		catch (Exception e) {

			Workbench.getInstance().error("error "+e.getMessage()+" has occured.");
			e.printStackTrace();
		}
	}


	/**
	 * @param project
	 * @throws Exception 
	 */
	@SuppressWarnings("unchecked")
	public void validateProject(WorkspaceAIB project) throws Exception {

		this.project = project;

		if(project == null) {

			throw new IllegalArgumentException("No ProjectGUISelected!");
		}
		else {

			try { 

				if(!ProjectServices.isMetabolicDataAvailable(project.getName()))
					throw new IllegalArgumentException("Please load metabolic data before adding the e-Biomass equations.");

				if(this.isProtein){

					if(ModelSequenceServices.checkGenomeSequences(project.getName(), SequenceType.PROTEIN)) {

						this.proteinSequences =  ModelSequenceServices.getGenomeFromDatabase(this.project.getName(), SequenceType.PROTEIN).entrySet().stream()
								.collect(Collectors.toMap(Map.Entry::getKey, e -> (ProteinSequence)e.getValue()));
					}
					else {
						throw new IllegalArgumentException("Please import protein fasta file ('.faa') to the project for calculating protein contents.");
					}
				}

				if(this.isDNA){

					File genomeFile = new File(FileUtils.getWorkspaceTaxonomyFolderPath(
							this.project.getName(), 
							this.project.getTaxonomyID()).concat(FileExtensions.GENOMIC_FNA.getName()));

					if(genomeFile.exists()){

						AmbiguousDNACompoundSet dnaCompoundSet = new AmbiguousDNACompoundSet();

						dnaCompoundSet.addNucleotideCompound("R", "Y");
						dnaCompoundSet.addNucleotideCompound("Y", "R");
						dnaCompoundSet.addNucleotideCompound("S", "W");
						dnaCompoundSet.addNucleotideCompound("W", "S");
						dnaCompoundSet.addNucleotideCompound("K", "M");
						dnaCompoundSet.addNucleotideCompound("M", "K");

						FastaReader<DNASequence, NucleotideCompound> fastaProxyReader =
								new FastaReader<DNASequence, NucleotideCompound>(
										genomeFile,
										new GenericFastaHeaderParser<DNASequence, NucleotideCompound>(),
										new FileProxyDNASequenceCreator(
												genomeFile,
												dnaCompoundSet,
												new FastaSequenceParser()
												)
										);
						Map<String, DNASequence> mapDNASequence = 	fastaProxyReader.process();

						this.dnaSequences = mapDNASequence.entrySet().stream().collect(Collectors.toMap(Map.Entry::getKey, e -> (AbstractSequence<NucleotideCompound>)e.getValue()));
					}
					else {
						throw new IllegalArgumentException("please import protein fasta file ('genomic.fna') to the project for calculating dna contents");
					}
				}


				if(this.isRNA){

					if(ModelSequenceServices.checkGenomeSequences(project.getName(), SequenceType.TRNA) && ModelSequenceServices.checkGenomeSequences(project.getName(), SequenceType.RRNA)){

						this.tRnaSequences =  ModelSequenceServices.getGenomeFromDatabase(this.project.getName(), SequenceType.TRNA).entrySet().stream()
								.collect(Collectors.toMap(Map.Entry::getKey, e -> (AbstractSequence<NucleotideCompound>)e.getValue()));

						this.rRnaSequences =  ModelSequenceServices.getGenomeFromDatabase(this.project.getName(), SequenceType.RRNA).entrySet().stream()
								.collect(Collectors.toMap(Map.Entry::getKey, e -> (AbstractSequence<NucleotideCompound>)e.getValue()));

					}
					else {

						throw new IllegalArgumentException("please import RNA fasta file ('rna_from_genomic.fna') to the project, for calculating tRNA and rRNA contents");
					}

					//mRNA
					if(ModelSequenceServices.checkGenomeSequences(project.getName(), SequenceType.CDS_DNA)){

						this.mRnaSequences =  ModelSequenceServices.getGenomeFromDatabase(this.project.getName(), SequenceType.CDS_DNA).entrySet().stream()
								.collect(Collectors.toMap(Map.Entry::getKey, e -> (AbstractSequence<NucleotideCompound>)e.getValue()));
					}
					else{
						throw new IllegalArgumentException("please import 'cds_from_genomic.fna' file to the project, for calculating mRNA contents");
					}
				}
			} 
			catch(FileNotFoundException e){

				e.printStackTrace();
			} catch (CompoundNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
		}
	}

	/**
	 * @param contents
	 */
	public void checkProteinContents(double contents) {

		if(contents<0 || contents>1)
			throw new IllegalArgumentException("the contents should be higher than 0 and lower than 1");
		this.proteinContents = contents;
	}
	public void checkDNA_Contents(double contents) {

		if(contents<0 || contents>1)
			throw new IllegalArgumentException("The contents should be higher than 0 and lower than 1!");
		this.dnaContents = contents;
	}
	public void checkRNA_Contents(double contents) {

		if(contents<0 || contents>1)
			throw new IllegalArgumentException("The contents should be higher than 0 and lower than 1!");
		this.rnaContents = contents;
	}
	public void check_mRNA_Contents(double contents) {

		if(contents<0 || contents>1)
			throw new IllegalArgumentException("The contents should be higher than 0 and lower than 1!");
		this.mRNA_Contents = contents;
	}
	public void check_tRNA_Contents(double contents) {

		if(contents<0 || contents>1)
			throw new IllegalArgumentException("The contents should be higher than 0 and lower than 1!");
		this.tRNA_Contents = contents;
	}
	public void check_rRNA_Contents(double contents) {

		if(contents<0 || contents>1)
			throw new IllegalArgumentException("The contents should be higher than 0 and lower than 1!");
		this.rRNA_Contents = contents;
	}


	/**
	 * Validation parameter for assigning the isDna use to the class field.
	 * 
	 * @param isDNA
	 */
	public void useDna(boolean isDNA) {

		this.isDNA = isDNA;

	}

	/**
	 * Validation parameter for assigning the idProtein use to the class field.
	 * 
	 * @param isDNA
	 */
	public void useProtein(boolean isProtein) {

		this.isProtein = isProtein;

	}

	/**
	 * Validation parameter for assigning the isRNA use to the class field.
	 * 
	 * @param isRNA
	 */
	public void useRna(boolean isRna) {

		this.isRNA = isRna;

	}

	/**
	 * @param compartmentID
	 * @throws Exception 
	 */
	public void checkBiomassCompartment(String compartment) throws Exception {

		CompartmentContainer container = CompartmentsVerifier.checkInteriorCompartment(compartment, this.project.getName());

		if(container==null) 
			Workbench.getInstance().warn("No external compartmentID defined!");
		else
			this.compartment = container.getName();

	}
}
