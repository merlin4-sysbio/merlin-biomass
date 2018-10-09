package pt.uminho.ceb.biosystems.merlin;

/**
 * @author Oscar Dias
 *
 */
public class Enumerators {

	/**
	 * @author Oscar Dias
	 *
	 */
	public enum ReturnType {

		MMol_GMacromolecule,
		MMol_GDW
	}
	
	public enum MetaboliteGroups {
		
		AA,
		DNA,
		RNA,
		OTHER,
		COFACTOR
	}
	
	public enum MetabolicDataSource{
		
		KEGG("KEGG"),
		
		MODEL_SEED("ModelSEED"){
			@Override
			public String toString(){
				return "ModelSEED metabolic data";
			}
		},
		
//		BIGG("BIGG"){
//			@Override
//			public String toString(){
//				return "BIGG smbl model";
//			}
//		}
		;
		
		private String source;
		
		private MetabolicDataSource(String metabolicDataSource){
			this.source = metabolicDataSource;
		}
		
		public String sourceName(){
			return this.source;
		}
		
		@Override
		public String toString(){
			return "KEGG metabolic data";
		}
	}
}
