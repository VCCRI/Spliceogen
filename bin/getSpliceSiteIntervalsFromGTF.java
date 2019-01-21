	//Steve Monger 17.09.18
	//Description: From GTF annotation, output splice site intervals for multi exon genes.
	//Exon start/end positions which are also transcript start/end positions are ignored, as are duplicate exons.
	//Output order: chr, geneStart, geneEnd, geneName, strand, donorStartPos(exon1..n-1), DonorEndPos(exon1..n-1),
	//AcceptorStartPos(exon2..n), AcceptorEndPos(exon2..n), donorExonID, acceptorExonID.
	import java.io.*;
	import java.util.Arrays;

	class getSpliceSiteIntervalsFromGTF {

	    public static int[] accStartPos = new int[10000];
	    public static int[] accEndPos = new int[10000];
	    public static int[] donStartPos = new int[10000];
	    public static int[] donEndPos = new int[10000];
	    public static int donIndex = 0;
	    public static int accIndex = 0;
	    public static String[] donNames = new String[10000];
	    public static String[] accNames = new String[10000];
	    public static String geneName = "";
	    public static int geneStart = 0;
	    public static int geneEnd = 0;
	    public static String chr = "";
	    public static String strand = "";
	    public static String prevLine = "firstLine";
	    public static boolean duplicateExon = false;

	public static void main(String[] args) {	
	    
	    String line = "";
	    try{
	        BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
	    while ((line = in.readLine()) != null) {
	        try {
	            if (!line.startsWith("#")) {
		            getIntervals(line);
	            }
	        }
	        catch (Exception e){
	        }    
	    }
	    in.close();
	    }catch (Exception e){
	        System.err.println("Error: ".concat(e.getMessage()));
	    }
	}

	public static void getIntervals(String line) {
		String[] sep = line.split("\t");
	    if (sep[2].equals("gene")) {
        	//remove the final exon end pos from donor output for the previous gene
	        if (prevLine.equals("exon")) {
                if (donIndex!=0) {
	                donStartPos[donIndex-1]=0;
                    donEndPos[donIndex-1]=0;
                    donNames[donIndex-1]="";
                    donIndex--;
	            }	            
	        }
	        //output previous
	        if (!prevLine.equals("firstLine")) {
	            printAndReset();
	        }
	        //update variables
	        geneName = sep[8].substring(sep[8].indexOf("gene_name \""));
	        geneName = geneName.substring(11, geneName.indexOf("\";"));
	        chr = sep[0];
		    geneStart = Integer.parseInt(sep[3]);
		    geneEnd = Integer.parseInt(sep[4]);
	        strand = sep[6];      
		    prevLine = "gene";
	    } else
	    if (sep[2].equals("transcript")) {
	        //if multiple transcripts, remove splice sites corresponding to transcript end positions
	        if (prevLine.equals("exon")) {
                if (donIndex!=0) {
                	//removing a donor
	                donStartPos[donIndex-1]=0;
                    donEndPos[donIndex-1]=0;
                    donNames[donIndex-1]="";
                    donIndex--;
	            }	            
	        }
		    prevLine = "transcript";
	    } else
	    if (sep[2].equals("exon")) {
	        String exonName = sep[8].substring(sep[8].indexOf("exon_id \""));
	        exonName = exonName.substring(9, exonName.indexOf("\";"));
	        int start = Integer.parseInt(sep[3]);
	        int end = Integer.parseInt(sep[4]);
	        //check if duplicate exon
	        duplicateExon = false;
	        for (int i=0; i<donIndex; i++) {
	            if (start==donStartPos[i] && end==donEndPos[i]) {
	                duplicateExon = true;
	            }
	        }
	        for (int i=0; i<accIndex; i++) {
	            if (start==accStartPos[i] && end==accEndPos[i]) {
	                duplicateExon = true;
	            }
	        }
	        if (duplicateExon) {
	        	prevLine = "duplicateExon";
	        } else {
		        //if positive strand
		        if (strand.equals("+")) {
	                //write donor
	                donStartPos[donIndex]=end-2;
	                donEndPos[donIndex]=end+6;
	                donNames[donIndex]=exonName;
	                donIndex++;
	                if (!(prevLine.equals("gene")||prevLine.equals("transcript"))) {
	                    //write acceptor
	                    accStartPos[accIndex]=start-20;
	                    accEndPos[accIndex]= start+2;
	                    accNames[accIndex]=exonName;
	                    accIndex++;
	                }
		        } else
		        //if negative strand
		        if (strand.equals("-")) {
	                //write donor
	                donStartPos[donIndex]=start-6;
	                donEndPos[donIndex]=start+2;
	                donNames[donIndex]=exonName;
	                donIndex++;
	                if (!(prevLine.equals("gene")||prevLine.equals("transcript"))) {
	                    //write acceptor
	                    accStartPos[accIndex]=end-2;
	                    accEndPos[accIndex]=end+20;
	                    accNames[accIndex]=exonName;
	                    accIndex++;
	                }
	            }
		        prevLine = "exon";
	        }
		}
	}
	public static void printAndReset() {
	    //prepare comma delimited strings for accStart, accEnd, donStart and donEnd
	    String accStart = ""; String accEnd = ""; String donStart = ""; String donEnd = ""; String accNameOutput = ""; String donNameOutput = "";
	    boolean multiExon = false;
	    for (int i=0; i<accIndex; i++) {
	        multiExon = true;
	        accStart = accStart.concat(Integer.toString(accStartPos[i])).concat(",");
	        accEnd = accEnd.concat(Integer.toString(accEndPos[i])).concat(",");
	        accNameOutput = accNameOutput.concat(accNames[i]).concat(",");
	    }
	    for (int i=0; i<donIndex; i++) {
            multiExon = true;
            donStart = donStart.concat(Integer.toString(donStartPos[i])).concat(",");
            donEnd = donEnd.concat(Integer.toString(donEndPos[i])).concat(",");
            donNameOutput = donNameOutput.concat(donNames[i]).concat(",");
        }
	    //output splice site intervals for annotated multi-exon genes
	    if (multiExon) {
	    	System.out.println(chr.concat("\t").concat(Integer.toString(geneStart)).concat("\t").concat(Integer.toString(geneEnd)).concat("\t").concat(geneName).concat("\t").concat("GENE").concat("\t").concat(strand).concat("\t").concat(donStart).concat("\t").concat(donEnd).concat("\t").concat(accStart).concat("\t").concat(accEnd).concat("\t").concat(donNameOutput.concat("\t").concat(accNameOutput)));
	    }
	    //reset variables
	    Arrays.fill(accStartPos, 0);
	    Arrays.fill(accEndPos, 0);
	    Arrays.fill(donStartPos, 0);
	    Arrays.fill(donEndPos, 0);
	    Arrays.fill(donNames, "");
	    Arrays.fill(accNames, "");
	    geneName = ""; geneStart = 0; geneEnd = 0; chr = ""; donIndex = 0; accIndex = 0;
	}
}



