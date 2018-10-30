//Steve Monger 17.09.18
//Description: From GTF annotation, output splice site intervals for multi exon genes.
//Output order: donorStartPos(exon1..n-1), DonorEndPos(exon1..n-1), AcceptorStartPos(exon2..n), AcceptorEndPos(exon2..n) 
import java.io.*;

class getSpliceSiteIntervalsFromGTF {

    public static int[] accStartPos = new int[1000];
    public static int[] accEndPos = new int[1000];
    public static int[] donStartPos = new int[1000];
    public static int[] donEndPos = new int[1000];
    public static int exonIndex = 0;
    public static String[] exonNames = new String[1000];
    public static String geneName = "";
    public static int geneStart = 0;
    public static int geneEnd = 0;
    public static String chr = "";
    public static String strand = "";
    public static boolean firstLine = true;

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
        //output previous
        if (!firstLine) {
            printAndReset();
        }
        geneName = sep[8].substring(sep[8].indexOf("gene_name \""));
        geneName = geneName.substring(11, geneName.indexOf("\";"));
        chr = sep[0];
	    geneStart = Integer.parseInt(sep[3]);
	    geneEnd = Integer.parseInt(sep[4]);
        strand = sep[6];
    } else
    if (sep[2].equals("exon")) {
        String exonName = sep[8].substring(sep[8].indexOf("exon_id \""));
        exonName = exonName.substring(9, exonName.indexOf("\";"));
        //check if duplicate exon
        for (int i=0; i<exonIndex; i++) {
            if (exonName.equals(exonNames[i])) {
                break;
            }
        }  
        exonNames[exonIndex]=exonName;
        int start = Integer.parseInt(sep[3]);
        int end = Integer.parseInt(sep[4]);
	    //if positive strand
	    if (strand.equals("+")) {
            //write donor
            donStartPos[exonIndex]=end-2;
            donEndPos[exonIndex]=end+6;
            //write acceptor
            accStartPos[exonIndex]=start-20;
            accEndPos[exonIndex]= start+2;
	    } else
	    //if negative strand
	    if (strand.equals("-")) {
            //write acceptor
            accStartPos[exonIndex]=end-2;
            accEndPos[exonIndex]=end+20;
            //write donor
            donStartPos[exonIndex]=start-6;
            donEndPos[exonIndex]=start+2;
            //donEndPos[exonIndex]=start+1;
        }
        exonIndex++;
	}
    firstLine = false;
}
public static void printAndReset() {
    //prepare comma delimited strings for accStart, accEnd, donStart and donEnd
    String accStart = ""; String accEnd = ""; String donStart = ""; String donEnd = ""; String exonOutput = exonNames[0].concat(",");
    for (int i=1; i<exonIndex; i++) {
        accStart = accStart.concat(Integer.toString(accStartPos[i])).concat(",");
        accEnd = accEnd.concat(Integer.toString(accEndPos[i])).concat(",");
        donStart = donStart.concat(Integer.toString(donStartPos[i-1])).concat(",");
        donEnd = donEnd.concat(Integer.toString(donEndPos[i-1])).concat(",");
        exonOutput = exonOutput.concat(exonNames[i]).concat(",");
    }
    //output splice site intervals for annotated multi-exon genes
    if (exonIndex>1) {
        System.out.println(chr.concat("\t").concat(Integer.toString(geneStart)).concat("\t").concat(Integer.toString(geneEnd)).concat("\t").concat(geneName).concat("\t").concat("GENE").concat("\t").concat(donStart).concat("\t").concat(donEnd).concat("\t").concat(accStart).concat("\t").concat(accEnd).concat("\t").concat(exonOutput));
    }
    //reset variables
    accStartPos = new int[1000]; accEndPos = new int[1000]; donStartPos = new int[1000]; donEndPos = new int[1000];
    exonNames = new String[1000]; geneName = ""; geneStart = 0; geneEnd = 0; chr = ""; exonIndex = 0;
}
}
