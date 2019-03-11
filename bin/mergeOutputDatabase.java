//Description: Merging GeneSplicer, MaxEntScan and ESRseq scores along with annotation information to be outputted together on one line per variant.
import java.io.*;
import java.util.Arrays;
import java.util.jar.Attributes.Name;
import java.lang.*;
import java.text.DecimalFormat;
public class mergeOutputDatabase {

	//store output lines in buffers to minimise I/O
    public static String[] avBuffer = new String[30000];
    public static String[] donBuffer = new String[30000];
    public static String[] accBuffer = new String[30000];
    public static String[] ssBuffer = new String[30000];
    public static String[] focusBuffer = new String[30000];
    public static int avBufferIndex = 0;
    public static int donBufferIndex = 0;
    public static int accBufferIndex = 0;
    public static int ssBufferIndex = 0;
    public static int focusBufferIndex = 0;

    public static void main (String[] args) {
        if (args.length < 1) {
            System.out.println("must provide vcf/bed file name as command line argument");
            System.exit(1);
        }
        String fileName=args[0];
        writeHeaders(fileName);
        //initialise score tracking variables
        int[] donStart = new int[10000];
        int[] donEnd = new int[10000];
        int[] accStart = new int[10000];
        int[] accEnd = new int[10000];
        String[] donNames = new String[10000];
        String[] accNames = new String[10000];
        double[] scores = new double[12];
        Arrays.fill(scores, -99.0);
        //Note "scores" array order:
        //mesDonRef(0);mesDonAlt(1);mesAccRef(2);mesAccAlt(3);gsDonRef(4);gsDonAlt(5);gsAccRef(6);gsAccAlt(7);ESEmaxRef(8);ESEmaxAlt(9);ESSminRef(10);ESSminAlt(11);
        String s = "";
        String[] geneID = { "", ".", ""}; //(0)chr, (1)name, (2)end
        String[] prevID = { "", "-99", "", "", "GENE"}; //(0)chr, (1)start, (2)ref, (3)alt, (4)type         
        //process lines
        try {
            //File file = new File("/home/steven/teeTest.txt");         
            //BufferedReader br = new BufferedReader(new FileReader(file));
            //while ((s = br.readLine()) != null) {
            BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
            while ((s = in.readLine()) != null && s.length() != 0) {
            String[] split = s.split("\\s+");
            String chr = split[0];
            int start = Integer.parseInt(split[1]);
            String ref = split[2].toUpperCase();
            String alt = split[3].toUpperCase();
            String type = split[4];
            //check if current line is next variant
            boolean nextVariant = false;
            if (!(chr.equals(prevID[0])&&start==Integer.parseInt(prevID[1]))) {
                nextVariant = true;
            } else
            if (!(ref.equals(prevID[2])&&alt.equals(prevID[3]))) {
                nextVariant = true;
            } else
            if (chr.equals("xxx")) {
                nextVariant = true;
            }
            //print scores of previous variant, if current line is next variant and previous line wasn't a gene annotation
            if(nextVariant && !prevID[4].equals("GENE")) {
            	//update overlapping genes and splice sites
            	geneID = updateOverlappingGenes(geneID, prevID);
            	String withinSS = checkForOverlappingSpliceSite(prevID, donNames, accNames, donStart, accStart, donEnd, accEnd);
            	//format output columns
                String[] out = new String[21];
                out[0] = prevID[0]; out[1] = Integer.toString(Integer.parseInt(prevID[1]));
                out[3] = prevID[2]; out[4] = prevID[3]; out[5] = geneID[1]; out[6] = withinSS;
                //adjust end value for indels
                int endPos = Integer.parseInt(prevID[1])+prevID[2].length()-1;
                if (prevID[3].equals("*")) {
                    endPos++;
                }
                out[2]= Integer.toString(endPos);
                //fill scores
                for (int i=0; i<12; i++) {
                    if (scores[i]==-99 || scores[i]==0) {
                        out[i+7]=".";
                    } else {
                        out[i+7]=Double.toString(scores[i]);
                    }
                }
                //calculate donor/acceptor creating logistic regression scores
                String[] lrScores = calculateLogRegScores(scores).split("\\s+");
                if (!out[6].contains("ENSE")) {
                    out[19]= lrScores[0];                
                    out[20]= lrScores[1];
                } else {
                    out[19]= ".";                
                    out[20]= ".";
                }
                //add output line to file buffers
                outputVariantToBuffers(out);
                //reset scores
                Arrays.fill(scores, -99.0);
                //append to file and empty buffers every ~25000 lines
                if (avBufferIndex > 25000) {
                    appendToFiles(fileName);
                    resetOutputArrays();
                }                
            }
        //MaxEntScan
    if (type.equals("MESDON")) {
        scores[0]=Double.parseDouble(split[5]);
        scores[1]=Double.parseDouble(split[6]);
    } else
    if(type.equals("MESACC")){
        scores[2]=Double.parseDouble(split[5]);
        scores[3]=Double.parseDouble(split[6]);
    } else
    //GeneSplicer
    if (type.equals("GSREF")) {
        if (split[6].equals("donor")) {
            scores[4]=Double.parseDouble(split[5]);
        } else
        if (split[6].equals("acceptor")) {
            scores[6]=Double.parseDouble(split[5]);
        } 
    } else   
    if (type.equals("GSALT")) {
        if (split[6].equals("donor")) {
            scores[5]=Double.parseDouble(split[5]);
        } else
        if (split[6].equals("acceptor")) {
            scores[7]=Double.parseDouble(split[5]);
        } 
    } else
    //ESR
        if (type.equals("ESR")) {
            scores[8]=Double.parseDouble(split[5]);
            scores[9]=Double.parseDouble(split[6]);
            scores[10]=Double.parseDouble(split[7]);
            scores[11]=Double.parseDouble(split[8]);
        }
        //GENE
        if (type.equals("GENE")) {
            //overlapping gene annotations
        	if (checkForOverlappingGenes(geneID, chr, start)) {
                geneID[1] = geneID[1].concat(";").concat(split[3]);
                geneID[2] = geneID[2].concat(";").concat(split[2]);
                geneID[0] = split[0];
                if (split.length>11) {
                    //update splice site pos/name arrays with info from new overlapping gene
                    String[] donStartStr=split[6].split(",");
                    String[] donEndStr=split[7].split(",");
                    String[] accStartStr=split[8].split(",");
                    String[] accEndStr=split[9].split(",");
                    String[] donNamesStr=split[10].split(",");
                    String[] accNamesStr=split[11].split(",");
                    int accIndex = 0;
                    int donIndex = 0;
                    while (donStart[donIndex]>0) {
                        donIndex++;
                    }
                    while (accStart[accIndex]>0) {
                        accIndex++;
                    }
                    for (int k = 0; k < donStartStr.length; k++) {
                        donStart[donIndex+k] = Integer.parseInt(donStartStr[k]);
                        donEnd[donIndex+k] = Integer.parseInt(donEndStr[k]);
                        donNames[donIndex+k] = donNamesStr[k];
                    }
                    for (int k = 0; k < accStartStr.length; k++) {
                        accStart[accIndex+k] = Integer.parseInt(accStartStr[k]);
                        accEnd[accIndex+k] = Integer.parseInt(accEndStr[k]);
                        accNames[accIndex+k] = accNamesStr[k];
                    }
                }
        //non-overlapping genes
        } else {
        	//reset arrays
            Arrays.fill(donStart, 0);
            Arrays.fill(donEnd, 0);
            Arrays.fill(accStart, 0);
            Arrays.fill(accEnd, 0);
            Arrays.fill(donNames, null);
            Arrays.fill(accNames, null);
            geneID[1] = split[3];
            geneID[2] = split[2];
            geneID[0] = split[0];
            if (split.length>11) {
                //reset and update splice site pos/name arrays
                String[] donStartStr=split[6].split(",");
                String[] donEndStr=split[7].split(",");
                String[] accStartStr=split[8].split(",");
                String[] accEndStr=split[9].split(",");
                String[] donNamesStr=split[10].split(",");
                String[] accNamesStr=split[11].split(",");
                for (int k=0; k<donStartStr.length; k++) {
                    donStart[k]=Integer.parseInt(donStartStr[k]);
                    donEnd[k]=Integer.parseInt(donEndStr[k]);
                    donNames[k]=donNamesStr[k];
                }
                for (int k=0; k<accStartStr.length; k++) {
                    accStart[k]=Integer.parseInt(accStartStr[k]);
                    accEnd[k]=Integer.parseInt(accEndStr[k]);
                    accNames[k]=accNamesStr[k];
                }
            }        	            	
            }
        }
                //update info for comparison with next line
                prevID[0] = chr;
                prevID[1] = Integer.toString(start);
                prevID[2] = ref;
                prevID[3] = alt;
                prevID[4] = type;
        }
        //final append to file
        if (avBufferIndex > 0) {
            appendToFiles(fileName);
        }
} catch (Exception e) {
	e.printStackTrace();
}
    }
        
    public static String calculateLogRegScores (double[] s) {
        //impute missing values
        //maxEntScan
        for (int i=0; i<4; i++) {
            if (s[i]==-99.0 | s[i]==0.0) {
                s[i] = -20.0;
            }
        }                    
        //gsDonRef
        if (s[4]==-99.0 | s[4]==0.0) {
            if (s[5]==-99.0 | s[5]==0.0) {
                s[4]=0.0;
                s[5]=0.0;
            } else {
                s[4]=-3.0;
            }
        } else
        //gsDonAlt
        if (s[7]==-99.0 | s[7]==0.0) {
            s[7]=-3.0;
        }
        //gsAccRef
        if (s[6]==-99.0 | s[6]==0.0) {
            if (s[7]==-99.0 | s[7]==0.0) {
                s[6]=0.0;
                s[7]=0.0;
            } else {
                s[6]=-3.0;
            }
        } else
        //gsAccAlt
        if (s[7]==-99.0 | s[7]==0.0) {
            s[7]=-3.0;
        }
        double mesDonChange = s[1] - s[0];
        double mesAccChange = s[3] - s[2];
        double gsDonChange = s[5] - s[4];
        double gsAccChange = s[7] - s[6];
        // calculate p = 1/(1+e^-(a + b1X1 + b2X2 + ... + bnXn))
        double pDon = 1/(1 + Math.exp(-(-0.9865 + (mesDonChange * 0.1129) + (gsDonChange * 0.01151) + (s[1] * 0.2076) + (s[5] * 0.4350) )));
        double pAcc = 1/(1 + Math.exp(-(-1.665 + (mesAccChange * 0.3323) + (gsAccChange * 0.05084) + (s[3] * 0.1877) + (s[7] * 0.1730) )));
        DecimalFormat df = new DecimalFormat("0.00");
        pDon = Double.valueOf(df.format(pDon));
        pAcc = Double.valueOf(df.format(pAcc));
        String ret = Double.toString(pDon) + "\t" + Double.toString(pAcc);
        return ret;
    }

public static void appendToFiles (String fileName) {
    String annovarName = "output/"+fileName+"_out.txt";
    String accName = "temp/"+fileName+"_acceptorCreating_unsorted.txt";
    String donName = "temp/"+fileName+"_donorCreating_unsorted.txt";
    String withinSSname = "temp/"+fileName+"_withinSS_unsorted.txt";
    String focussedName = "output/"+fileName+"_focussed.txt";
    try {
        //write to _out
        FileWriter fwAV = new FileWriter(annovarName, true);
        BufferedWriter avWriter = new BufferedWriter(fwAV);
        for (int i=0; i<avBufferIndex; i++) {
            avWriter.write(avBuffer[i]+"\n");
        }
        avWriter.close();
        //write to _focussed
        FileWriter fwFC = new FileWriter(focussedName, true);
        BufferedWriter fcWriter = new BufferedWriter(fwFC);
        for (int i=0; i<focusBufferIndex; i++) {
            fcWriter.write(focusBuffer[i]+"\n");
        }
        fcWriter.close();
    } catch (IOException e) {
        System.out.println(e.getMessage());
    }
}

public static void writeHeaders (String fileName) {
    String annovarName = "output/"+fileName+"_out.txt";
    String accName = "temp/"+fileName+"_acceptorCreating_unsorted.txt";
    String donName = "temp/"+fileName+"_donorCreating_unsorted.txt";
    String withinSSname = "temp/"+fileName+"_withinSS_unsorted.txt";
    String focussedName = "output/"+fileName+"_focussed.txt";
    try {
        //write to _out
        FileWriter fwAV = new FileWriter(annovarName);
        BufferedWriter avWriter = new BufferedWriter(fwAV);
        avWriter.write("#CHR\tSTART\tEND\tREF\tALT\tGENE\twithinSite\tmesDonRef\tmesDonAlt\tmesAccRef\tmesAccAlt\tgsDonRef\tgsDonAlt\tgsAccRef\tgsAccAlt\tESEmaxRef\tESEmaxAlt\tESSminRef\tESSminAlt\tdonCreateP\taccCreateP"+"\n");
        avWriter.close();
        //write to _focussed
        FileWriter fwFC = new FileWriter(focussedName);
        BufferedWriter fcWriter = new BufferedWriter(fwFC);
        fcWriter.write("#CHR\tSTART\tEND\tREF\tALT\tGENE\twithinSite\tmesDonRef\tmesDonAlt\tmesAccRef\tmesAccAlt\tgsDonRef\tgsDonAlt\tgsAccRef\tgsAccAlt\tdonCreateP\taccCreateP\n");
        fcWriter.close();
    } catch (IOException e) {
        System.out.println(e.getMessage());
    }
} 

public static void resetOutputArrays() {
    avBuffer = new String[30000];
    donBuffer = new String[30000];
    accBuffer = new String[30000];
    ssBuffer = new String[30000];
    focusBuffer = new String[30000];
    avBufferIndex = 0; donBufferIndex = 0; accBufferIndex = 0; ssBufferIndex = 0; focusBufferIndex = 0;
    System.gc();
}

//String[] geneID = { "", ".", ""}; //(0)chr, (1)name, (2)end
//String[] prevID = { "", "-99", "", "", "GENE"}; //(0)chr, (1)start, (2)ref, (3)alt, (4)type 

public static String[] updateOverlappingGenes(String[] geneID, String[] prevID) {
    //check if previous variant falls within an annotated splice site motif
	//if previous genes overlapped, check individually whether they overlap this variant
    if (geneID[1].contains(";")) {
        String[] geneNameSplit = geneID[1].split(";");
        String[] geneEndSplit = geneID[2].split(";");
        geneID[1]="";
        geneID[2]="";
        int geneIndex = Math.min(geneEndSplit.length, geneNameSplit.length);
        for (int i=0; i<geneIndex; i++) {
            if (Integer.parseInt(geneEndSplit[i]) >= Integer.parseInt(prevID[1])) {
                if (!geneID[1].equals("")) {
                    geneID[1] = geneID[1].concat(";");
                    geneID[2] = geneID[2].concat(";");
                }
                geneID[1] = geneID[1].concat(geneNameSplit[i]);
                geneID[2] = geneID[2].concat(geneEndSplit[i]);
            }
        }
    }
    if ( geneID[1].equals("") || !geneID[0].equals(prevID[0])) {
        geneID[1]= ".";
        geneID[2]= ".";
    }
    return geneID;
}

	public static String checkForOverlappingSpliceSite( String[] prevID, String[] donNames, String[] accNames, int[] donStart, int[] accStart, int[] donEnd, int[] accEnd) {
    String withinSS="";
    int startPos = Integer.parseInt(prevID[1]);
    int endPos = Integer.parseInt(prevID[1])+prevID[2].length()-1;
    if (prevID[3].equals("*")) {
        endPos++;
    }
    for (int i=0; i<donNames.length; i++) {
        if (donNames[i]==null) {
            break;
        }
        if (!withinSS.contains(donNames[i])) {
            if (endPos>=donStart[i]&&startPos<=donEnd[i]) {
                if (!withinSS.equals("")) {
                        withinSS=withinSS.concat(",");
                }
            withinSS = withinSS.concat(donNames[i]).concat("_donor");
            }
        }
    }
    for (int i=0; i<accNames.length; i++) {
        if (accNames[i]==null) {
            break;
        }
        if (!withinSS.contains(accNames[i])) {
            if (endPos>=accStart[i]&&startPos<=accEnd[i]) {
                if (!withinSS.equals("")) {
                    withinSS=withinSS.concat(",");
                }
                withinSS = withinSS.concat(accNames[i]).concat("_acceptor");
            }
        }
    }
    if (withinSS.equals("")) {
        withinSS=".";
    }
    return withinSS;
	}
	
	public static void outputVariantToBuffers(String[] out) {	
    //write full line for annovar
    String avLine = "";
    for (int i=0; i<20; i++) {
        avLine = avLine+out[i]+"\t";
    }
    avLine = avLine+out[20];
    avBuffer[avBufferIndex] = avLine;
    avBufferIndex++;
    //write withinSS
    if (out[6].contains("ENSE")) {
        String ssLine = "";
        for (int i=0; i<15; i++) {
            ssLine = ssLine+out[i]+"\t";
        }
        ssLine = ssLine+out[19]+"\t"+out[20];
        ssBuffer[ssBufferIndex] = ssLine;
        ssBufferIndex++;
        focusBuffer[focusBufferIndex] = ssLine;
        focusBufferIndex++;
    } else {
        //write donCreating
        if (Double.parseDouble(out[19])>=0.7) {
            String donLine = "";
            for (int i=0; i<15; i++) {
                donLine= donLine+out[i]+"\t";
            }
            donLine = donLine+out[19]+"\t"+out[20];
            donBuffer[donBufferIndex] = donLine;
            donBufferIndex++;
            focusBuffer[focusBufferIndex] = donLine;
            focusBufferIndex++;
        }
        //write accCreating
        if (Double.parseDouble(out[20])>=0.7) {
            String accLine = "";
            for (int i=0; i<15; i++) {
                accLine= accLine+out[i]+"\t";
            }
            accLine = accLine+out[19]+"\t"+out[20];
            accBuffer[accBufferIndex] = accLine;
            accBufferIndex++;
            focusBuffer[focusBufferIndex] = accLine;
            focusBufferIndex++;
        }
    }              
}
	
	public static boolean checkForOverlappingGenes(String[] geneID, String chr, int start) {
        //check for overlapping genes
        int overlapTest = 0;
        //check individually for already overlapping genes
        if (geneID[1].concat(geneID[2]).contains(";")) {
            String[] geneEndSplit = geneID[2].split(";");
            String[] geneNameSplit = geneID[1].split(";");
            String updatedGeneName = "";
            String updatedGeneEnd = "";
            int geneIndex = Math.min(geneEndSplit.length, geneNameSplit.length);
            for (int i=0; i<geneIndex; i++) {
                if (Integer.parseInt(geneEndSplit[i]) >= start) {
                    overlapTest = Integer.parseInt(geneEndSplit[i]);
                    if (!updatedGeneName.equals("")) {
                        updatedGeneName = updatedGeneName.concat(";");
                        updatedGeneEnd = updatedGeneEnd.concat(";");
                    }
                    updatedGeneName = updatedGeneName.concat(geneNameSplit[i]);
                    updatedGeneEnd = updatedGeneEnd.concat(geneEndSplit[i]);
                }
            }
            geneID[1] = updatedGeneName;
            geneID[2] = updatedGeneEnd;
        } else if (geneID[2].equals(".")) {    
            overlapTest = 0;
            geneID[2]="0";
        } else if (!geneID[2].equals("")) {    
            overlapTest = Integer.parseInt(geneID[2]);
        }
        if (start <= overlapTest && geneID[0].equals(chr)) {
        	return true;
        }
        return false;
	}
	}



