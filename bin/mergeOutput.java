//Description: Merging GeneSplicer, MaxEntScan and ESRseq scores along with annotation information to be outputted together on one line per variant.
import java.io.*;
import java.util.Arrays;
import java.util.jar.Attributes.Name;
import java.lang.*;
import java.text.DecimalFormat;
public class mergeOutput {

    //store output lines in buffers to minimise I/O
    public static String[] avBuffer = new String[30000];
    public static String[] gainBuffer = new String[30000];
    public static String[] ssBuffer = new String[30000];
    public static int gainBufferIndex = 0;
    public static int avBufferIndex = 0;
    public static int ssBufferIndex = 0;

    public static void main (String[] args) {
        if (args.length < 3) {
            System.out.println("Too few input args...exiting");
            System.exit(1);
        }
        String fileName=args[0];
        String chrAdd=args[1];
        String chrRemove=args[2];
        writeHeaders(fileName);
        //initialise score tracking variables
        int[] donStart = new int[10000];
        int[] donEnd = new int[10000];
        int[] accStart = new int[10000];
        int[] accEnd = new int[10000];
        String[] donNames = new String[10000];
        String[] accNames = new String[10000];
        //Note array elements:
        //scores[]...mesDonRef(0);mesDonAlt(1);mesAccRef(2);mesAccAlt(3);gsDonRef(4);gsDonAlt(5);gsAccRef(6);gsAccAlt(7);ESEmaxRef(8);ESEmaxAlt(9);ESSminRef(10);ESSminAlt(11);
        //scores[]...mesAccRef2ndScore(12);mesAccAlt2ndScore(13);mesAccRefPos1(14);mesAccAltPos1(15);mesAccRefPos2(16);mesAccAltPos2(17);
        //geneID[]...chr(0), name(1), end(2)
        //prevID[]...chr(0), start(1), ref(2), alt(3), type(4)
        double[] scores = new double[18];
        Arrays.fill(scores, -99.0);
        String[] geneID = { "", ".", ""};
        String[] prevID = { "", "-99", "", "", "GENE"};
        String s = "";
        //process lines
        try {
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
                    String[] out = new String[23];
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
                    String[] lrScores = calculateLogRegScores(scores, out).split("\\s+");
                    out[19]= lrScores[0];
                    out[20]= lrScores[1];
                    out[21]= lrScores[2];
                    out[22]= lrScores[3];
                    if (chrAdd.equals("inputAdd=chr")) {
                        out[0]="chr"+out[0];
                    } else if (chrRemove.equals("inputRemove=chr")) {
                        out[0]=out[0].substring(3);
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
                        scores[12]=Double.parseDouble(split[7]);
                        scores[13]=Double.parseDouble(split[8]);
                        scores[14]=Double.parseDouble(split[9]);
                        scores[15]=Double.parseDouble(split[10]);
                        scores[16]=Double.parseDouble(split[11]);
                        scores[17]=Double.parseDouble(split[12]);
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
                    if (ifGenesOverlap(geneID, chr, start)) {
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

    public static String calculateLogRegScores (double[] s, String[] out) {
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
        double pDonGain = -1; double pAccGain = -1; double pDonLoss = -1; double pAccLoss = -1;
        // using within/outside splice site logistic regression models, calculate p = 1/(1+e^-(a + b1X1 + b2X2 + ... + bnXn))
        if (!out[6].contains("ENSE")) {
            pDonGain = 1/(1 + Math.exp(-(-0.9865 + (mesDonChange * 0.1129) + (gsDonChange * 0.01151) + (s[1] * 0.2076) + (s[5] * 0.4350) )));
            pAccGain = 1/(1 + Math.exp(-(-1.665 + (mesAccChange * 0.3323) + (gsAccChange * 0.05084) + (s[3] * 0.1877) + (s[7] * 0.1730) )));
        }
        if (out[6].contains("donor")){
            //pDonLoss = 1/(1 + Math.exp(-(-1.6678 + (mesDonChange * -0.6752) )));
            pDonLoss = 1/(1 + Math.exp(-(-1.1028 + (mesDonChange * -0.4240) )));
            double denovoScore = s[13];
            double ogScore = s[3];
            //if highest alt score is at a different site (accounting for indels)
            if (Math.abs(out[3].length() - out[4].length()) < Math.abs(s[14] - s[15]) ) {
                denovoScore = s[3];
                ogScore = s[13];
            }
            double denovoChange = denovoScore - s[12];
            pDonGain = 1/(1 + Math.exp(-(99.705 + (s[2] * -0.36584) + (s[12] * 12.3598) + (denovoScore * -5.7997) + (ogScore * 0.41298) )));
        }
        if (out[6].contains("acceptor")){
            //pAccLoss = 1/(1 + Math.exp(-(-0.6209 + (mesAccChange * -0.5100) )));
            pAccLoss = 1/(1 + Math.exp(-(-0.2355 + (mesAccChange * -0.2283) )));
            double denovoScore = s[13];
            double ogScore = s[3];
            //if highest alt score is at a different site (accounting for indels)
            if (Math.abs(out[3].length() - out[4].length()) < Math.abs(s[14] - s[15]) ) {
                denovoScore = s[3];
                ogScore = s[13];
            }
            double denovoChange = denovoScore - s[12];
            //pAccGain = 1/(1 + Math.exp(-(-0.18441 + (s[2] * -0.13376) + (s[12] * -0.54157) + (denovoScore * 0.61567) + (ogScore * -0.069948) )));
            pAccGain = 1/(1 + Math.exp(-(1.2445 + (s[2] * 0.55004) + (s[12] * -0.06490) + (denovoScore * 0.33992) + (ogScore * -0.000419) )));
        }
        String pDonGainStr = Double.toString(pDonGain);
        String pAccGainStr = Double.toString(pAccGain);
        String pDonLossStr = Double.toString(pDonLoss);
        String pAccLossStr = Double.toString(pAccLoss);
        //round scores to 2 decimal places
        try {
            pDonGainStr = String.format("%.02f",pDonGain);
            pAccGainStr = String.format("%.02f",pAccGain);
            pDonLossStr = String.format("%.02f",pDonLoss);
            pAccLossStr = String.format("%.02f",pAccLoss);
        } catch (Exception e) {
            //use unrounded score
        }
        if (pDonGainStr.equals("-1.00")) {
            pDonGainStr = ".";
        }
        if (pAccGainStr.equals("-1.00")) {
            pAccGainStr = ".";
        }
        if (pDonLossStr.equals("-1.00")) {
            pDonLossStr = ".";
        }
        if (pAccLossStr.equals("-1.00")) {
            pAccLossStr = ".";
        }
        String ret = pDonGainStr + "\t" + pAccGainStr + "\t" + pDonLossStr + "\t" + pAccLossStr;
        return ret;
    }

    public static void appendToFiles (String fileName) {
        String annovar_file = "output/"+fileName+"_out.txt";
        String gain_file_unsorted = "temp/"+fileName+"_gain_unsorted.txt";
        String ss_file_unsorted = "temp/"+fileName+"_loss_unsorted.txt";
        try {
            //write to _out
            FileWriter fwAV = new FileWriter(annovar_file, true);
            BufferedWriter avWriter = new BufferedWriter(fwAV);
            for (int i=0; i<avBufferIndex; i++) {
                avWriter.write(avBuffer[i]+"\n");
            }
            avWriter.close();
            //write to _acceptorCreating
            FileWriter fwGain = new FileWriter(gain_file_unsorted, true);
            BufferedWriter gainWriter = new BufferedWriter(fwGain);
            for (int i=0; i<gainBufferIndex; i++) {
                gainWriter.write(gainBuffer[i]+"\n");
            }
            gainWriter.close();
            //write to _withinSS
            FileWriter fwSS = new FileWriter(ss_file_unsorted, true);
            BufferedWriter ssWriter = new BufferedWriter(fwSS);
            for (int i=0; i<ssBufferIndex; i++) {
                ssWriter.write(ssBuffer[i]+"\n");
            }
            ssWriter.close();
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }
    }

    public static void writeHeaders (String fileName) {
        String annovar_file = "output/"+fileName+"_out.txt";
        String gain_file = "temp/"+fileName+"_gain_unsorted.txt";
        String ss_file = "temp/"+fileName+"_loss_unsorted.txt";
        try {
            //write to _out
            FileWriter fwAV = new FileWriter(annovar_file);
            BufferedWriter avWriter = new BufferedWriter(fwAV);
            avWriter.write("#CHR\tSTART\tEND\tREF\tALT\tGENE\twithinSite\tmesDonRef\tmesDonAlt\tmesAccRef\tmesAccAlt\tgsDonRef\tgsDonAlt\tgsAccRef\tgsAccAlt\tESEmaxRef\tESEmaxAlt\tESSminRef\tESSminAlt\tdonGainP\taccGainP\tdonLossP\taccLossP"+"\n");
            avWriter.close();
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }
    } 

    public static void resetOutputArrays() {
        avBuffer = new String[30000];
        gainBuffer = new String[30000];
        ssBuffer = new String[30000];
        avBufferIndex = 0; gainBufferIndex = 0; ssBufferIndex = 0;
        System.gc();
    }

    public static String[] updateOverlappingGenes(String[] geneID, String[] prevID) {
        //check if previous variant falls within an annotated splice site motif
        //if previous genes overlapped, check individually whether they overlap this variant and update
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
        for (int i=0; i<22; i++) {
            avLine = avLine+out[i]+"\t";
        }
        avLine = avLine+out[22];
        avBuffer[avBufferIndex] = avLine;
        avBufferIndex++;
        //write withinSS
        if (out[6].contains("ENSE")) {
            String ssLine = out[0]+"\t"+out[1]+"\t"+out[3]+"\t"+out[4]+"\t"+out[5]+"\t"+out[6]+"\t"+out[21]+"\t"+out[22];
            Double lossMax = 0.00 ;
            try {
                if (!out[21].equals(".")) {
                    lossMax = Double.parseDouble(out[21]);
                }
                if (!out[22].equals(".")) {
                    if (lossMax < Double.parseDouble(out[22])) {
                        lossMax = Double.parseDouble(out[22]);
                    }
                }
            } catch (Exception e) {}
            ssLine = ssLine+"\t"+lossMax;
            ssBuffer[ssBufferIndex] = ssLine;
            ssBufferIndex++;
        } else {
            //write gain
            String gainLine = out[0]+"\t"+out[1]+"\t"+out[3]+"\t"+out[4]+"\t"+out[5]+"\t"+out[19]+"\t"+out[20];
            Double gainMax = 0.00 ;
            try {
                if (!out[19].equals(".")) {
                    gainMax = Double.parseDouble(out[19]);
                }
                if (!out[20].equals(".")) {
                    if (gainMax < Double.parseDouble(out[20])) {
                        gainMax = Double.parseDouble(out[20]);
                    }
                }
            } catch (Exception e) {}
            if (gainMax>=0.7) {
                gainLine = gainLine+"\t"+Double.toString(gainMax);
                gainBuffer[gainBufferIndex] = gainLine;
                gainBufferIndex++;
            }
        }              
    }

    public static boolean ifGenesOverlap(String[] geneID, String chr, int start) {
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
