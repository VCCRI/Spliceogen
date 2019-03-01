//Description: Merging GeneSplicer, MaxEntScan and ESRseq scores along with annotation information to be outputted together on one line per variant.
import java.io.*;
import java.util.Arrays;
import java.util.jar.Attributes.Name;
import java.lang.*;
import java.text.DecimalFormat;
public class mergeOutput {
		
			public static void main (String[] args) {
                if (args.length < 1) {
                    System.out.println("must provide input file as command line argument");
                    System.exit(1);
                }
                String fileName=args[0];
                printHeaders(fileName);
			    try {
			    BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
			    //initialise score tracking variables
			    String prevChr = "";
			    int prevStart = -99;
			    String prevRef = "";
			    String prevAlt = "";
			    double[] scores = new double[12];
			    int[] donStart = new int[10000];
			    int[] donEnd = new int[10000];
			    int[] accStart = new int[10000];
			    int[] accEnd = new int[10000];
			    String[] donNames = new String[10000];
			    String[] accNames = new String[10000];
			    Arrays.fill(scores, -99.0);
			    //Note "scores" array order:
			    //mesDonRef(0);mesDonAlt(1);mesAccRef(2);mesAccAlt(3);gsDonRef(4);gsDonAlt(5);gsAccRef(6);gsAccAlt(7);ESEmaxRef(8);ESEmaxAlt(9);ESSminRef(10);ESSminAlt(11);
			    String s = "";
			    String geneName = ".";
			    String geneChr = "";
			    String prevType = "GENE";
			    String geneEnd = "";
			    //process lines
			    while ((s = in.readLine()) != null && s.length() != 0) {
			        String[] split = s.split("\\s+");
			        try {
			        String chr = split[0];
			        int start = Integer.parseInt(split[1]);
			        String ref = split[2].toUpperCase();
			        String alt = split[3].toUpperCase();
			        String type = split[4];
			        //check if current line is next variant
			        boolean nextVariant = false;
			        if (!(chr.equals(prevChr)&&start==prevStart)) {
			            nextVariant = true;
			        } else
			        if (!(ref.equals(prevRef)&&alt.equals(prevAlt))) {
			            nextVariant = true;
			        } else
			        if (chr.equals("xxx")) {
			            nextVariant = true;
			        }
			        //print scores of previous variant, if current line is next variant and previous line wasn't a gene annotation
			        if(nextVariant && !prevType.equals("GENE")) {
		                //check if previous variant falls within an annotated splice site motif
			        	if (geneName.contains(";")) {
			        		String[] geneNameSplit = geneName.split(";");
			        		String[] geneEndSplit = geneEnd.split(";");
			        		geneName="";
		                    geneEnd="";
			                int geneIndex = Math.min(geneEndSplit.length, geneNameSplit.length);
			        		for (int i=0; i<geneIndex; i++) {
			        			if (Integer.parseInt(geneEndSplit[i]) >= prevStart) {
			        				if (!geneName.equals("")) {
			        					geneName = geneName.concat(";");
			        					geneEnd = geneEnd.concat(";");
			        				}
			        				geneName = geneName.concat(geneNameSplit[i]);
			        				geneEnd = geneEnd.concat(geneEndSplit[i]);
			        			}
			        		}
		                }
		                if ( geneName.equals("") || !geneChr.equals(prevChr)) {
		                    geneName= ".";
		                    geneEnd= ".";
		                }
		                String withinSS="";
		                int startPos = prevStart;
		                int endPos = prevStart+prevRef.length()-1;
		                if (prevAlt.equals("*")) {
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
		                String[] id = { prevChr, Integer.toString(prevStart), Integer.toString(endPos), prevRef, prevAlt, geneName, withinSS };
		                printScores(id, scores, fileName);
		                //Reset scores
			            Arrays.fill(scores, -99.0);
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
			            //check for overlapping genes
			            int overlapTest = 0;
			            //check individually for already overlapping genes
			            if (geneName.concat(geneEnd).contains(";")) {
			                String[] geneEndSplit = geneEnd.split(";");
			                String[] geneNameSplit = geneName.split(";");
			                String updatedGeneName = "";
			                String updatedGeneEnd = "";
			                int geneIndex = Math.min(geneEndSplit.length, geneNameSplit.length);
			                for (int i=0; i<geneIndex; i++) {
			                    if (Integer.parseInt(geneEndSplit[i]) >= Integer.parseInt(split[1])) {
			                    	overlapTest = Integer.parseInt(geneEndSplit[i]);
			                    	if (!updatedGeneName.equals("")) {
			                    		updatedGeneName = updatedGeneName.concat(";");
			                    		updatedGeneEnd = updatedGeneEnd.concat(";");
			                    	}
			                        updatedGeneName = updatedGeneName.concat(geneNameSplit[i]);
			                        updatedGeneEnd = updatedGeneEnd.concat(geneEndSplit[i]);
			                    }
			                }
			                geneName = updatedGeneName;
			                geneEnd = updatedGeneEnd;
			            } else if (geneEnd.equals(".")) {    
			                overlapTest = 0;
			                geneEnd="0";
			            } else if (!geneEnd.equals("")) {    
			                overlapTest = Integer.parseInt(geneEnd);
			            }
			            // if overlapping genes           
			            if (Integer.parseInt(split[1]) <= overlapTest && geneChr.equals(split[0])) {
			                //handle overlapping gene annotations
			                geneName = geneName.concat(";").concat(split[3]);
			                geneEnd = geneEnd.concat(";").concat(split[2]);
			                geneChr = split[0];
			                if (split.length>6) {
			                    //update splice site pos/name arrays with info from new overlapping gene
			                    String[] donNamesStr=split[10].split(",");
			                    String[] accNamesStr=split[11].split(",");
			                    String[] donStartStr=split[6].split(",");
			                    String[] donEndStr=split[7].split(",");
			                    String[] accStartStr=split[8].split(",");
			                    String[] accEndStr=split[9].split(",");
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
			            	geneName = split[3];
			            	geneEnd = split[2];
			            	geneChr = split[0];
			            	if (split.length>6) {
			            		//reset and update splice site pos/name arrays
			            		Arrays.fill(donStart, 0);
			            		Arrays.fill(donEnd, 0);
			            		Arrays.fill(accStart, 0);
			            		Arrays.fill(accEnd, 0);
			            		Arrays.fill(donNames, null);
			            		Arrays.fill(accNames, null);
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
			        prevChr = chr;
			        prevStart = start;
			        prevRef = ref;
			        prevAlt = alt;
			        prevType = type;
			        //report line ID
			        } catch (Exception e) {
			        	e.printStackTrace();
			        	System.out.println("Merging error in line: " +s);
			        }
			    }
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
/*
    public static void printScores(String[] id , double[] s) {
        String[] scores = new String[12];
        for (int i=0; i<12; i++) {
            scores[i]=Double.toString(s[i]);
            if (scores[i].equals("-99.0") || scores[i].equals("0.0")) {
                scores[i]=".";
            }
        }
        for (int i=0; i<7; i++) {
            System.out.print(id[i]+"\t");   
        }
        for (int i=0; i<12; i++) {
            System.out.print(scores[i]+"\t");
        }
        System.out.println(calculateLogRegScores(s));        
    }
*/
    public static void printScores(String[] id , double[] s, String fileName) {
        String line = "";
        String[] out = new String[21];
        for (int i=0; i<7; i++) {
            out[i]=id[i];   
        }
        for (int i=0; i<12; i++) {
            if (s[i]==-99 || s[i]==0) {
                out[i+7]=".";
            } else {
                out[i+7]=Double.toString(s[i]);
            }
        }
        String[] lrScores = calculateLogRegScores(s).split("\\s+");
        out[19]= lrScores[0];                
        out[20]= lrScores[1];                
    //append to files
    String annovarName = "output/"+fileName+"_out.txt";
    String accName = "temp/"+fileName+"_acceptorCreating_unsorted.txt";
    String donName = "temp/"+fileName+"_donorCreating_unsorted.txt";
    String withinSSname = "temp/"+fileName+"_withinSS_unsorted.txt";
    try {
        //write full line for annovar
        FileWriter fwAV = new FileWriter(annovarName, true);
        BufferedWriter avWriter = new BufferedWriter(fwAV);
        for (int i=0; i<out.length-1; i++) {
            avWriter.write(out[i]+"\t");
        }
        avWriter.write(out[out.length-1]+"\n");
        avWriter.close();
        //write withinSS
        if (out[6].contains("ENSE")) {
            FileWriter fwSS = new FileWriter(withinSSname, true);
            BufferedWriter ssWriter = new BufferedWriter(fwSS);
            for (int i=0; i<16; i++) {
                ssWriter.write(out[i]+"\t");
            }
            double maxScoreDecrease = -99;
            //if donor
            if (out[6].contains("_donor")) {
                if (!out[7].equals(".") && !out[8].equals(".")) {
                    maxScoreDecrease = s[0] - s[1]; 
                }
            }
            //if acceptor
            if (out[6].contains("_acceptor")) {
                if (!out[9].equals(".") && !out[10].equals(".")) {
                    if (maxScoreDecrease < s[2] - s[3] ) {
                        maxScoreDecrease = s[2] - s[3]; 
                    }
                }
            }
            ssWriter.write(maxScoreDecrease+"\n");
            ssWriter.close();
        } else {
            //write donCreating
            if (Double.parseDouble(out[19])>=0.7) {
                FileWriter fwDon = new FileWriter(donName, true);
                BufferedWriter donWriter = new BufferedWriter(fwDon);
                for (int i=0; i<6; i++) {
                    donWriter.write(out[i]+"\t");
                }
                donWriter.write(out[7]+"\t"+out[8]+"\t"+out[11]+"\t"+out[12]+"\t"+out[19]+"\n");
                donWriter.close();
            }
            //write accCreating
            if (Double.parseDouble(out[20])>=0.7) {
                FileWriter fwAcc = new FileWriter(accName, true);
                BufferedWriter accWriter = new BufferedWriter(fwAcc);
                for (int i=0; i<6; i++) {
                    accWriter.write(out[i]+"\t");
                }
                accWriter.write(out[9]+"\t"+out[10]+"\t"+out[13]+"\t"+out[14]+"\t"+out[20]+"\n");
                accWriter.close();
            }
        }
    } catch (IOException e) {
        System.out.println(e.getMessage());
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
			        //double pDon = 1/(1 + Math.exp(-(-0.8440 + (mesDonChange * 0.09236) + (gsDonChange * -0.05291) + (s[1] * 0.2053) + (s[5] * 0.4548) )));
			        //double pAcc = 1/(1 + Math.exp(-(-1.615 + (mesAccChange * 0.3256) + (gsAccChange * 0.09301) + (s[3] * 0.1879) + (s[7] * 0.1311) )));
                    DecimalFormat df = new DecimalFormat("0.00");
                    pDon = Double.valueOf(df.format(pDon));
                    pAcc = Double.valueOf(df.format(pAcc));
                    String ret = Double.toString(pDon) + "\t" + Double.toString(pAcc);
                    return ret;
                }


    public static void printHeaders(String fileName) {
    String annovarName = "output/"+fileName+"_out.txt";
    String accName = "output/"+fileName+"_acceptorCreating.txt";
    String donName = "output/"+fileName+"_donorCreating.txt";
    String withinSSname = "output/"+fileName+"_withinSS.txt";
    try {
        //annovar
        FileWriter fwAV = new FileWriter(annovarName);
        BufferedWriter avWriter = new BufferedWriter(fwAV);
        avWriter.write("#CHR\tSTART\tEND\tREF\tALT\tGENE\twithinSite\tmesDonRef\tmesDonAlt\tmesAccRef\tmesAccAlt\tgsDonRef\tgsDonAlt\tgsAccRef\tgsAccAlt\tESEmaxRef\tESEmaxAlt\tESSminRef\tESSminAlt\tdonCreateP\taccCreateP"+"\n");
        avWriter.close();
        //withinSS
        FileWriter fwSS = new FileWriter(withinSSname);
        BufferedWriter ssWriter = new BufferedWriter(fwSS);
        ssWriter.write("#CHR\tSTART\tEND\tREF\tALT\tGENE\twithinSite\tmesDonRef\tmesDonAlt\tmesAccRef\tmesAccAlt\tgsDonRef\tgsDonAlt\tgsAccRef\tgsAccAlt"+"\n");
        ssWriter.close();
        //donCreating
        FileWriter fwDon = new FileWriter(donName);
        BufferedWriter donWriter = new BufferedWriter(fwDon);
        donWriter.write("#CHR\tSTART\tEND\tREF\tALT\tGENE\tmesDonRef\tmesDonAlt\tgsDonRef\tgsDonAlt\tdonCreateP"+"\n");
        donWriter.close();
        //accCreating
        FileWriter fwAcc = new FileWriter(accName);
        BufferedWriter accWriter = new BufferedWriter(fwAcc);
        accWriter.write("#CHR\tSTART\tEND\tREF\tALT\tGENE\tmesAccRef\tmesAccAlt\tgsAccRef\tgsAccAlt\taccCreateP"+"\n");
        accWriter.close();
    } catch (IOException e) {
        System.out.println(e.getMessage());
    }
}    

}
