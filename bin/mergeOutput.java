//Description: Merging GeneSplicer, MaxEntScan and ESRseq scores along with annotation information to be outputted together on one line per variant.
import java.io.*;
import java.util.Arrays;
import java.util.jar.Attributes.Name;
import java.lang.*;
public class mergeOutput {
		
			public static void main (String[] args) {
			    try {
			    BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
				System.out.println("#CHR\tSTART\tEND\tREF\tALT\tGENE\twithinSite\tmesDonRef\tmesDonAlt\tmesAccRef\tmesAccAlt\tgsDonRef\tgsDonAlt\tgsAccRef\tgsAccAlt\tESEmaxRef\tESEmaxAlt\tESSminRef\tESSminAlt");
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
		                    if (endPos>=donStart[i]&&startPos<=donEnd[i]) {
		                        if (!withinSS.equals("")) {
		                            withinSS=withinSS.concat(",");
		                        }
		                        withinSS = withinSS.concat(donNames[i]).concat("_donor");
		                    }
		                }
		                for (int i=0; i<accNames.length; i++) {
		                	if (accNames[i]==null) {
		                		break;
		                	}
		                    if (endPos>=accStart[i]&&startPos<=accEnd[i]) {
		                        if (!withinSS.equals("")) {
		                            withinSS=withinSS.concat(",");
		                        }
		                        withinSS = withinSS.concat(accNames[i]).concat("_acceptor");
		                    }
		                }
		                if (withinSS.equals("")) {
		                    withinSS=".";
		                }
		                String[] id = { prevChr, Integer.toString(prevStart), Integer.toString(endPos), prevRef, prevAlt, geneName, withinSS };
		                printScores(id, scores);
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
			            } else if (!geneEnd.equals(".")) {    
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
			        for (int i=0; i<11; i++) {
			            System.out.print(scores[i]+"\t");
			        }
			        System.out.println(scores[11]);        
			    }
			}
