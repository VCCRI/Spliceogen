//Description: Merging GeneSplicer, MaxEntScan and ESRseq scores along with annotation information to be outputted together on one line per variant.
import java.io.*;
import java.util.Arrays;
class mergeOutput {
public static void main (String[] args) {
    try {
    BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
    //print header
    System.out.println("#CHR\t#START\t#REF\t#ALT\t#GENE\t#withinSite\t#mesDonRef\t#mesDonAlt\t#mesAccRef\t#mesAccAlt\t#gsDonRef\t#gsDonAlt\t#gsAccRef\t#gsAccAlt\t#ESEmaxRef\t#ESEmaxAlt\t#ESSminRef\t#ESSminAlt");
    //initialise score tracking variables
    String prevChr = "";
    int prevStart = -99;
    String prevRef = "";
    String prevAlt = "";
    double[] scores = new double[12];
    String[] finalLineID = new String[6];
    double[] finalLineScores = new double[scores.length];
    int[] donStart = new int[1000];
    int[] donEnd = new int[1000];
    int[] accStart = new int[1000];
    int[] accEnd = new int[1000];
    String[] donNames = new String[1000];
    String[] accNames = new String[1000];
    Arrays.fill(scores, -99.0);
    Arrays.fill(finalLineScores, -99.0);
    //Note "scores" array order:
    //mesDonRef(0);mesDonAlt(1);mesAccRef(2);mesAccAlt(3);gsDonRef(4);gsDonAlt(5);gsAccRef(6);gsAccAlt(7);ESEmaxRef(8);ESEmaxAlt(9);ESSminRef(10);ESSminAlt(11);
    String s = "";
    String geneName = "";
    String geneChr = "";
    String prevType = "";
    int geneEnd = 0;
    //process lines
    while ((s = in.readLine()) != null && s.length() != 0) {
        String[] split = s.split("\\s+");
        String chr = split[0];
        int start = Integer.parseInt(split[1]);
        String ref = split[2];
        String alt = split[3];
        String type = split[4];
        //check if current line is next variant
        boolean nextVariant = false;
        if (!(chr.equals(prevChr)&&start==prevStart)) {
            nextVariant = true;
        } else
        if (!(ref.equals(prevRef)&&alt.equals(prevAlt))) {
            nextVariant = true;
        }
        //print scores of previous variant, if current line is next variant and previous line wasn't a gene annotation
        if(nextVariant && !prevType.equals("GENE")) {
            if (prevStart!=-99) {
                if (!(prevStart <= geneEnd && geneChr.equals(prevChr))) {
                    geneName= ".";
                }
                //check if previous variant falls within an annotated splice site motif
                String withinSS="";
                int startPos = prevStart;
                int endPos = prevStart+prevRef.length()-1;
                if (prevAlt.equals("*")) {
                    endPos++;
                }
                for (int i=0; i<donNames.length; i++) {
                    if (endPos>=donStart[i]&&startPos<=donEnd[i]) {
                        if (!withinSS.equals("")) {
                            withinSS=withinSS.concat(",");
                        }
                        withinSS = withinSS.concat(donNames[i]).concat("_donor");
                    }
                }
                for (int i=0; i<accNames.length; i++) {
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
                String[] id = { prevChr, Integer.toString(prevStart), prevRef, prevAlt, geneName, withinSS };
                printScores(id, scores);
            }
        }
        //update info if current line is next variant
        if(nextVariant) {
            prevChr=chr;
            prevStart=start;
            prevRef=ref;
            prevAlt=alt;
            Arrays.fill(scores, -99.0);
        }
        //update scores
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
            geneName = split[3];
            geneEnd = Integer.parseInt(split[2]);
            geneChr = split[0];
            if (split.length>6) {
                String[] donStartStr=split[6].split(",");
                String[] donEndStr=split[7].split(",");
                String[] accStartStr=split[8].split(",");
                String[] accEndStr=split[9].split(",");
                for (int k=0; k<donStartStr.length; k++) {
                    donStart[k]=Integer.parseInt(donStartStr[k]);
                    donEnd[k]=Integer.parseInt(donEndStr[k]);
                }
                for (int k=0; k<accStartStr.length; k++) {
                    accStart[k]=Integer.parseInt(accStartStr[k]);
                    accEnd[k]=Integer.parseInt(accEndStr[k]);
                }
                donNames=split[10].split(",");
                accNames=split[11].split(",");
            }
        }
        //update final line
        finalLineID[0] = chr; finalLineID[1]=Integer.toString(start);
        finalLineID[2]=ref; finalLineID[3]=alt;
        finalLineID[4]="."; finalLineID[5]=".";
        if (start <= geneEnd && geneChr.equals(chr)) {
            finalLineID[4]= geneName;
        }
        for (int i=0; i<finalLineScores.length; i++) {
            finalLineScores[i]= scores[i];
        }
        prevType = type;
    }
    //print final variant
    if (!(prevType.equals("GENE"))) {
        printScores(finalLineID, finalLineScores);
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
        for (int i=0; i<6; i++) {
            System.out.print(id[i]+"\t");   
        }
        for (int i=0; i<11; i++) {
            System.out.print(scores[i]+"\t");
        }
        System.out.println(scores[11]);        
    }

}
