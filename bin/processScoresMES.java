//Description: Processing raw output of MaxEntScan. Takes the many scores generated for most variants
import java.io.*;
import java.util.Arrays;
class processScoresMES {

    public static void main(String[] args) {	
        try {
            BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
            //score tracking variables
            String previousID = "";
            String previousDonOrAcc = "";
            String s = "";
            double highestRef = -99;
            double highestAlt = -99;
            double currentScore = -99;
            int currentPos = -99;
            double[] scores = new double[4];
            int[] positions = new int[4];
            Arrays.fill(scores, -99);
            Arrays.fill(positions, -99);
            //0=highestRef,1=highestAlt,2=2ndHighestRef,3=2ndHighestAlt,4=hrPos,5=haPos,6=2rPos,7=2aPos
            //process lines
            while ((s = in.readLine()) != null) {
                //read score and variant info
                String[] splitScore = s.split("\\s+");
                s = splitScore[0];
                currentScore = Double.parseDouble(splitScore[1]);
                currentPos = Integer.parseInt(s.substring(s.lastIndexOf(";")+1));
                int leading = s.indexOf('>');
                s=s.substring(leading+1);
                String refOrAlt = "REF";
                String donOrAcc = "MESACC";
                if (s.contains("ALT")) {
                    refOrAlt = "ALT";
                }
                if (s.contains("DON")) {
                    donOrAcc = "MESDON";
                }
                //if very first line
                if (previousID.equals("")) {
                    previousID = s;	
                    previousDonOrAcc = donOrAcc;
                    //if current line is next variant 
                } else if (!(s.substring(0,s.lastIndexOf(";")-3).equals(previousID.substring(0,previousID.lastIndexOf(";")+-3)))) {
                        String[] sep = previousID.split(";");
                        String chr = sep[0];
                        int start = Integer.parseInt(sep[1]);
                        String ref = sep[2];
                        int end = start + ref.length() -1;
                        String alt = sep[3];
                        int pos = Integer.parseInt(sep[5]);
                        alt = alt.replace("(", "");
                        System.out.print(chr+"\t"+start+"\t"+ref+"\t"+alt+"\t"+donOrAcc+"\t");
                        for (int i=0; i<4; i++) {
                            System.out.print(Double.toString(scores[i])+"\t");
                        }
                        for (int i=0; i<3; i++) {
                            System.out.print(Integer.toString(positions[i])+"\t");
                        }
                        System.out.println(Integer.toString(positions[3]));
                        Arrays.fill(scores, -99);
                        Arrays.fill(positions, -99);
                        previousID = s;	
                        previousDonOrAcc = donOrAcc;
                }
                //update score
                if (refOrAlt.equals("REF")) {
                    if (currentScore > scores[0]) {
                        scores[2] = scores[0];
                        scores[0] = currentScore;
                        positions[2] = positions[0];
                        positions[0] = currentPos;
                    } else if (currentScore > scores[2] && currentScore != scores[0]) {
                        scores[2] = currentScore;
                        positions[2] = currentPos;
                    }
                } else if (refOrAlt.equals("ALT")) {
                    if (currentScore > scores[1]) {
                        scores[3] = scores[1];
                        scores[1] = currentScore;
                        positions[3] = positions[1];
                        positions[1] = currentPos;
                    } else if (currentScore > scores[3] && currentScore != scores[1]) {
                        scores[3] = currentScore;
                        positions[3] = currentPos;
                    }
                }	
            }
            //print final variant
            String[] sep = previousID.split(";");
            String chr = sep[0];
            int start = Integer.parseInt(sep[1]);
            String ref = sep[2];
            int end = start + ref.length() -1;
            String alt = sep[3];
            int pos = Integer.parseInt(sep[5]);
            alt = alt.replace("(", "");
            System.out.print(chr+"\t"+start+"\t"+ref+"\t"+alt+"\t"+previousDonOrAcc+"\t");
            for (int i=0; i<4; i++) {
                System.out.print(Double.toString(scores[i])+"\t");
            }
            for (int i=0; i<3; i++) {
                System.out.print(Integer.toString(positions[i])+"\t");
            }
            System.out.println(Integer.toString(positions[3]));
            in.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
