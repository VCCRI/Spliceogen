//Description: Processing raw output of MaxEntScan. Takes the many scores generated for most variants
import java.io.*;
class processScoresMES {

public static void main(String[] args) {	
    try {
    BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
	//score tracking variables
	String previousID = "";
    String s = "";
	double highestRef = -99;
	double highestAlt = -99;
    double currentScore = -99;
    //process lines
    while ((s = in.readLine()) != null) {
	    //read score and variant info
        String[] splitScore = s.split("\\s+");
        s = splitScore[0];
	    currentScore = Double.parseDouble(splitScore[1]);
	    int leading = s.indexOf('h');
        s=s.substring(leading-1);
    	String refOrAlt = s.substring(s.length()-3);
        String donOrAcc = "MES"+s.substring(s.length()-6, s.length()-3);
    	//if not very first line
    	if (!(previousID.equals(""))) {
    	    //if current line is next variant 
    	    if (!(s.substring(0, s.length()-3).equals(previousID.substring(0,previousID.length()-3)))) {
                String[] sep = previousID.split(";");
                String chr = sep[0];
                int start = Integer.parseInt(sep[1]);
                String ref = sep[2];
                int end = start + ref.length() -1;
                String alt = sep[3];
                alt = alt.replace("(", "");
    		    System.out.println(chr+"\t"+start+"\t"+ref+"\t"+alt+"\t"+donOrAcc+"\t"+highestRef+"\t"+highestAlt);
    	        highestRef = -99;
                highestAlt = -99;
    	    }
        }
    	//update score
    	if (refOrAlt.equals("REF")) {
    		if (currentScore > highestRef) {
    			highestRef = currentScore;
    		}
    	} else
    	if (refOrAlt.equals("ALT")) {
    		if (currentScore > highestAlt) {
    			highestAlt = currentScore;
    		}
    	}	
    	previousID = s;	
    }
   	//print final variant
    String[] sep = s.split(";");
    String chr = sep[0];
    int start = Integer.parseInt(sep[1]);
    String ref = sep[2];
    int end = start + ref.length() -1;
    String alt = sep[3];
    alt = alt.replace("(", "");
    String donOrAcc = "MES"+s.substring(s.length()-6, s.length()-3);
   	System.out.println(chr+"\t"+start+"\t"+ref+"\t"+alt+"\t"+donOrAcc+"\t"+highestRef+"\t"+highestAlt);
    in.close();
    } catch (Exception e) {
        e.printStackTrace();
    }
}
}
