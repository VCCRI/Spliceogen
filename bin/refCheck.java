//Description: for generating bed intervals as input to bedtools getfasta
import java.io.*;

class refCheck  {

public static void main(String[] args) {	
    try {
    BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
    String s = "";
    //process lines
    while ((s = in.readLine()) != null) {
        String[] split = s.split("\\s+");
        int prefixStart= Integer.parseInt(split[2]) - 91;
        int suffixEnd = Integer.parseInt(split[2]) + 89 + split[4].length();
   	    System.out.println(split[1]+"\t"+Integer.toString(prefixStart)+"\t"+Integer.toString(suffixEnd)+"\t"+split[1]+";"+split[2]+";"+split[4]+";"+split[5]+";"+split[3]+";\t1\t"+split[3]);
    }
    in.close();
    } catch (Exception e) {
        System.err.println("Error: " + e.getMessage());
    }
}
}
