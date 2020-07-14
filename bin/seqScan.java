//Steve Monger 21.06.18
//Preparation of ref and alt strings for input to MaxEntScan from FASTA file
//Scans for any GT/AG dinucleotides such that ref/alt falls anywhere within the 9/23bp donor/acceptor motif
//Handles SNPs and indels
import java.io.*;

class seqScan {

    public static String[] mesOutputDon = new String[15000];
    public static String[] mesOutputAcc = new String[40000];
    public static String[] gsOutput = new String[30000];
    public static String[] ESRoutput = new String[10000];
    public static int gsIndex = 0;
    public static int mesIndexDon = 0;
    public static int mesIndexAcc = 0;
    public static int ESRoutputIndex = 0;
    public static double[] ESEtable = createESRarrays("data/ESE.txt");
    public static double[] ESStable = createESRarrays("data/ESS.txt");
    public static boolean useESR= true;
    public static String fileID="";
    public static boolean unsplitVCFfound = false;

    public static void main(String[] args) {	

        String header = "";
        String s = "";
        if (args.length < 3) {
            System.out.println("must provide input file and -useESR/skipESR and fileID as command line argument");
            System.exit(1);
        }
        if (args[1]=="-skipESR") {
            useESR = false;
        }
        fileID=args[2];
        try{
            fillESRtables();
            //open input file
            FileInputStream fstream = new FileInputStream(args[0]);
            DataInputStream in = new DataInputStream(fstream);
            BufferedReader br = new BufferedReader(new InputStreamReader(in));
            //process file 2 lines at a time
            while ((header = br.readLine()) != null) {
                s = br.readLine();
                try {
                    scanStrings(header, s, args[2]);
                }
                catch (Exception e){
                    e.printStackTrace();   
                }    
            }
            in.close();
        }catch (Exception e){
            System.err.println("Error: " + e.getMessage());
        }
        //write output to files
        appendToFiles();
    }

    public static void scanStrings (String header, String s, String fileID) {
        /* SCANNING EXPLANATION
           For MaxEntScan, a region of 45bp with variant at middle nucleotide 23 is needed (anywhere ref/alt can contribute to the motif).
           Scan for acceptor AG from nucleotide pos 19 to 41, and for donor GT from pos 18 to 26 and generate query strings for scoring
           The additional 68bp flanking sequence is necessary for GeneSplicer and filling out MaxEntScan alt strings shortened by deletions.
           For indels, alt strings have a different length. Adjustment "correct1" is made so that only ref/alt equivalent positions are scanned.
correct1: only needed on ref string in event of deletions on antisense strand, because start value is sense oriented eg. GAT>GT is represented as 5'-GA>G-3' not 3'-TA>T-5'.
correct2: adjusts for empty alt strings, arising from deletions denoted for example as ref=A, alt=-or*
         */
        int flank = 68;
        int correct1 = 0;
        s = s.toLowerCase();
        String id = header.substring(0, header.indexOf(":"));
        String[] sep = id.split(";");
        String chr = sep[0].substring(1);
        String start = sep[1];
        String ref = sep[2].toUpperCase();
        String alt = sep[3].toUpperCase();
        String strand = sep[4];
        String mesAdd = "";
        double maxRefESE = 0;
        double maxAltESE = 0;
        double minRefESS = 0;
        double minAltESS = 0;
        //ignore deletions larger than flanking region
        if (ref.length()>35 || alt.length()>35) {
            return;
        }
        if (strand.equals("+")) {
            strand = "POS";
        }
        //index correction for deletions notated as alt="*"
        String refName = ref;
        String altName = alt;
        int correct2=0;
        if (alt.equals("-") || alt.equals(".")) {
            alt="*";
        }
        if (alt.equals("*")) {
            alt="";
            correct2=1;
        }
        //catch unsplit VCFs (ie. multiple alt alleles represented on one comma separated line)
        if (alt.contains(",")) {
            if (!unsplitVCFfound) {
                System.out.println("Note: unsplit VCF/BED record(s) were detected, ie. variants with multiple alt alleles represented on a single line. Such records are susceptible to false negatives, so we recommend splitting these records onto multiple lines");
                unsplitVCFfound = true;
            }
            String[] splitAlt = alt.split("");
            alt = splitAlt[0];
        }
        //if negative strand
        if (strand.equals("-")) {
            //reverse complement alt
            id = id.concat("(-)");
            ref = getRC(ref);
            alt = getRC(alt);
            correct1 = ref.length()-1;
            strand = "NEG";
        }
        //generate alt string
        String altSeq = s.substring(0, 22 + flank) + alt + s.substring(22 + flank + ref.length());		
        altSeq = altSeq.toLowerCase();
        //retrieve reference allele from fasta to check for mismatch with provided ref allele
        String fastaRef = s.substring(22 + flank, 22 + flank + ref.length()).toUpperCase();
        //ignore variants containing invalid fasta chars such as "n" which cause MaxEntScan to exit prematurely
        String nTest = altSeq.substring(17+flank-correct1-18, 41+flank-correct1+6+alt.length()).concat(ref).toUpperCase();
        //check fasta reference matches provided ref allele
        if (!fastaRef.equals(ref)) {
            appendOmmittedRefMismatch(id);    
        } else {
            //genesplicer fasta strings
            gsOutput[gsIndex]=">".concat(strand).concat(chr).concat(":").concat(start).concat(":").concat(refName).concat(":").concat(altName).concat(":GSREF");
            gsIndex++;
            gsOutput[gsIndex]=s;
            gsIndex++;
            gsOutput[gsIndex]=">".concat(strand).concat(chr).concat(":").concat(start).concat(":").concat(refName).concat(":").concat(altName).concat(":GSALT");
            gsIndex++;
            gsOutput[gsIndex]=altSeq;
            gsIndex++;
            //ESRseq
            if(useESR) {
                for (int i=17+flank; i<23+flank; i++) {
                    int refIndex = getESRindex(s.substring(i, i+6));
                    int altIndex = getESRindex(altSeq.substring(i, i+6));
                    if (ESEtable[refIndex]>maxRefESE) {
                        maxRefESE = ESEtable[refIndex];
                    }
                    if (ESEtable[altIndex]>maxAltESE) {
                        maxAltESE = ESEtable[altIndex];
                    }
                    if (ESStable[refIndex]<minRefESS) {
                        minRefESS = ESStable[refIndex];
                    }
                    if (ESStable[altIndex]<minAltESS) {
                        minAltESS = ESStable[altIndex];
                    }
                    if (i==22+flank) {
                        //ESR deletions
                        if (ref.length()>1 || correct2==1) {
                            int j=1;
                            do {
                                refIndex = getESRindex(s.substring(i+j, i+j+6));
                                if (ESEtable[refIndex]>maxRefESE) {
                                    maxRefESE = ESEtable[refIndex];
                                }
                                if (ESStable[refIndex]<minRefESS) {
                                    minRefESS = ESStable[refIndex];
                                }
                                j++;
                            }
                            while (j<ref.length());
                        }
                        //ESR insertions
                        if (alt.length()>1) {
                            int j=1;
                            do {
                                altIndex = getESRindex(altSeq.substring(i+j, i+j+6));
                                if (ESEtable[altIndex]>maxAltESE) {
                                    maxAltESE = ESEtable[altIndex];
                                }
                                if (ESStable[altIndex]<minAltESS) {
                                    minAltESS = ESStable[altIndex];
                                }
                                j++;
                            }
                            while (j<ref.length());
                        }
                    }
                }
            }
            //MaxEntScan
            if (containsIllegalFastaChar(nTest, nTest.length())) {
                appendOmmittedInvalidFasta(id);
            } else {
                //scan ref sequence
                for (int i=17+flank-correct1; i<41+flank-correct1+ref.length()-1; i++) {
                    //if acceptor
                    if (s.substring(i, i+2).equals("ag") && i!=17+flank-correct1) {
                        //output ref and alt acceptor strings
//test
                        int pos = flank-i-correct1+41; 
                        mesOutputAcc[mesIndexAcc]=s.substring(i-18, i+5).concat(id).concat("ACC").concat("REF;").concat(Integer.toString(pos));
                        mesIndexAcc++;
                        mesOutputAcc[mesIndexAcc]=altSeq.substring(i-18, i+5).concat(id).concat("ACC").concat("ALT;").concat(Integer.toString(pos));
                        mesIndexAcc++;
                        //for indels, output offset alt acceptor string
                        if (ref.length()>1 || correct2==1) {
                            int j=1;
                            do {                
                                mesOutputAcc[mesIndexAcc]=altSeq.substring(i-18-j, i+5-j).concat(id).concat("ACC").concat("ALT;").concat(Integer.toString(pos+j));
//test
                                mesIndexAcc++;
                                j++;
                            }
                            while (j<ref.length()+correct2 && i>=j+18);
                        } else
                            if (alt.length()>1) {
                                int j=1;
                                do {
//test
                                    mesOutputAcc[mesIndexAcc]=altSeq.substring(i-18+j, i+5+j).concat(id).concat("ACC").concat("ALT;").concat(Integer.toString(pos+j));
                                    mesIndexAcc++;
                                    j++;
                                }
                                while (j<alt.length());
                            }
                    }	
                    //if donor
                    if (s.substring(i, i+2).equals("gt") && i<26+correct1+flank+ref.length()-1) {
                        //output ref and alt donor strings
                        int pos = flank-i-correct1+32; 
                        mesOutputDon[mesIndexDon]=s.substring(i-3, i+6).concat(id).concat("DON").concat("REF;").concat(Integer.toString(pos));
                        mesIndexDon++;
                        mesOutputDon[mesIndexDon]=altSeq.substring(i-3, i+6).concat(id).concat("DON").concat("ALT;").concat(Integer.toString(pos));
                        mesIndexDon++;
                        //for indels, output offset alt donor string
                        if (ref.length()>1 || correct2==1) {
                            int j=1;
                            do {                
                                mesOutputDon[mesIndexDon]=altSeq.substring(i-3-j, i+6-j).concat(id).concat("DON").concat("ALT;").concat(Integer.toString(pos+j));
                                mesIndexDon++;
                                j++;
                            }
                            while (j<ref.length()+correct2  && i >= j+3);
                        } else
                            if (alt.length()>1) {
                                int j=1;
                                do {
                                    mesOutputDon[mesIndexDon]=altSeq.substring(i-3+j, i+6+j).concat(id).concat("DON").concat("ALT;").concat(Integer.toString(pos+j));
                                    mesIndexDon++;
                                    j++;
                                }
                                while (j<alt.length());
                            }
                    }	
                }
                //scan alt sequence
                for (int i=17+flank; i<41+flank+alt.length()-1; i++) {
                    //if acceptor
                    if (altSeq.substring(i, i+2).equals("ag") && i!=17+flank) {
                        //output ref and alt acceptor strings
//test
                        int pos = flank-i-correct1+41;
                        mesOutputAcc[mesIndexAcc]=altSeq.substring(i-18, i+5).concat(id).concat("ACC").concat("ALT;").concat(Integer.toString(pos));
                        mesIndexAcc++;
                        mesOutputAcc[mesIndexAcc]=s.substring(i-18, i+5).concat(id).concat("ACC").concat("REF;").concat(Integer.toString(pos));
                        mesIndexAcc++;
                        //for indels, output offset ref acceptor string
                        if (alt.length()>1) {
                            int j=1;
                            do {
//test
                                mesOutputAcc[mesIndexAcc]=s.substring(i-18-j, i+5-j).concat(id).concat("ACC").concat("REF;").concat(Integer.toString(pos+j));
                                mesIndexAcc++;
                                j++;
                            }
                            while (j<alt.length() && i >= j+18);
                        }
                        if (ref.length()>1||correct2==1) {
                            int j=1;
                            do {
                                mesOutputAcc[mesIndexAcc]=s.substring(i-18+j, i+5+j).concat(id).concat("ACC").concat("REF;").concat(Integer.toString(pos+j));
                                mesIndexAcc++;
                                j++;
                            }
                            while (j<ref.length()+correct2);
                        }
                    }
                    //if donor
                    if (altSeq.substring(i, i+2).equals("gt") && i<26+correct1+flank+alt.length()-1) {
                        int pos = flank-i-correct1+32; 
                        //output ref and alt donor strings
                        mesOutputDon[mesIndexDon]=altSeq.substring(i-3, i+6).concat(id).concat("DON").concat("ALT;").concat(Integer.toString(pos));
                        mesIndexDon++;
                        mesOutputDon[mesIndexDon]=s.substring(i-3, i+6).concat(id).concat("DON").concat("REF;").concat(Integer.toString(pos));
                        mesIndexDon++;
                        //for indels, output offset ref donor string
                        if (alt.length()>1) {
                            int j=1;
                            do {
                                mesOutputDon[mesIndexDon]=s.substring(i-3-j, i+6-j).concat(id).concat("DON").concat("REF;").concat(Integer.toString(pos+j));
                                mesIndexDon++;
                                j++;
                            }
                            while (j<alt.length() && i >= j+3);
                        }
                        if (ref.length()>1 || correct2==1) {
                            int j=1;
                            do {
                                mesOutputDon[mesIndexDon]=s.substring(i-3+j, i+6+j).concat(id).concat("DON").concat("REF;").concat(Integer.toString(pos+j));
                                mesIndexDon++;
                                j++;                    
                            }
                            while (j<ref.length()+correct2);
                        }
                    }
                }
            }	
        }
        //Format ESR scores and ouput
        String[] ESRoutputScores = { Double.toString(maxRefESE), Double.toString(maxAltESE), Double.toString(minRefESS), Double.toString(minAltESS) };
        for (int k=0; k<4; k++) {
            if (ESRoutputScores[k].equals("0.0")) {
                ESRoutputScores[k]=".";
            }
        }    
        ESRoutput[ESRoutputIndex]=chr+"\t"+start+"\t"+refName+"\t"+altName+"\t"+"ESR"+"\t"+Double.toString(maxRefESE)+"\t"+Double.toString(maxAltESE)+"\t"+Double.toString(minRefESS)+"\t"+Double.toString(minAltESS);
        ESRoutputIndex++;
        //append arrays to file then reset
        if (mesIndexAcc > 25000 || mesIndexDon > 9000 ) {
            appendToFiles();
            resetOutputArrays();
        } else
            if (gsIndex > 13000 || ESRoutputIndex > 3250 ) {
                appendToFiles();
                resetOutputArrays();
            }
    }

    //get reverse complement
    public static String getRC (String query) {
        String rc = "";
        for (int i=query.length()-1; i>-1; i--) {
            if (query.charAt(i)=='A') {
                rc = rc.concat("T");
            }
            if (query.charAt(i)=='C') {
                rc = rc.concat("G");
            }
            if (query.charAt(i)=='G') {
                rc = rc.concat("C");
            }
            if (query.charAt(i)=='T') {
                rc = rc.concat("A");
            }
        }
        return rc;
    }

    //generate reference arrays for ESRseq scores, with all potential hexamers mapping to a unique index assigned by getESRindex()
    public static double[] createESRarrays (String fileName) {
        double[] ESRscores = new double[4096];
        try {
            File ESRfile = new File(fileName);
            BufferedReader br = new BufferedReader(new FileReader(ESRfile));
            String line = "";
            while ((line = br.readLine()) != null) {
                String[] split = line.split("\\s+");
                ESRscores[getESRindex(split[0])]=Double.parseDouble(split[1]);
            }
            br.close();
        }catch (Exception e){
            System.err.println("Error: " + e.getMessage());
        }
        return ESRscores;
    }

    public static void fillESRtables() {
        try {
            //ESE
            File ESEfile = new File("data/ESE.txt");
            BufferedReader ESEbr = new BufferedReader(new FileReader(ESEfile));
            String line = "";
            while ((line = ESEbr.readLine()) != null) {
                String[] split = line.split("\\s+");
                ESEtable[getESRindex(split[0])]=Double.parseDouble(split[1]);
            }
            ESEbr.close();
            //ESS
            File ESSfile = new File("data/ESS.txt");
            BufferedReader ESSbr = new BufferedReader(new FileReader(ESSfile));
            line = "";
            while ((line = ESSbr.readLine()) != null) {
                String[] split = line.split("\\s+");
                ESStable[getESRindex(split[0])]=Double.parseDouble(split[1]);
            }
            ESSbr.close();
        }catch (Exception e){
            System.err.println("Error: " + e.getMessage());
        }
    }

    //convert hexamer sequence to integer position in ESR reference score array
    public static int getESRindex (String seq) {
        seq = seq.toUpperCase();
        double val=0;
        for (int i=0; i<6; i++) {
            char c = seq.charAt(i);
            switch(c) {
                case 'A': val = val + 0*Math.pow(4,i);
                          break;
                case 'C': val = val + 1*Math.pow(4,i);
                          break;
                case 'G': val = val + 2*Math.pow(4,i);
                          break;
                case 'T': val = val + 3*Math.pow(4,i);
                          break;
            }
        }
        return (int)val;
    }

    public static void resetOutputArrays() {
        mesOutputAcc = new String[40000];
        mesOutputDon = new String[15000];
        gsOutput = new String[30000];
        ESRoutput = new String[10000];
        gsIndex = 0; mesIndexDon = 0; mesIndexAcc = 0; ESRoutputIndex = 0;
        System.gc();
    }

    public static boolean containsIllegalFastaChar(String s, int i) {
        s = s.toUpperCase();
        for (int k=0; k<i; k++) {
            switch(s.charAt(k))
            {
                case 'A': break;
                case 'C': break;
                case 'G': break;
                case 'T': break;
                default: return true;
            }
        }
        return false;
    }

    public static void appendToFiles () {
        String gsName = "temp/"+fileID+"gsInput.FASTA";
        String mesAccName = "temp/"+fileID+"mesAcceptorInput.txt";
        String mesDonName = "temp/"+fileID+"mesDonorInput.txt";
        String ESRname = "temp/"+fileID+"ESRoutput.txt";
        try {
            //write GeneSplicer input to file
            FileWriter fwGS = new FileWriter(gsName, true);
            BufferedWriter gsWriter = new BufferedWriter(fwGS);
            for (int i=0; i<gsIndex; i++) {
                gsWriter.write(gsOutput[i]+"\n");
            }
            gsWriter.close();
            //write maxEntScan Acceptor input to file
            FileWriter fwMESacc = new FileWriter(mesAccName, true);
            BufferedWriter mesAccWriter = new BufferedWriter(fwMESacc);
            for (int i=0; i<mesIndexAcc; i++) {
                if (!containsIllegalFastaChar(mesOutputAcc[i], 23)) {
                    mesAccWriter.write(mesOutputAcc[i]+"\n");
                }
            }
            mesAccWriter.close();
            //write maxEntScan Donor input to file
            FileWriter fwMESdon = new FileWriter(mesDonName, true);
            BufferedWriter mesDonWriter = new BufferedWriter(fwMESdon);
            for (int i=0; i<mesIndexDon; i++) {
                if (!containsIllegalFastaChar(mesOutputDon[i], 9)) {
                    mesDonWriter.write(mesOutputDon[i]+"\n");
                }
            }
            mesDonWriter.close();
            //write ESR output to file
            FileWriter fwESR = new FileWriter(ESRname, true);
            BufferedWriter ESRwriter = new BufferedWriter(fwESR);
            for (int i=0; i<ESRoutputIndex; i++) {
                ESRwriter.write(ESRoutput[i]+"\n");
            }
            ESRwriter.close();
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }
    }

    public static void appendOmmittedInvalidFasta(String id) {
        try {
            String fileName = "output/"+fileID+"mesOmmitted.txt";
            FileWriter fw = new FileWriter(fileName, true);
            BufferedWriter writer = new BufferedWriter(fw);
            writer.write(id+"\n");
            writer.close();
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }
    }

    public static void appendOmmittedRefMismatch(String id) {
        try {
            String fileName = "output/"+fileID+"refMismatch.txt";
            FileWriter fw = new FileWriter(fileName, true);
            BufferedWriter writer = new BufferedWriter(fw);
            writer.write(id+"\n");
            writer.close();
        } catch (IOException e) {
            System.out.println(e.getMessage());
        }
    }

}
