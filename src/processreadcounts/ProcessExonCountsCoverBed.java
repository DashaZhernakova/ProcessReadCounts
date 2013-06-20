
package processreadcounts;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;


public class ProcessExonCountsCoverBed {

    public void process(String fName) throws FileNotFoundException, IOException{
        System.out.println("Starting to sum up reads for transcripts in " + fName);
        System.out.println("Writing output to " + fName.replace(".txt", "_out.txt"));
        BufferedReader exCounts=new BufferedReader(new FileReader(fName));
        BufferedWriter trCounts=new BufferedWriter(new FileWriter(fName.replace(".txt", "_out.txt")));
        String line="", curTranscr="", transcr="";
        float count;
        float curCount=0;
        while ((line = exCounts.readLine()) != null){
            String[] fields = line.split("\t");
            count = Float.valueOf(fields[1])*Float.valueOf(fields[4]);
            //count = Float.valueOf(fields[1]);
            transcr = fields[0].split(":")[0];
            if (transcr.equals(curTranscr))
                curCount += count;
            else{
                if (!curTranscr.isEmpty())
                    trCounts.write(curTranscr + "\t" + Float.toString(curCount)+"\n");
                curCount=count;
                curTranscr=transcr;
            }
        }
        trCounts.write(curTranscr + "\t" + Float.toString(curCount));

        trCounts.close();
        exCounts.close();
        
        //System.out.println("Finished summing up");
    }
    
    public void processExonWise(String fName) throws FileNotFoundException, IOException{
        System.out.println("Starting to sum up reads for exons in " + fName);
        System.out.println("Writing output to " + fName.replace(".txt", "_exonwise_out.txt"));
        BufferedReader in=new BufferedReader(new FileReader(fName));
        BufferedWriter out=new BufferedWriter(new FileWriter(fName.replace(".txt", "_exonwise_out.txt")));
        String line="";
         
        String[] fields;
        while ( (line = in.readLine()) != null ){
           fields = line.split("\t");
           out.write(fields[0] + "\t" + String.valueOf(Float.valueOf(fields[1])*Float.valueOf(fields[4])) + "\n");
           
        }
        in.close();
        out.close();
        //System.out.println("Finished summing up");
    }
    public static void main(String[] args) throws IOException {
        //System.out.println(args[0]);
        //System.out.println(args[1]);
        ProcessExonCountsCoverBed p = new ProcessExonCountsCoverBed();
        //p.process("/Users/dashazhernakova/Documents/UMCG/cluster/tmp.txt");
        
        if (args.length > 1){
            if (args[0].equals("exonwise"))    
                p.processExonWise(args[1]);
        }
        else if (args.length == 1)
                p.process(args[0]);
        else{
                System.out.println("no input file given");
                //p.process("/Users/dashazhernakova/Documents/UMCG/hapmapData/tophat/1382_1_1/covBed_cut_new.txt");
        }   
    }
}
