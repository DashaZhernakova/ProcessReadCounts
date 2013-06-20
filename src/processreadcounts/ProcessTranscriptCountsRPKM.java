
package processreadcounts;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author dashazhernakova
 */
public class ProcessTranscriptCountsRPKM {
    String bed12Fname;
    String dirName;
    String outFile;
    String fnamePattern;
    String alnFnamePattern;
    ArrayList<String> allTranscr;
    String[] sampleNames;
    float[] readsPerSample;
    boolean[] featureExpressed;
    int trSize;
    int samplesSize;
    boolean normalize;
    boolean normalizeByLength;
    String convertGeneName;
    HashMap<String,Integer> geneLengths;
    ProcessTranscriptCountsRPKM(String dir, String fname, String bed12, String output, String aln){
        bed12Fname=bed12;
	fnamePattern = fname;
        dirName=dir;
	outFile=output;
        allTranscr = new ArrayList<String>();
        alnFnamePattern = aln;
    }
  
    
    /**
     * gets all genes from bed12
     * @return
     * @throws FileNotFoundException
     * @throws IOException 
     */
    public void getAllGenes() throws FileNotFoundException, IOException{

	System.out.println("getting all genes from " + bed12Fname);
        
        allTranscr = new ArrayList();
        geneLengths = new HashMap<String, Integer>();
        TextFile bed12 = new TextFile(bed12Fname, false);
        String[] els = bed12.readLineElems(TextFile.tab);
        String cur_gene = els[3], lengths = els[10], gene;
        int gene_len;
        while( (els = bed12.readLineElems(TextFile.tab)) != null){
            gene = els[3];
            allTranscr.add(gene);
            if (cur_gene.equals(gene)){
                lengths += els[10];
            }
            else{
                gene_len = 0;
                for (String len : lengths.split(","))
                    gene_len += Integer.parseInt(len);
                geneLengths.put(cur_gene, gene_len);
                cur_gene = els[3];
                lengths = els[10];
            }
        }
        bed12.close();
    }
    
    /**
     * gets filenames ending with fnamePattern in dir dirName
     * @return 
     */
    public ArrayList<String> getFileNames(){
	System.out.println("getting all files with reads per transcripts counts. Files names finish with " + fnamePattern);
        ArrayList<String> file_names = new ArrayList<String>();
        File dir = new File(dirName);
        for (File ch : dir.listFiles()) {
            if (ch.isDirectory())
                for (File child : ch.listFiles()){
                   if ( (child.getName().equals(fnamePattern)) ){
                       file_names.add(child.getPath());
                  }
            }
        }
         return file_names;
    }
    
    private boolean avgHigherThanThreshold(float thres, float[] values){
        float sum = 0;
        for (float count : values){
            sum+=count;
        }
        if (sum/values.length > thres)
            return true;
        return false;
    }
    
    public void printTable(float[][] table) throws IOException{
	System.out.println("printing the resulting expression table to " + outFile);
        
        TextFile outReads = new TextFile(outFile, true);
        String curCount="", numReads="",trName="";
        float reads=0;
        int len;
        outReads.write("probe\t");
        outReads.writelnTabDelimited(sampleNames);
        for (int tr = 0; tr < trSize; tr++){            
            
            //if ( (featureExpressed[tr]) && (avgHigherThanThreshold(0, table[tr])) ){
              if (featureExpressed[tr]){
                String tr_id = allTranscr.get(tr);
                if (geneLengths.containsKey(tr_id)){
                    len = geneLengths.get(tr_id);
                    outReads.write(tr_id);
                    for (int sam = 0; sam < samplesSize; sam++){

                        reads = (table[tr][sam] * 1000000000) / readsPerSample[sam] / len;
                        outReads.write("\t" + String.format("%.8f", reads));//reads
                        //if (tr_id.equals("7SK"))
                        //    System.out.println(table[tr][sam] + "\t" + len + "\t" + reads);

                    }
                }
                else{
                    System.out.println("No gene " + tr_id + " in annotation file!");
              }
                
                outReads.writeln();
            }
        }
        outReads.close();
    }
    
    
    public String[][] fill(String[][] table){
        for (int i = 0; i < table.length; i++)
            for (int j = 0; j < table[i].length; j++)
                if ((i == 0) || (j ==0))
                    table[i][j] = "";
                else
                    table[i][j] = "0.0";
        return table;
    }
    
    
    /**
     * gets number of mapped reads from filtered sam file
     * @param path - path to tophat output
     * @return
     * @throws IOException 
     */
    public int getNumMappedReads(String path){
        int numReads = 0;
        try {
            TextFile filtSam;
            if (alnFnamePattern.isEmpty())
                filtSam = new TextFile(path.replace(fnamePattern, "accepted_hits.filtered.sam"), false);
            else
                filtSam = new TextFile(path.replace(fnamePattern, alnFnamePattern), false);
            String line ="";
            
            while ( (line = filtSam.readLine()) != null){
                if (! line.startsWith("@")){
                    numReads++;

                }
            }
            filtSam.close();
            //System.out.println("number of mapped reads: " + numReads);
            
        } catch (IOException ex) {
            Logger.getLogger(processreadcounts.ProcessTranscriptCounts.class.getName()).log(Level.SEVERE, null, ex);
        }
        return numReads;
    }
    /**
     * writes feature counts into expression table. New version - with float[][]
     * @throws IOException 
     */
    public void readCounts() throws IOException{
        System.out.println("started processing " + dirName);
	TextFile in = null;
	String line = "";
        String[] splLine;
        
        int index = 0;
        
        //all genes from bed12, sorted
        //allTranscr = getAllTranscripts();
        getAllGenes();
        Collections.sort(allTranscr);
        
        //all files to use
        ArrayList<String> fNames = getFileNames();
        
        trSize = allTranscr.size(); //number of transcripts
        samplesSize = fNames.size(); //number of samples
        //table with only counts
        float[][] table = new float[trSize][samplesSize];
        readsPerSample = new float[samplesSize];
        sampleNames = new String[samplesSize];
        featureExpressed = new boolean[trSize];
        //table=fill(table);//fill with empty strings
        
        int curId = 0;
        System.out.println("started generating the expression table");
        String trId;
        for (String fName : fNames){
            int rCount=0;
            in = new TextFile(fName, false);
            //getting sample name from the file name
            System.out.println("processing file: " + fName);
            String sampleId="";
            String[] splName = fName.split("/");
            sampleId = splName[splName.length - 2];
            
            sampleNames[curId] = sampleId;
            
            //get number of mapped reads from filtered sam file
            readsPerSample[curId] = getNumMappedReads(fName);
            //in.readLine();
            while ((line = in.readLine()) != null){
                    splLine = line.split("\t");
                    trId = splLine[0];
                    if (trId != null){
                        float count = Float.valueOf(splLine[1]);
                        //System.out.println(trId);

                        //get the transcript position in all transcripts list 
                        index = Collections.binarySearch(allTranscr, trId);
                        
                        if (index >= 0){ //if transcript found
                            if (count > 0){
                                featureExpressed[index] = true; //set this feature as expressed
                                rCount+=count; //junk
                                table[index][curId] += count;
                            }
                        }   
                    }
            }
            curId++;
            System.out.println("overall number of reads=" + rCount);
            System.out.println("overall number of mapped reads=" + readsPerSample[curId - 1]);
            in.close();
            
        }
        //print the result
        printTable(table);
       
    }
    
    public static void main(String[] args) throws IOException {
        ProcessTranscriptCountsRPKM p;
        //if (args.length > 3){
            //p = new ProcessTranscriptCounts(args[0], args[1], args[2], args[3]);
        p = new processreadcounts.ProcessTranscriptCountsRPKM("/Users/dashazhernakova/Documents/UMCG/data/lincRNA_Sebo/mappedData/", 
                "intersectBed_gene_counts.txt.gz", 
                "/Users/dashazhernakova/Documents/UMCG/hg19/exonic_genes_v69_stranded_cut.bed12", 
                "/Users/dashazhernakova/Documents/tmp.txt",
                "ikzf1.sam");
        //p = new ProcessTranscriptCounts("/Users/dashazhernakova/Documents/UMCG/pickrell_data/tophat_out/", "covBed_exon_counts_out.txt", "/Users/dashazhernakova/Documents/UMCG/hg19/annotation_transcr_hg19.txt", "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/tmp_expression.txt");
                
        p.readCounts();
        //p.readCounts();
       // }
        //else
          //  System.out.println("not enough arguments given");
        
    }
}




