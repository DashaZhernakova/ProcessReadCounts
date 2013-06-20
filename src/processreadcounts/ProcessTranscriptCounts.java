
package processreadcounts;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import processtmap.GeneNameConverter;
import umcg.genetica.io.text.TextFile;

/**
 * processes files of the type:
 * trId\tcount
 * trId\tcount
 * 
 * generates a table of ids vs trCounts
 * ...
 * @author dashazhernakova
 */
public class ProcessTranscriptCounts {
    
    String annotationFname;
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
    HashMap<String,Integer> trLengths;
    ProcessTranscriptCounts(String dir, String fname, String annot, String output, boolean norm, String conv){
        annotationFname=annot;
	fnamePattern = fname;
        dirName=dir;
	outFile=output;
        allTranscr = new ArrayList<String>();
        normalize = norm;
        convertGeneName = conv;
    }
    
    ProcessTranscriptCounts(String dir, String fname, String annot, String output, String aln, boolean norm, String conv){
        annotationFname=annot;
	fnamePattern = fname;
        dirName=dir;
	outFile=output;
        allTranscr = new ArrayList<String>();
        alnFnamePattern = aln;
        normalize = norm;
        convertGeneName = conv;
    }
    ProcessTranscriptCounts(String dir, String fname, String annot, String output, boolean norm){
        annotationFname=annot;
	fnamePattern = fname;
        dirName=dir;
	outFile=output;
        allTranscr = new ArrayList<String>();
        normalize = norm;
        convertGeneName = null;
    }
    
    ProcessTranscriptCounts(String dir, String fname, String annot, String output, String aln, boolean norm){
        annotationFname=annot;
	fnamePattern = fname;
        dirName=dir;
	outFile=output;
        allTranscr = new ArrayList<String>();
        alnFnamePattern = aln;
        normalize = norm;
        convertGeneName = null;
    }

    ProcessTranscriptCounts(String dir, String fname, String annot, String output, String aln, boolean norm, boolean lenNorm){
        annotationFname=annot;
	fnamePattern = fname;
        dirName=dir;
	outFile=output;
        allTranscr = new ArrayList<String>();
        alnFnamePattern = aln;
        normalize = norm;
        convertGeneName = null;
        normalizeByLength = lenNorm;
    }
    
    /**
     * gets all transcripts from Ensembl
     * @return
     * @throws FileNotFoundException
     * @throws IOException 
     */
    public void getAllTranscripts() throws FileNotFoundException, IOException{

	System.out.println("Getting all transcripts from " + annotationFname);
        TextFile in = new TextFile(annotationFname, false);
        allTranscr = new ArrayList();
        String line = in.readLine();
        int len = 0, cnt = 0;
        String[] spl;
        if (! normalizeByLength){
            while ((line = in.readLine()) != null){
                allTranscr.add(line.split("\t")[1]);
                cnt++;
            }
        }
        else{
            trLengths = new HashMap<String, Integer>();
            while ((line = in.readLine()) != null){
                //System.out.println(line);
                spl = line.split("\t");
                allTranscr.add(spl[1]);
                len = Integer.parseInt(spl[5]) - Integer.parseInt(spl[4]);
                trLengths.put(spl[1], len);
                cnt++;
            }
        }
        in.close();
        System.out.println(cnt + " transcripts / genes read.\n");
        //return allTranscr;
    }
    
    /**
     * gets filenames ending with fnamePattern in dir dirName
     * @return 
     */
    public ArrayList<String> getFileNames(){
	System.out.println("Getting all files with reads per transcripts counts. Files names finish with " + fnamePattern);
        ArrayList<String> file_names = new ArrayList<String>();
        File dir = new File(dirName);
        for (File ch : dir.listFiles()) {
            if (ch.isDirectory())
                for (File child : ch.listFiles()){
                   if (child.getName().equals(fnamePattern)){
                       file_names.add(child.getPath());
                  }
            }
        }
        System.out.println(file_names.size() + " files will be processed.\n");
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
        System.out.println("Normalize=" + normalize);
	//BufferedWriter outReads = new BufferedWriter (new FileWriter(outFile));
	TextFile outReads = new TextFile(outFile, true);
        String curCount="", numReads="",trName="";
        float reads=0;
       
        outReads.write("probe\t");
        outReads.writelnTabDelimited(sampleNames);
        for (int tr = 0; tr < trSize; tr++){            

            //if ( (featureExpressed[tr]) && (avgHigherThanThreshold(0, table[tr])) ){
              if (featureExpressed[tr]){
                outReads.write(allTranscr.get(tr));
                String tr_id = allTranscr.get(tr);
                
                for (int sam = 0; sam < samplesSize; sam++){
                    
                        if (normalize)
                            reads = (table[tr][sam] * 1000000) / readsPerSample[sam];
                        else
                            reads = table[tr][sam];
                        //reads = table[tr][sam];
                        if (! normalize){
                            outReads.write("\t" + String.format("%.3f", reads));//reads
                        }
                        else
                            outReads.write("\t" + String.format("%.8f", reads));//reads
                    
                }
                outReads.writeln();
            }
        }
        outReads.close();
    }
    
    public void printTableNormByLength(float[][] table) throws IOException{
	System.out.println("printing the resulting expression table to " + outFile);
        System.out.println("Normalize=" + normalize);
	//BufferedWriter outReads = new BufferedWriter (new FileWriter(outFile));
	TextFile outReads = new TextFile(outFile, true);
        String curCount="", numReads="",trName="";
        String tr_id;
        float reads=0;
        outReads.write("probe\t");
        outReads.writelnTabDelimited(sampleNames);
        for (int tr = 0; tr < trSize; tr++){            
            //if ( (featureExpressed[tr]) && (avgHigherThanThreshold(1, table[tr])) ){
            if (featureExpressed[tr]){    
                tr_id = allTranscr.get(tr);
                outReads.write(tr_id);
                if (normalize){
                    for (int sam = 0; sam < samplesSize; sam++){
                        reads = (table[tr][sam] * 1000000) / (readsPerSample[sam] * trLengths.get(tr_id));
                        outReads.write("\t" + String.format("%.8f", reads));//reads

                    }
                }
                else{
                    for (int sam = 0; sam < samplesSize; sam++){
                        
                        reads = (table[tr][sam] * 1000000) / trLengths.get(tr_id);
                        outReads.write("\t" + String.format("%.8f", reads));//reads
                        
                    }
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
            //TextFile filtSam = new TextFile(path.replace(fnamePattern, "accepted_hits.filtered.sam"), false);
            //filtSam.open();
            //BufferedReader filtSam = new BufferedReader(new FileReader(path.replace(fnamePattern, "accepted_hits.filtered.sam")));
            TextFile filtSam;
            if (alnFnamePattern.isEmpty())
                filtSam = new TextFile(path.replace(fnamePattern, "accepted_hits.filtered.sam"), false);
            else if (alnFnamePattern.endsWith(".sam")){
                filtSam = new TextFile(path.replace(fnamePattern, alnFnamePattern), false);
                String line ="";

                while ( (line = filtSam.readLine()) != null){
                    if (! line.startsWith("@")){
                        numReads++;
                    }
                }
                filtSam.close();
            }
            else if (alnFnamePattern.endsWith("idxstats")){
                TextFile stats = new TextFile(path.replace(fnamePattern, alnFnamePattern), false);
                
                String[] els;
                while ((els = stats.readLineElems(TextFile.tab)) != null){
                    numReads += Integer.parseInt(els[2]);
                }
                stats.close();
            }
            //System.out.println("number of mapped reads: " + numReads);
            
        } catch (IOException ex) {
            Logger.getLogger(ProcessTranscriptCounts.class.getName()).log(Level.SEVERE, null, ex);
            System.out.println("No alignment file or summary file (" + alnFnamePattern + ") for " + path);
        }
        return numReads;
    }
    
    private String[] readExpressionLine(String line){
        String[] spl;
        String count;
        if (fnamePattern.endsWith(".gtf")){
            spl = line.split("\"");
            count = spl[6].split(" ")[2].replace(";", "");
            return new String[]{spl[1], count};
        }
        return line.split("\t");
        
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
        
        //all transcripts from ensembl, sorted
        //allTranscr = getAllTranscripts();
        getAllTranscripts();
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
        GeneNameConverter converter = null; 
        if (convertGeneName != null)
           converter = new GeneNameConverter(convertGeneName);
        
        int curId = 0;
        System.out.println("Started generating the expression table");
        String trId;
        for (String fName : fNames){
            int rCount=0;
            in = new TextFile(fName, false);
            //getting sample name from the file name
            System.out.println("\nProcessing file: " + fName);
            String sampleId="";
            String[] splName = fName.split("/");
            sampleId = splName[splName.length - 2];
            
            sampleNames[curId] = sampleId;
            
            //get number of mapped reads from filtered sam file
            if (normalize)
                readsPerSample[curId] = getNumMappedReads(fName);
            //in.readLine();
            while ((line = in.readLine()) != null){
                //read the expression from simple txt or from gtf   
                splLine = readExpressionLine(line);
                float count = Float.valueOf(splLine[1]);
                
                //skip not expressed transcripts
                if (count == 0)
                    continue;

                if (convertGeneName != null)
                    trId = converter.getHugoName(splLine[0]);
                else
                    trId = splLine[0];
                
                if (trId != null){
                    //get the transcript position in all transcripts list 
                    index = Collections.binarySearch(allTranscr, trId);

                    if (index >= 0){ //if transcript found
                        featureExpressed[index] = true; //set this feature as expressed
                        rCount+=count; //junk
                        table[index][curId] += count;
                        
                    }   
                }
            }
            curId++;
            System.out.println("overall number of reads=" + rCount);
            System.out.println("overall number of mapped reads=" + readsPerSample[curId - 1]);
            in.close();
            
        }
        //print the result
        if (! normalizeByLength)
            printTable(table);
        else
            printTableNormByLength(table);
    }
    
    
    public static void main(String[] args) throws IOException {
        ProcessTranscriptCounts p;
        //if (args.length > 3){
            //p = new ProcessTranscriptCounts(args[0], args[1], args[2], args[3]);
        p = new ProcessTranscriptCounts("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/mappedData_masked/", 
                "reads_unique_hits.sorted.flux.gtf", 
                "/Users/dashazhernakova/Documents/UMCG/hg19/annotations/annotation_transcr_hg19_geneIds.txt", 
                "/Users/dashazhernakova/Documents/UMCG/data/geuvadis/expression_table/test.txt",
                "reads_unique_hits.sorted.idxstats",
                true);
        //p = new ProcessTranscriptCounts("/Users/dashazhernakova/Documents/UMCG/pickrell_data/tophat_out/", "covBed_exon_counts_out.txt", "/Users/dashazhernakova/Documents/UMCG/hg19/annotation_transcr_hg19.txt", "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/tmp_expression.txt");
                
        p.readCounts();
        //p.readCounts();
       // }
        //else
          //  System.out.println("not enough arguments given");
        
    }
}

