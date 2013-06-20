
package processreadcounts;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author dashazhernakova
 */
public class MakeRPKMsfromRawCountsTable {
    HashMap<String, Integer> geneLengths;
    HashMap<String, Integer> numMapped;
    HashMap<Integer, String> posToSampleId;
    
    public void makeHashes(String bed12Name, String flagstatPattern, String dir) throws IOException{
        getLengths(bed12Name);
        getNumMappedFlagstat(flagstatPattern, dir);
        
    }
    /*
     * Annotation file in bed12 should be sorted by gene name!!!
     */
    private void getLengths(String bed12Name) throws IOException{
        geneLengths = new HashMap<String, Integer>();
        TextFile bed12 = new TextFile(bed12Name, false);
        String[] els = bed12.readLineElems(TextFile.tab);
        String cur_gene = els[3], lengths = els[10], gene;
        int gene_len;
        while( (els = bed12.readLineElems(TextFile.tab)) != null){
            gene = els[3];
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
    
    private void getNumMappedFlagstat(String flagstatPattern, String dirName) throws IOException{
        numMapped = new HashMap<String, Integer>();
        TextFile flagstat;
        File dir = new File(dirName);
        for (File ch : dir.listFiles()) {
            if (ch.isDirectory()){
                for (File child : ch.listFiles()){
                   if ( (child.getName().equals(flagstatPattern))){
                       flagstat = new TextFile(child.getPath(), false);
                       flagstat.readLine();
                       flagstat.readLine();
                       int num = Integer.parseInt(flagstat.readLineElems(TextFile.space)[0]);
                       String sampleId = ch.getName();
                       numMapped.put(sampleId, num);
                       flagstat.close();
                   }
                }
            }
        }
    }
    
    private void getPosToSampleId(String[] header){
        posToSampleId = new HashMap<Integer, String>();
        for (int i = 1; i < header.length; i++){
            posToSampleId.put(i, header[i]);
        }
    }
    public void changeTable(String inFname, String outFname) throws IOException{
        TextFile in = new TextFile(inFname, false);
        TextFile out = new TextFile(outFname, true);
        
        String[] els = in.readLineElems(TextFile.tab);
        getPosToSampleId(els);
        out.writelnTabDelimited(els);
        
        String gene;
        float rpkm;
        int len;
        while ( (els = in.readLineElems(TextFile.tab)) != null ){
            gene = els[0];
            len = geneLengths.get(gene);
            out.write(gene);
            for (int i = 1; i < els.length; i++){
                rpkm = (1000000000*Float.parseFloat(els[i]))/(numMapped.get(posToSampleId.get(i))*(long)len);
                out.write("\t" + rpkm);
            }
            out.writeln();
        }
        in.close();
        out.close();
    }
    public static void main(String[] args) throws IOException {
        MakeRPKMsfromRawCountsTable m = new MakeRPKMsfromRawCountsTable();
        m.makeHashes("/Users/dashazhernakova/Documents/UMCG/hg19/exonic_genes_v69_stranded_cut.bed12", 
                "accepted_hits.filtered.flagstat", 
                "/Users/dashazhernakova/Documents/UMCG/data/lincRNA_Sebo/mappedData/");
        m.changeTable("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/tmp.txt", 
                "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/tmp2.txt");
    }
}
