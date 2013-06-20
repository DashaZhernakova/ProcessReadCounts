/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package processreadcounts;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author dashazhernakova
 */
public class CountReadsFromPairToBed {
    /**
     * reads output of pairToBed, determines the transcripts the reads can be assigned to 
     * @throws IOException 
     */
    public void count() throws IOException{
        TextFile tf = new TextFile("/Users/dashazhernakova/Documents/UMCG/hapmapData/tophat/1382_1_1/6.pairToBed_out.bed", false);
        tf.open();
        String[] fields;
        Map<String,ArrayList<String>> read_map = new HashMap<String, ArrayList<String>>();
        //first read
        fields = tf.readLineElems(TextFile.tab);
        String cur_read = fields[6];
        read_map.put(cur_read, new ArrayList<String>());
        ArrayList<String> cur_trList = read_map.get(cur_read);
        cur_trList.add(fields[13].replaceAll(":exon:[0-9]*", ""));
                            
        while ( (fields = tf.readLineElems(TextFile.tab)) != null){
            if (cur_read.equals("IL26_1382:1:12:857:252"))
                System.out.println("read");
            if (cur_read.equals(fields[6]))
                cur_trList.add(fields[13].replaceAll(":exon:[0-9]*", ""));
            else{
                //cur_trList = getCandidateTranscripts(cur_trList);
                if (getCandidateTranscripts(cur_trList) > 0)
                    System.out.println(" > 1: " + cur_read + " : " + fields[13]);
                cur_read = fields[6];
                if (cur_read.equals("IL26_1382:1:12:857:252"))
                    System.out.println("read");
                read_map.put(cur_read, new ArrayList<String>());
                cur_trList = read_map.get(cur_read);
                cur_trList.add(fields[13].replaceAll(":exon:[0-9]*", ""));
            }
                
            
        }
        
        tf.close();
    }

    //private ArrayList<String> getCandidateTranscripts(ArrayList<String> cur_trList) {
    private int getCandidateTranscripts(ArrayList<String> cur_trList) {
        ArrayList<String> new_trList = new ArrayList<String>();
        Collections.sort(cur_trList);
        String cur_tr = cur_trList.get(0);
        int cur_count = 0;
        String tr;
        for (int i = 1 ; i < cur_trList.size(); i++){
            tr = cur_trList.get(i);
            if (tr.equals(cur_tr))
                cur_count ++;
            else{
                if (cur_count > 0)
                    new_trList.add(cur_tr);
                if (cur_count > 1){
                    System.out.println(" > 1: " + cur_tr);
                    return cur_count;
                    
                }
                cur_count = 0;
                cur_tr = tr;
            }
        }
        if (cur_count > 0)
            new_trList.add(cur_tr);
        return 0;
        //return new_trList;
            
    }
    public static void main(String[] args) throws IOException {
        CountReadsFromPairToBed c = new CountReadsFromPairToBed();
        c.count();
    }
}
