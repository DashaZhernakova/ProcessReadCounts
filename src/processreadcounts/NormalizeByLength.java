
package processreadcounts;

import java.io.IOException;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author dashazhernakova
 */
public class NormalizeByLength {
    HashMap<String, Integer> lengths;
    public NormalizeByLength(String annotFileName) throws IOException{
        TextFile annot = new TextFile(annotFileName, false);
        String[] els = annot.readLineElems(TextFile.tab);
        while ( (els = annot.readLineElems(TextFile.tab)) != null ){
            lengths.put(els[1], Integer.valueOf(els[5]) - Integer.valueOf(els[4]) + 1);
        }
        
        annot.close();
    }
    public void normalize(String countsFname, String outFname) throws IOException{
        TextFile counts = new TextFile(countsFname, false);
        TextFile out = new TextFile(outFname, true);
        String[] els;
        while ( (els = counts.readLineElems(TextFile.tab)) != null ){
            out.writeln(els[0] +"\t" + Integer.toString(Integer.parseInt(els[1])/lengths.get(els[0])));
        }
        counts.close();
        out.close();
    }
    public static void main(String[] args) throws IOException {
        NormalizeByLength n = new NormalizeByLength("/Users/dashazhernakova/Documents/UMCG/hg19/annotations/annotation_genes_hg19_v69.txt");
        n.normalize("/Users/dashazhernakova/Documents/UMCG/globinReduction/2476-700/geneCounts_sorted.txt",
                "/Users/dashazhernakova/Documents/UMCG/globinReduction/2476-700/geneCounts_normByLen.txt");
    }
}
