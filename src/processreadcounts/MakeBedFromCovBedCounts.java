package processreadcounts;


import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;
import processtmap.GeneNameConverter;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author dashazhernakova
 */
public class MakeBedFromCovBedCounts {
    Map<String, String> id2pos;
    public void readAnnotation(String fName) throws IOException{
        TextFile annot = new TextFile(fName, false);
        id2pos = new TreeMap<String, String>();
        String[] els = annot.readLineElems(TextFile.tab);
        while ( (els = annot.readLineElems(TextFile.tab)) != null ){
            id2pos.put(els[1], els[3] + "\t" + els[4] + "\t" + els[5]);
        }
        annot.close();
    }
    
    public void convert(String fName) throws IOException{
        TextFile in = new TextFile(fName, false);
        TextFile out = new TextFile(fName.replace(".txt", ".bed"), true);
        GeneNameConverter converter = new GeneNameConverter();
        converter.makeIdToHUGOmap();
        String[] spl;
        String hugo;
        while( (spl = in.readLineElems(TextFile.tab)) != null ){
            hugo = converter.getHugoName(spl[0]);
            if (hugo != null)
            if (id2pos.containsKey(hugo))
                out.writeln(id2pos.get(hugo)+ "\t" + hugo + "\t" + spl[1]);
        }
        in.close();
        out.close();
        
    }
    public static void main(String[] args) throws IOException {
        MakeBedFromCovBedCounts m = new MakeBedFromCovBedCounts();
        m.readAnnotation("/Users/dashazhernakova/Documents/UMCG/hg19/annotation_genes_hg19.txt");
        m.convert("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/simulation/6_NA18486_yale/cut_sample.genes.results.txt");
    }
}
