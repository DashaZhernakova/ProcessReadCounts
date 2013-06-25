
package processreadcounts;

import java.io.IOException;

/**
 *
 * @author dashazhernakova
 */
public class ProcessReadCounts {

    /**
     *
     */
    
    public static void usage(){
        System.out.println("Options:");
        System.out.println("\t--mode\tmakeExpressionTable\n\t\tprocessExonCountsCovBed");
        System.out.println("\t--mode makeExpressionTable:\n\t\t"
                + "--in (Input folder. Should contain folders with alignment files and read counts files for each sample)\n\t\t"
                + "--annot (Annotation file of the format: Platform\tProbeName(as in readCounts file)\tGeneName\tProbeChr\tProbe start\t Probe end)\n\t\t"
                + "--out (Path to output expression table file)\n\t\t"
                + "--pattern (Filename of the file with read counts)\n\t\t"
                + "--alnFnamePattern (Filename of the file with mapped reads)\n\t\t"
                + "--normalize (true if you want the expression values normalized by the total number of mapped reads per sample)\n\t\t"
                + "--convertGeneName (skip it)\n\t\t"
                + "--normalizeByGeneLength (true if you want the expression values normalized by the gene length)\n\t\t"
                + "--fileList (File with file paths to process)");
    }
    public static void main(String[] args) throws IOException {
        
        String mode = null;
        ProcessTranscriptCounts p;
        
        int i = 0;
        for (i = 0; i < args.length; i++) {
	    String arg = args[i];
	    String val = null;

	    if (i + 1 < args.length) {
		val = args[i + 1];
	    }

	    if (arg.equals("--mode")) {
		mode = val;
                //System.out.println("mode");
                break;
	    }
            
	}
        
        if (mode == null) 
	    System.out.println("ERROR: Please supply --mode");
        else if (mode.equals("makeExpressionTable")){
            String annot = null, pattern = null, dir = null, outFile = null, alnPattern = null, convert = null, fileList = null;
            boolean normalize = true, normLen = false, rpkm = false;
            //System.out.println("mode=make_expression_table");
            for (int j = i; j < args.length; j++) {
                String arg = args[j];
                String val = null;

                if (j + 1 < args.length) {
                    val = args[j + 1];
                }

                if (arg.equals("--annot")) {
                    annot = val;
                }
                if (arg.equals("--in")) {
                    dir = val;
                }
                if (arg.equals("--out")) {
                    outFile = val;
                }
                if (arg.equals("--pattern")) {
                    pattern = val;
                }
                if (arg.equals("--alnFnamePattern")) {
                    alnPattern = val;
                }
                if (arg.equals("--fileList")) {
                    fileList = val;
                }
                if (arg.equals("--normalize")) {
                    normalize = Boolean.valueOf(val);
                }
                if (arg.equals("--convertGeneName")) {
                    convert = val;
                }
                if (arg.equals("--normalizeByGeneLength")){
                    normLen = Boolean.valueOf(val);
                }
                if (arg.equals("--rpkm")){
                    rpkm = Boolean.valueOf(val);
                }
                }
            if ( (dir == null) || (pattern == null) || (annot == null) || (outFile == null)){
                System.out.println("Not enough arguments!");
                usage();
            }
            System.out.println("Used settings:\n--in " + dir + "\n--out " + outFile +"\n--annot " + annot + "\n--pattern " + pattern + "\n--normalize " + normalize + "\n--convertGeneName " + convert + "\n--alnFnamePattern " + alnPattern + "\n--normalizeByGeneLength " + normLen + "\n");
            if (rpkm){
                ProcessTranscriptCountsRPKM pr = new ProcessTranscriptCountsRPKM(dir, pattern, annot, outFile, alnPattern);
                pr.readCounts();
            }else{
                if (fileList == null){
                    if (alnPattern == null){
                        if (convert == null){
                            p = new ProcessTranscriptCounts(dir, pattern, annot, outFile, normalize);
                            System.out.println("alnPattern == null, convert == null");
                        }
                        else{
                            p = new ProcessTranscriptCounts(dir, pattern, annot, outFile, normalize, convert);
                            System.out.println("alnPattern == null");
                        }
                    }
                    else{
                        if (convert == null){
                            System.out.println("convert == null");
                            p = new ProcessTranscriptCounts(dir, pattern, annot, outFile, alnPattern, normalize, normLen);
                        }
                        else{
                            p = new ProcessTranscriptCounts(dir, pattern, annot, outFile, alnPattern, normalize, convert);
                            System.out.println("nothing null");
                        }
                    }
                    p.run();
                }
                else{
                    p = new ProcessTranscriptCounts(fileList, annot, outFile);
                }
            }
        }
        
        else if (mode.equals("processExonCountsCovBed")){
            String inDir = "";
            boolean exonwise = false;
            ProcessExonCountsCoverBed pr = new ProcessExonCountsCoverBed();
            
            for (int j = i; j < args.length; j++) {
                String arg = args[j];
                String val = null;
                if (j + 1 < args.length) {
                    val = args[j + 1];
                }
                if (arg.equals("--in"))
                    inDir = val;
                if (arg.equals("exonwise"))
                    exonwise = true;    
            }
            if (exonwise)
                pr.processExonWise(inDir);
            else
                pr.process(inDir);
        }
    }
}
