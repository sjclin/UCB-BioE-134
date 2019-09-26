package org.ucb.c5.composition.checkers;

import org.ucb.c5.sequtils.RevComp;
import java.util.ArrayList;
import java.util.List;

/**
 * Checks to see if any Internal Ribosome Binding Sites existt within a string
 * By Naveen Kumaran
 */
public class InternalRBSChecker {

    private List<String> sdSequences;
    private RevComp revcomp;

    public void initiate() throws Exception {

        // create a list of RBSOptions to search for in desired code.
        // Found from a random internet website (idk how credible the site was).
        // Any RBS sequences found later can be added below.
        sdSequences = new ArrayList<>();
        sdSequences.add("AGGAGG");
        sdSequences.add("AGCAGG");
        sdSequences.add("AGGA");
        sdSequences.add("GGAG");
        sdSequences.add("AAGGAGGT");
        sdSequences.add("GGAGGTG");
        sdSequences.add("AAGGAGGG");
        sdSequences.add("TAAGGAGG");
        sdSequences.add("GGGTTGAT");
        sdSequences.add("GTAAGGTGA");

        revcomp = new RevComp();
        revcomp.initiate();
    }

    public boolean run(String inSeq) throws Exception {
        //Uppercase and include reverse complement
        String upseq = inSeq.toUpperCase();
        String revseq = revcomp.run(upseq);
        String combiSeq = upseq + "aaaaaaaaaaaaaaaaaaaaaaaaaaaaa" + revseq;

        ArrayList<Integer> startCodonPositions = new ArrayList<>();

        //Run through sequence and look for start start codons.
        //note the indexes of the start codon into the arraylist.
        for (int i = 0; i < combiSeq.length() - 3; i++) {
            String codon = combiSeq.substring(i, i + 3);
            if (codon.equals("ATG")) {
                startCodonPositions.add(i);
            }
        }

        for (int Position : startCodonPositions) {
            String utr = combiSeq.substring(Position - 12, Position - 4);
            for (String sd : sdSequences) {
                if (utr.contains(sd)) {
                    return false;
                }

                //JCA:  Silencing this block, the logic ix not correct.  The RC needs to be handled differently (also added)
//                if (CurrentSeq.contains(revcomp.run(sd))) {
//                    return false;
//                }
            }
        }

        return true;
    }

    public static void main(String[] args) throws Exception {
        InternalRBSChecker checker = new InternalRBSChecker();
        checker.initiate();

        //Strong rbs
        String seq1 = "tttttttgggctaacAGGAGGaattaaccATG";
        boolean result = checker.run(seq1);
        System.out.println("Strong (false):     " + result + " on " + seq1);

        //Strong rbs, wrong spacing
        String seq2 = "tttttttgggctaacAGGAGGaATG";
        boolean result2 = checker.run(seq2);
        System.out.println("Bad spacing (true): " + result2 + "  on " + seq2);

        //Strong rbs, wrong spacing
        String seq3 = "tttttttgggctaacAGGAGGaaaaaaattaaccATG";
        boolean result3 = checker.run(seq3);
        System.out.println("Bad spacing (true): " + result3 + "  on " + seq3);
    }
}
