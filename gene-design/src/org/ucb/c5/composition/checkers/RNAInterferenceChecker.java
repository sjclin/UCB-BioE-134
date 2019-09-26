package org.ucb.c5.composition.checkers;

import org.ucb.c5.sequtils.RevComp;
import org.ucb.c5.utils.FileUtils;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

/**
 * Check RNA interference from CDS in the composition
 *
 * @author Kevin Hong
 */
public class RNAInterferenceChecker {

    private RevComp revcomp;

    //Decide the number of base complements to cause interference
    private final int length = 15;

    //Create Hashset containing unique interfering sequences from CDS
    public final Set<String> interferenceSeq = new HashSet<>();

    public void initiate() throws Exception {
        revcomp = new RevComp();
        revcomp.initiate();

        //Create ArrayList to collect all CDS
        ArrayList<String> allSeq = new ArrayList<>();

        //Read RBS and CDS data file
        String CDSData = FileUtils.readResourceFile("composition/data/coli_genes.txt");

        //Split it into lines
        String[] CDSLines = CDSData.split("\\r|\\r?\\n");

        //Scan through each line and collect 5' UTR and CDS sequences
        for (int i = 0; i < CDSLines.length; i++) {
            String line = CDSLines[i];

            //Split the line between tabs
            String[] tabs = line.split("\t");

            //Add to List
            allSeq.add(revcomp.run(tabs[5] + tabs[6]));
        }

        //Collect all unique interference sequences
        for (String seq : allSeq) {

            //Scan through sequence
            for (int i = 0; i < seq.length() - length + 1; i++) {
                String interferingSeq = seq.substring(i, i + length);

                //Add to interference sequence hashset if unique
                if (!interferingSeq.contains(seq)) {
                    interferenceSeq.add(interferingSeq);
                }
            }
        }
    }

    public boolean run(String constructSeq) throws Exception {
        //Uppercase all base pairs
        constructSeq = constructSeq.toUpperCase();

        if (!constructSeq.matches("[ATCG]+")) {
            throw new IllegalArgumentException("RNAInterferenceChecker passed a non [ATCG]+ argument");
        }

        //Ensure sequence is long enough
        if (constructSeq.length() <= length) {
            throw new IllegalArgumentException("RNAInterferenceChecker passed a too-short argument");
        }

        //Create ArrayList for sequences in the construct that may be subject to interference
        Set<String> susceptibleSeq = new HashSet<>();

        //Scan through construct sequence and look for unique sequences of decided length
        for (int i = 0; i < constructSeq.length() - length + 1; i++) {
            String seq = constructSeq.substring(i, i + length);

            //Add sequence to list if unique
            if (!susceptibleSeq.contains(seq)) {
                susceptibleSeq.add(seq);
            }
        }

        //Check if construct sequence may be susceptible for interference
        for (String seq : susceptibleSeq) {

            for (int i = 0; i < seq.length() - length + 1; i++) {
                String seqInQuestion = seq.substring(i, i + length);

                //Check if construct sequence is susceptible for interference
                if (interferenceSeq.contains(seqInQuestion)) {
                    return false;
                }
            }
        }
        return true;
    }

    public static void main(String[] args) throws Exception {

        //Create example construct sequence
        String falseSeq = "TACCAGCTCATAGAGGTCATTCACATCCTGACCATTCAGTTGAGCAAAATAGTTCTTCAGTGCCTGTTTAACCGAGTCACGC";
        String nottrueSeq = "GAATCATGGCCTGCACGGCAATATAATGGACTTCGACATGGCAATAACGCCTCGTTTCTACGTAATAGTATAAACATAAGCAGCCATGCTG";

        //Run and output result
        RNAInterferenceChecker RIC = new RNAInterferenceChecker();
        RIC.initiate();
        boolean result = RIC.run(nottrueSeq);
        System.out.println("RNAInterferenceChecker: " + result);
    }
}
