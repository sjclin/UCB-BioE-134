package org.ucb.c5.composition;

import org.ucb.c5.composition.checkers.ForbiddenSequenceChecker;
import org.junit.BeforeClass;
import org.junit.Test;
import org.ucb.c5.composition.model.*;
import org.ucb.c5.sequtils.HairpinCounter;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertFalse;

public class TestJunctions {

    private static CompositionToDNA c2d;
    private static ForbiddenSequenceChecker checker;

    @BeforeClass
    public static void setUpClass() throws Exception {
        c2d = new CompositionToDNA();
        c2d.initiate();
        checker = new ForbiddenSequenceChecker();
        checker.initiate();
    }

    private boolean forbiddenSequenceChecker(String dnaSeq) {
        return (dnaSeq.contains("AAAAAAAA")
                || dnaSeq.contains("TTTTTTTT")
                || dnaSeq.contains("CCCCCCCC")
                || dnaSeq.contains("GGGGGGGG")
                || dnaSeq.contains("ATATATAT")
                || dnaSeq.contains("CAATTG")
                || dnaSeq.contains("GAATTC")
                || dnaSeq.contains("GGATCC")
                || dnaSeq.contains("AGATCT")
                || dnaSeq.contains("ACTAGT")
                || dnaSeq.contains("TCTAGA")
                || dnaSeq.contains("GGTCTC")
                || dnaSeq.contains("CGTCTC")
                || dnaSeq.contains("CACCTGC")
                || dnaSeq.contains("CTGCAG")
                || dnaSeq.contains("CTCGAG")
                || dnaSeq.contains("GCGGCCGC")
                || dnaSeq.contains("AAGCTT"));
    }

    @Test
    /*See if there is any forbidden sequence or secondary structure in between
     * the end of the promoter and the start of RBS
     * @Anruo Shen
     * */
    public void testPromoterRBSjunction() throws Exception {
        //initiate basic information
        String promoter = "ttatgacaacttgacggctacatcattcactttttcttcacaa";
        String terminator = "TGCCTGGCGGCAGTAGCGCGGTGGTCCCACCTGACCCCATGCC";
        String protein = "MYPFIRTARMTVCAKKHVHLTRDAAEQLLADIDRRLDQLLPVE";
        String CheckFBSeq = new String();
        String CheckHP = new String();
        int l = promoter.length();
        List<String> proteins = new ArrayList<>();
        proteins.add(protein);
        Composition comp = new Composition(Host.Ecoli, promoter, proteins, terminator);
        
        //execute the target function
        Construct dna = c2d.run(comp);

        List<Transcript> mRNA = dna.getmRNAs();
        for (Transcript t : mRNA) {
            RBSOption rbs = t.getRbs();

            //get the sequence for checking respectively
            CheckFBSeq = promoter.substring(l - 7) + rbs.getRbs().substring(0, 8);
            CheckHP = promoter + rbs.getRbs();
        }
        
        //check the forbidden sequence
        boolean result1 = forbiddenSequenceChecker(CheckFBSeq);
        
        //check the secondary structure
        boolean result2 = false;
        HairpinCounter hpc = new HairpinCounter();
        hpc.initiate();
        
        double threshold = 5;
        if (hpc.run(CheckHP) > threshold) {
            result2 = true;
        }

        assertFalse(result1 || result2);
    }

    /*See if there is any forbidden sequence or seconday structure
     *in between the end of cds and the start of the terminator
     * @Anruo Shen
     * */
    @Test
    public void testCDSTerminatorjunction() throws Exception {
        //initiate basic information
        String promoter = "GTACCAGTTACGCAGTAGCTAGCATGCTAGCTATGCGTAGTCGATCTAGCTGATCGTAGCTGCTG";
        String terminator = "AGTCGTAGTCGATGCATTTACGATCGCTAGCCCTGTGTGAAAACGTAGCTGTGCTAGTCG";
        String protein = "MKSGKW";
        String CheckFBSeq = new String();
        String CheckHP = new String();
        List<String> proteins = new ArrayList<>();
        proteins.add(protein);
        Composition comp = new Composition(Host.Ecoli, promoter, proteins, terminator);

        //execute the target function
        Construct dna = c2d.run(comp);
        List<Transcript> mRNA = dna.getmRNAs();
        for (Transcript t : mRNA) {
            RBSOption rbs = t.getRbs();
            int l = rbs.getCds().length();

            //get the sequence for checking respectively
            CheckFBSeq = rbs.getCds().substring(l - 7) + terminator.substring(0, 8);
            CheckHP = rbs.getCds() + terminator;
        }

        //check the forbidden sequence
        boolean result1 = forbiddenSequenceChecker(CheckFBSeq);

        //check the secondary structure
        boolean result2 = false;
        HairpinCounter hpc = new HairpinCounter();
        hpc.initiate();
        
        double threshold = 5;
        if (hpc.run(CheckHP) > threshold) {
            result2 = true;
        }
        assertFalse(result1 || result2);
    }
}
