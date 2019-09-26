package org.ucb.c5.composition.checkers;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 * @author Joanne Chang student joanne91218
 */
public class ForbiddenSequenceCheckerTest {

    /**
     * CheckerTest for ForbiddenSequenceChecker.  It gives forbidden sequences
     * that should be caught, and random sequences that should not.
     *
     * @throws Exception
     */
    @Test
    public void testForbiddenSites() throws Exception {
        ForbiddenSequenceChecker checker = new ForbiddenSequenceChecker();
        checker.initiate();

        //12 seqs with forbidden sequences
        String[] seq1= {"TTGACAATTgaattcCGAACTAGTATAAT",
                        "TTTTTTTTTTGTCGAGAAATTTATAAT",
                        "TTGACATTACCGTCTCGAGCGCCTATAAT",
                        "TTGACAGCGAACGCTTCAGACTGCAGAT",
                        "TTGACTGCAGTTGTAACTTATATAAT",
                        "TTGACAAGAAGCGGCCGCTCAATTATAAT",
                        "TTGACATTATGACACCTGCTTATTATAAT",
                        "TTGACACTCGAGCACAGGCTCTATAAT",
                        "TTGACATCGGGGGGGGTTTTACCATGGTCGTTATAAT",
                        "TTGACAAAGTCGATTTTTTTTTTCCTTCGATTATAAT",
                        "TTGACACGGTCTCATTCACTAGGTTATAAT",
                        "TTGACAGTCTAGAGTCTGAACAAGGAGATCTTAAT"};
        
        System.out.println(">>Testing problematic sequences (false)");
        for(String seq : seq1) {
            boolean result = checker.run(seq);
            System.out.println("result: " + result + " on " + seq);
            assert(result == false);
        }
        
        //12 random sequences
        String[] seq2= {"AAACTGTAATCCACCACAAGTCAAGCCAT",
                        "GCCTCTCTGAGGACGCCGTATGAATTAATA",
                        "GTAAACTTTGCGCGGGTTCACTGCGATCC",
                        "TTCAGTCTCGTCCAAGGGCACAATCGAAT",
                        "ATCCCCCGAAGTTTAGCAGGTCGTGAGGT",
                        "TCATGGAGGCTCTCGTTCATCCCGTGGGA",
                        "ATCAAGGCTTCGCCTTGATAAAGCACCCCG",
                        "TCGGGTGTAGCAGAGAAGACGCCTACTGA",
                        "TTGTGCGATCCCTCCACCTCAGCTAAGGT",
                        "GCTACCAATATTTAGTTTTTTAGCCTTGC",
                        "ACAGACCTCCTACTTAGATTGCCACGCAT",
                        "GTTCCGCTGGCGATCCATCGTTGGCGGCCG"};
        
        System.out.println("\n>>Testing random sequences (true)");
        for(String seq : seq2) {
            boolean result = checker.run(seq);
            System.out.println("result: " + result + " on " + seq);
            assert(result == true);
        }
    }
}
