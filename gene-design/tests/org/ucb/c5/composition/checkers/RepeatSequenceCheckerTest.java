package org.ucb.c5.composition.checkers;

import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.*;

public class RepeatSequenceCheckerTest {

    @Test
    /**
     * Test whether RepeatSequenceChecker can correctly detect repetitive
     * sequence
     *
     * @author Kevin Hong
     */
    public void testRepeatSequenceChecker() throws Exception {

        RepeatSequenceChecker RSC = new RepeatSequenceChecker();

        //Create goodlist with repetitive sequences examples, and badlist with repetitive sequences
        List<String> goodlist = new ArrayList<>();
        List<String> badlist = new ArrayList<>();

        //Populate goodlist
        goodlist.add("AGCTAGCTAGCTAaCTtGCT");
        goodlist.add("CTGAttatgcatcAttatTGAtttgcacatt");
        goodlist.add("GTACCTTGacggtacgttcgATGCCAT");
        goodlist.add("GTCctgaactggtctgctactgatcgtaattaaggctctag");
        goodlist.add("AGCTGCATGAGCTGCATG");

        //Populate badlist
        badlist.add("TGTACTGAACTGTTGGAAGTCACACCTATGTACTGAACTGTTGGAAGTCACACCTA");
        badlist.add("AGTACGTACGTTGACAGTACGTACGTTGAC");
        badlist.add("AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT");
        badlist.add("actgctaggctagcactgctaggctagcactgctaggctagc");
        badlist.add("gtagctgctAgtagctgctA");

        //Test and output boolean result
        for (String seq : goodlist) {
            boolean result = RSC.run(seq);
            assertTrue(result);
        }

        for (String seq : badlist) {
            boolean result = RSC.run(seq);
            assertFalse(result);
        }
    }

}
