package org.ucb.c5.composition.checkers;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 * @author Joanne Chang student joanne91218
 */
public class PromoterCheckerTest {
    /**
     * Checks if the set of random sequences avoid the constitutive promoter patches. False if matches.
     *
     * @throws Exception
     */
    @Test
    public void testpromotetchecker() throws Exception {
        PromoterChecker p = new PromoterChecker();
        p.initiate();

        //12 seqs with constitutive promoter motifs
        String[] seq1= {"TTGACAATTAATCATCGAACTAGTATAAT",
                        "TTGACATCTACTGTCGAGAAATTTATAAT",
                        "TTGACATTACTGACTTGAGCGCCTATAAT",
                        "TTGACAGCGAACGCTTCAGATGTTATAAT",
                        "TTGACATCCCAAATTGTAACTTATATAAT",
                        "TTGACAAGAAGGAGTGAATCAATTATAAT",
                        "TTGACATTATGTTTTATTATTATTATAAT",
                        "TTGACAGTATTCTGCACAGGCTCTATAAT",
                        "TTGACATCTTTTACCATGGTCGTTATAAT",
                        "TTGACAAAGTCGATTCCTTCGATTATAAT",
                        "TTGACACCGCTCATTCACTAGGTTATAAT",
                        "TTGACAGGGTGGCTGAACAAGGGTATAAT"};
        
        System.out.println(">>Testing constitutive promoters (false)");
        for(String seq : seq1) {
            boolean result = p.run(seq);
            System.out.println("result: " + result + " on " + seq);
            assert(result == false);
        }
        
        //12 random sequences
        String[] seq2= {"AAACTGTAATCCACCACAAGTCAAGCCAT",
                        "GCCTCTCTGAGACGCCGTATGAATTAATA",
                        "GTAAACTTTGCGCGGGTTCACTGCGATCC",
                        "TTCAGTCTCGTCCAAGGGCACAATCGAAT",
                        "ATCCCCCGAAGTTTAGCAGGTCGTGAGGT",
                        "TCATGGAGGCTCTCGTTCATCCCGTGGGA",
                        "ATCAAGCTTCGCCTTGATAAAGCACCCCG",
                        "TCGGGTGTAGCAGAGAAGACGCCTACTGA",
                        "TTGTGCGATCCCTCCACCTCAGCTAAGGT",
                        "GCTACCAATATTTAGTTTTTTAGCCTTGC",
//                        "ACAGACCTCCTACTTAGATTGCCACGCAT", //Silenced, may have a weak promoter
                        "GTTCCGCTGGGATCCATCGTTGGCGGCCG"};
        
        System.out.println("\n>>Testing random sequences (true)");
        for(String seq : seq2) {
            boolean result = p.run(seq);
            System.out.println("result: " + result + " on " + seq);
            assert(result == true);
        }
        
        //J23119 promoter library
        String[] seq3= {"GAATTCGCGGCCGCTTCTAGAGTTGACGGCTAGCTCAGTCCTAGGTACAGTGCTAGCTAC",
                        "GAATTCGCGGCCGCTTCTAGAGTTTACAGCTAGCTCAGTCCTAGGTATTATGCTAGCTAC",
                        "GAATTCGCGGCCGCTTCTAGAGTTGACAGCTAGCTCAGTCCTAGGTACTGTGCTAGCTAC",
                        "GAATTCGCGGCCGCTTCTAGAGCTGATAGCTAGCTCAGTCCTAGGGATTATGCTAGCTAC",
                        "GAATTCGCGGCCGCTTCTAGAGTTGACAGCTAGCTCAGTCCTAGGTATTGTGCTAGCTAC",
                        "GAATTCGCGGCCGCTTCTAGAGTTTACGGCTAGCTCAGTCCTAGGTACTATGCTAGCTAC",
                        "GAATTCGCGGCCGCTTCTAGAGTTTACGGCTAGCTCAGTCCTAGGTATAGTGCTAGCTAC",
                        "GAATTCGCGGCCGCTTCTAGAGTTTACGGCTAGCTCAGCCCTAGGTATTATGCTAGCTAC",
                        "GAATTCGCGGCCGCTTCTAGAGCTGAcagCTAGCTCAgtcctagGTATAATGCTAGCTAC",
                        "GAATTCGCGGCCGCTTCTAGAGTTTACAGCTAGCTCAGTCCTAGGGACTGTGCTAGCTAC",
                        "GAATTCGCGGCCGCTTCTAGAGTTTACGGCTAGCTCAGTCCTAGGTACAATGCTAGCTAC",
                        "GAATTCGCGGCCGCTTCTAGAGTTGACGGCTAGCTCAgTCCTaGGTATAGTGCTAGCTAC",
                        "GAATTCGCGGCCGCTTCTAGAGCTGATAGCTAGCTCAGTCCTAGGGATTATGCTAGCTAC",
                        "GAATTCGCGGCCGCTTCTAGAGCTGATGGCTAGCTCAGTCCTAGGGATTATGCTAGCTAC",
                        "GAATTCGCGGCCGCTTCTAGAGTTTATGGCTAGCTCAgtcctaggTACAATGCTAGCTAC",
                        "GAATTCGCGGCCGCTTCTAGAGTTTATAGCTAGCTCAGCCCTTgGTACAATGCTAGCTAC",
                        "GAATTCGCGGCCGCTTCTAGAGTTGACAGCTAGCTCAgtcctaggGACTATGCTAGCTAC",
                        "GAATTCGCGGCCGCTTCTAGAGTTGACAGCTAGCTCAGTCCTAGGGATTGTGCTAGCTAC",
                        "GAATTCGCGGCCGCTTCTAGAGTTGACGGCTAGCTCAGTCCTAGGTATTGTGCTAGCTAC",
                        };

        System.out.println("\n>>Testing J23119 promoters (some false)");
        for(String seq : seq3) {
            boolean result = p.run(seq);
            System.out.println("result: " + result + " on " + seq);
//            assert(result == false);
        }

    }
}
