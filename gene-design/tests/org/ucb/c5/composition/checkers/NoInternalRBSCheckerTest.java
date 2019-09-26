/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.ucb.c5.composition.checkers;
import org.junit.BeforeClass;
import org.junit.Test;
/*Program checks if there are no internal rbs.
 *if there are internal RBS, returns false.
 *if there aren't internal RBS, returns true.
 * t's are used in sequences in order to make sure that there are no false positives or true negatives with respect to the assertions. 
 * @author Niki Shakouri
 */
public class NoInternalRBSCheckerTest {
    
    private static NoInternalRBSChecker checker;

    @BeforeClass
    public static void setUpClass() throws Exception {

        checker = new NoInternalRBSChecker();
        checker.initiate();
    }

    @Test
    public void IsDNASequence() throws Exception {
        //Input sequence is a DNA sequence (has only A,C,T,G as characters), should return True (since this means the sequence is invalid and can't have RBS).
        String DNA1 ="tttttttttttttMtttttttttttttATG";
        String DNA2 ="MtttttttttttttttttttttttttttATG";
        String DNA3 ="tttttttttttttttttttttttttttATGM";
        
        try {
            checker.run(DNA1);
            assert(false);
        }  catch(Exception err) {}
        
        try {
            checker.run(DNA2);
            assert(false);
        }  catch(Exception err) {}
        
        try {
            checker.run(DNA3);
            assert(false);
        }  catch(Exception err) {}
        
        assert(true);
        
    }
    
    @Test
    public void NoInternalRBS() throws Exception {
        //Input sequence doesn't have internal rbs, should return True.
        //These DNA sequences are from E. Coli File
        //String Ecoli1 = "ATGTTCGAACAACGCGTAAATTCTGACGTACTGACCGTTTCTACCGTTAACTCTCAGGATCAGGTAACCCAAAAACCCCTGCGTGACTCGGTTAAACAGGCACTGAAGAACTATTTTGCTCAACTGAATGGTCAGGATGTGAATGACCTCTATGAGCTGGTACTGGCTGAAGTAGAACAGCCCCTGTTGGACATGGTGATGCAATACACCCGTGGTAACCAGACCCGTGCTGCGCTGATGATGGGCATCAACCGTGGTACGCTGCGTAAAAAATTGAAAAAATACGGCATGAACTAA";
        //String Ecoli2 = "ATGAGCGAAGCACTTAAAATTCTGAACAACATCCGTACTCTTCGTGCGCAGGCAAGAGAATGTACACTTGAAACGCTGGAAGAAATGCTGGAAAAATTAGAAGTTGTTGTTAACGAACGTCGCGAAGAAGAAAGCGCGGCTGCTGCTGAAGTTGAAGAGCGCACTCGTAAACTGCAGCAATATCGCGAAATGCTGATCGCTGACGGTATTGACCCGAACGAACTGCTGAATAGCCTTGCTGCCGTTAAATCTGGCACCAAAGCTAAACGTGCTCAGCGTCCGGCAAAATATAGCTACGTTGACGAAAACGGCGAAACTAAAACCTGGACTGGCCAAGGCCGTACTCCAGCTGTAATCAAAAAAGCAATGGATGAGCAAGGTAAATCCCTCGACGATTTCCTGATCAAGCAATAA";
        //boolean result1 = checker.run(Ecoli1);
        //boolean result2 = checker.run(Ecoli2);
        //assert(result1 == true);
        //assert(result2 == true);
    }
    
    @Test
    public void HasInternalRBS() throws Exception {
        //Input sequence has internal rbs, should returns False.
        String BBa_B0029 = "TCTAGAGTTCACACAGGAAACCTACTAGATG";
        String BBa_B0034 = "TCTAGAGAAAGAGGAGAAATACTAGATG";
        //String BBa_B0034TTG = "TCTAGAGAAAGAGGAGAAATACTAGTTG";
        String BBa_B0033 = "TCTAGAGTCACACAGGACTACTAGATG";
        String BBa_B0031 = "TCTAGAGTCACACAGGAAACCTACTAGATG";
        boolean result1 = checker.run(BBa_B0029);
        boolean result2 = checker.run(BBa_B0034);
        //boolean result3 = checker.run(BBa_B0034TTG);
        boolean result4 = checker.run(BBa_B0033);
        boolean result5 = checker.run(BBa_B0031);
        assert(result1 == false);
        assert(result2 == false);
        //assert(result3 == false);
        assert(result4 == false);
        assert(result5 == false);
        
    }
          
  
    @Test
    public void ContainsInternalRBS() throws Exception {
        //Input sequence has internal rbs S after the start codon, not acceptable locations, should return True (meaning it does not find RBS).
        //Should return TRUE
        String DNA = "tttttttttttttATGtttttAAGGAGGGtttttttt";
        boolean result = checker.run(DNA);
        assert(result == true);
    }

    @Test
    public void InternalRBSBetweenStartCodons() throws Exception {
        //Input sequence has an internal RBS between two start codons, can be traced, should return False.
        String DNA = "tttttttttttttATGtttttAAGGAGGGtttATGtttttttt";
        boolean result = checker.run(DNA);
        assert(result == false);
    }

    @Test
    public void InternalRBSBeforeAndAfterStartCodon() throws Exception {
        //Input sequence has an internal RBS before and after start codon, but since it has one before, should return False.
        String DNA = "ttttttAGGAGGtttttttATGtttttAAGGAGGGttttttttttt";
        boolean result = checker.run(DNA);
        assert(result == false);
    }

    @Test
    public void InternalRBSNextToStartCodon() throws Exception {
        //Input sequence is an internal RBS that is right next to start codon (can't be found) should return True.
        String DNA = "ttttttttttttttttttAGGAGGATGttttttttttt";
        boolean result = checker.run(DNA);
        assert(result == true);
    
}
    @Test
    public void InternalRBSFarFromStartCodon() throws Exception {
        //Input sequence is an internal RBS that is far from start codon (can't be found) should return True.
        String DNA = "AGGAGGttttttttttttttATGttttttttttt";
        boolean result = checker.run(DNA);
        assert(result == true);
    
    }
    @Test
    public void RandomSequences() throws Exception {
    // Input sequence is random, has no internal RBS, should return True.
        String Random1 = "ATTTTTTTGGGGGGGCCCCCGGGGATG";
        String Random2 = "aaaaaatttttGCCCAATTTTGAAAaaaaATG";
        //String Random3 = "ATGATGATGATGATGATGATGATG";
        boolean result1 = checker.run(Random1);
        boolean result2 = checker.run(Random2);
        //boolean result3 = checker.run(Random3);
        assert(result1 == true);
        assert(result2 == true);
        //assert(result3 == true);
        

    }
    @Test
    public void ATGIncorrectStartLocation() throws Exception {
    // Input sequence does not have ATG in the correct location for a sliding frame to be computed (meaning sequence does not have atleast 14 base pairs to left of the start codon), should return True. 
    String DNA1 = "ATGttttAGGAGGttttttttttttt";
    String DNA2 = "AGGAGGtttttATGtttttttttttttt";
    boolean result1 = checker.run(DNA1);
    boolean result2 = checker.run(DNA2);
    assert(result1 == true);
    assert(result2 == true);
        
    }
    @Test
    public void NoATG() throws Exception {
    // Input sequence has no ATG,therfore can't have internal RBS, should return True.
    String DNA1 = "ttttAGGAGGttttttttttttt";
    String DNA2 = "AGGAGGAttttttttttttttttTGtt";
    boolean result1 = checker.run(DNA1);
    boolean result2 = checker.run(DNA2);
    assert(result1 == true);
    assert(result2 == true);
        
    }
}
