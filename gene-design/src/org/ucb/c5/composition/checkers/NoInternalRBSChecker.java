/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.ucb.c5.composition.checkers;

import org.ucb.c5.sequtils.RevComp;
import java.util.ArrayList;

/**
 * Checks a sequence for internal RBS
 *
 * Uses Position Weight Matrix (#7) from source below:
 * http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.221.9530&rep=rep1&type=pdf
 *
 * @author Niki Shakouri
 */
public class NoInternalRBSChecker {

    private double[][] pwm;

    public void initiate() {

// PWM matrix that depicts the sequence pattern of internal ribosome binding sites
        pwm = new double[][]{
            {0.38, 0.55, 0.09, 0.10, 0.60, 0.27},
            {0.31, 0.15, 0.07, 0.07, 0.10, 0.07},
            {0.07, 0.17, 0.11, 0.05, 0.14, 0.16},
            {0.24, 0.14, 0.73, 0.75, 0.13, 0.46}
        };

    }

    public boolean run(String seq) {
        seq = seq.toUpperCase();
        if (seq.matches("[ATCG]+")) {
        } else {
            throw new IllegalArgumentException("Must contain only A,C,T, or G");
        }

        //Computing score of each sliding frame (to determine if it has an internal RBS)
        int slidingframe = 6;
        double score = 0.00;
        int pwm_ncol = seq.length();

        //JCA:  Adjusted threshold
        double threshold = 2.1;

        ArrayList<Integer> startCodonPositions = new ArrayList<>();

        //Runs through sequence,finding and recording start codon (ATG)
        for (int i = 0; i <= seq.length() - 3; i++) {
            String codon = seq.substring(i, i + 3);
            if (codon.equals("ATG")) {
                startCodonPositions.add(i);
            }
        }
        //RBS is 8-14 base pairs away from start codon
        for (int p : startCodonPositions) {
            if (p - 14 >= 0) {
                for (int start = p - 14; start < p - 8; start++) {
                    score = 0.00;

                    String partseq = seq.substring(start, start + slidingframe);

                    for (int x = 0; x < partseq.length(); x++) {
                        int y = -1;
                        switch (partseq.charAt(x)) {
                            case 'A':
                                y = 0;
                                break;
                            case 'C':
                                y = 1;
                                break;
                            case 'T':
                                y = 2;
                                break;
                            case 'G':
                                y = 3;
                                break;
                        }
                        if (y == -1) {
                            continue;
                        }

                        score += pwm[y][x];

                    }
                //If the score is greater or equal to the threshold value, then you've found an internal RBS, should return False

                    if (score >= threshold) {
                        return false;
                    } else {
                        continue;

                    }
                }
            }
        }
        return true;
    }

    public static void main(String[] args) {
        NoInternalRBSChecker checker = new NoInternalRBSChecker();
        checker.initiate();

        //False = has RBS; True = doesn't have RBS
        //Strong RBS bins, should return false 
        String BBa_B0029 = "TCTAGAGTTCACACAGGAAACCTACTAGATG";
        String BBa_B0034 = "TCTAGAGAAAGAGGAGAAATACTAGATG";
        String BBa_B0033 = "TCTAGAGTCACACAGGACTACTAGATG";
        String BBa_B0031 = "TCTAGAGTCACACAGGAAACCTACTAGATG";
        String BBa_J61101 = "TCTAGAGAAAGACAGGACCCACTAGATG";
        String BBa_J61127 = "TCTAGAGAAAGAGTGGAACTACTAGATG";
        String BBa_J61100 = "TCTAGAGAAAGAGGGGACAAACTAGATG";
        String BBa_J61117 = "TCTAGAGAAAGACATGAGTTACTAGATG";
        String B2 = "TCTAGAGAAAGATCCGATGTACTAGATGAAACCCAAATAATACTAGT";

        //From E. Coli File, do not have RBS, should return true
        String Ecoli1 = "TGTTCGAACAACGCGTAAATTCTGACGTACTGACCGTTTCTACCGTTAACTCTCAGGATCAGGTAACCCAAAAACCCCTGCGTGACTCGGTTAAACAGGCACTGAAGAACTATTTTGCTCAACTGAATGGTCAGGATGTGAATGACCTCTATGAGCTGGTACTGGCTGAAGTAGAACAGCCCCTGTTGGACATGGTGATGCAATACACCCGTGGTAACCAGACCCGTGCTGCGCTGATGATGGGCATCAACCGTGGTACGCTGCGTAAAAAATTGAAAAAATACGGCATGAACTAA";
        String Ecoli2 = "TGAGCGAAGCACTTAAAATTCTGAACAACATCCGTACTCTTCGTGCGCAGGCAAGAGAATGTACACTTGAAACGCTGGAAGAAATGCTGGAAAAATTAGAAGTTGTTGTTAACGAACGTCGCGAAGAAGAAAGCGCGGCTGCTGCTGAAGTTGAAGAGCGCACTCGTAAACTGCAGCAATATCGCGAAATGCTGATCGCTGACGGTATTGACCCGAACGAACTGCTGAATAGCCTTGCTGCCGTTAAATCTGGCACCAAAGCTAAACGTGCTCAGCGTCCGGCAAAATATAGCTACGTTGACGAAAACGGCGAAACTAAAACCTGGACTGGCCAAGGCCGTACTCCAGCTGTAATCAAAAAAGCAATGGATGAGCAAGGTAAATCCCTCGACGATTTCCTGATCAAGCAATAA";
        String Ecoli3 = "ATGAGTAGTAAAGAACAGAAAACGCCTGAGGGGCAAGCCCCGGAAGAAATTATCATGGATCAGCACGAAGAGATTGAGGCAGTTGAGCCAGAAGCTTCTGCTGAGCAGGTGGATCCGCGCGATGAAAAAGTTGCGAATCTCGAAGCTCAGCTGGCTGAAGCCCAGACCCGTGAACGTGACGGCATTTTGCGTGTAAAAGCCGAAATGGAAAACCTGCGTCGTCGTACTGAACTGGATATTGAAAAAGCCCACAAATTCGCGCTGGAGAAATTCATCAACGAATTGCTGCCGGTGATTGATAGCCTGGATCGTGCGCTGGAAGTGGCTGATAAAGCTAACCCGGATATGTCTGCGATGGTTGAAGGCATTGAGCTGACGCTGAAGTCGATGCTGGATGTTGTGCGTAAGTTTGGCGTTGAAGTGATCGCCGAAACTAACGTCCCACTGGACCCGAATGTGCATCAGGCCATCGCAATGGTGGAATCTGATGACGTTGCGCCAGGTAACGTACTGGGCATTATGCAGAAGGGTTATACGCTGAATGGTCGTACGATTCGTGCGGCGATGGTTACTGTAGCGAAAGCAAAAGCTTAA";
        String Ecoli4 = "ATGACGGACAAATTGACCTCCCTTCGTCAGTACACCACCGTAGTGGCCGACACTGGGGACATCGCGGCAATGAAGCTGTATCAACCGCAGGATGCCACAACCAACCCTTCTCTCATTCTTAACGCAGCGCAGATTCCGGAATACCGTAAGTTGATTGATGATGCTGTCGCCTGGGCGAAACAGCAGAGCAACGATCGCGCGCAGCAGATCGTGGACGCGACCGACAAACTGGCAGTAAATATTGGTCTGGAAATCCTGAAACTGGTTCCGGGCCGTATCTCAACTGAAGTTGATGCGCGTCTTTCCTATGACACCGAAGCGTCAATTGCGAAAGCAAAACGCCTGATCAAACTCTACAACGATGCTGGTATTAGCAACGATCGTATTCTGATCAAACTGGCTTCTACCTGGCAGGGTATCCGTGCTGCAGAACAGCTGGAAAAAGAAGGCATCAACTGTAACCTGACCCTGCTGTTCTCCTTCGCTCAGGCTCGTGCTTGTGCGGAAGCGGGCGTGTTCCTGATCTCGCCGTTTGTTGGCCGTATTCTTGACTGGTACAAAGCGAATACCGATAAGAAAGAGTACGCTCCGGCAGAAGATCCGGGCGTGGTTTCTGTATCTGAAATCTACCAGTACTACAAAGAGCACGGTTATGAAACCGTGGTTATGGGCGCAAGCTTCCGTAACATCGGCGAAATTCTGGAACTGGCAGGCTGCGACCGTCTGACCATCGCACCGGCACTGCTGAAAGAGCTGGCGGAGAGCGAAGGGGCTATCGAACGTAAACTGTCTTACACCGGCGAAGTGAAAGCGCGTCCGGCGCGTATCACTGAGTCCGAGTTCCTGTGGCAGCACAACCAGGATCCAATGGCAGTAGATAAACTGGCGGAAGGTATCCGTAAGTTTGCTATTGACCAGGAAAAACTGGAAAAAATGATCGGCGATCTGCTGTAA";
        String Ecoli5 = "ATGACTGAATCTTTTGCTCAACTCTTTGAAGAGTCCTTAAAAGAAATCGAAACCCGCCCGGGTTCTATCGTTCGTGGCGTTGTTGTTGCTATCGACAAAGACGTAGTACTGGTTGACGCTGGTCTGAAATCTGAGTCCGCCATCCCGGCTGAGCAGTTCAAAAACGCCCAGGGCGAGCTGGAAATCCAGGTAGGTGACGAAGTTGACGTTGCTCTGGACGCAGTAGAAGACGGCTTCGGTGAAACTCTGCTGTCCCGTGAGAAAGCTAAACGTCACGAAGCCTGGATCACGCTGGAAAAAGCTTACGAAGATGCTGAAACTGTTACCGGTGTTATCAACGGCAAAGTTAAGGGCGGCTTCACTGTTGAGCTGAACGGTATTCGTGCGTTCCTGCCAGGTTCTCTGGTAGACGTTCGTCCGGTGCGTGACACTCTGCACCTGGAAGGCAAAGAGCTTGAATTTAAAGTAATCAAGCTGGATCAGAAGCGCAACAACGTTGTTGTTTCTCGTCGTGCCGTTATCGAATCCGAAAACAGCGCAGAGCGCGATCAGCTGCTGGAAAACCTGCAGGAAGGCATGGAAGTTAAAGGTATCGTTAAGAACCTCACTGACTACGGTGCATTCGTTGATCTGGGCGGCGTTGACGGCCTGCTGCACATCACTGACATGGCCTGGAAACGCGTTAAGCATCCGAGCGAAATCGTCAACGTGGGCGACGAAATCACTGTTAAAGTGCTGAAGTTCGACCGCGAACGTACCCGTGTATCCCTGGGCCTGAAACAGCTGGGCGAAGATCCGTGGGTAGCTATCGCTAAACGTTATCCGGAAGGTACCAAACTGACTGGTCGCGTGACCAACCTGACCGACTACGGCTGCTTCGTTGAAATCGAAGAAGGCGTTGAAGGCCTGGTACACGTTTCCGAAATGGACTGGACCAACAAAAACATCCACCCGTCCAAAGTTGTTAACGTTGGCGATGTAGTGGAAGTTATGGTTCTGGATATCGACGAAGAACGTCGTCGTATCTCCCTGGGTCTGAAACAGTGCAAAGCTAACCCGTGGCAGCAGTTCGCGGAAACCCACAACAAGGGCGACCGTGTTGAAGGTAAAATCAAGTCTATCACTGACTTCGGTATCTTCATCGGCTTGGACGGCGGCATCGACGGCCTGGTTCACCTGTCTGACATCTCCTGGAACGTTGCAGGCGAAGAAGCAGTTCGTGAATACAAAAAAGGCGACGAAATCGCTGCAGTTGTTCTGCAGGTTGACGCAGAACGTGAACGTATCTCCCTGGGCGTTAAACAGCTCGCAGAAGATCCGTTCAACAACTGGGTTGCTCTGAACAAGAAAGGCGCTATCGTAACCGGTAAAGTAACTGCAGTTGACGCTAAAGGCGCAACCGTAGAACTGGCTGACGGCGTTGAAGGTTACCTGCGTGCTTCTGAAGCATCCCGTGACCGCGTTGAAGACGCTACCCTGGTTCTGAGCGTTGGCGACGAAGTTGAAGCTAAATTCACCGGCGTTGATCGTAAAAACCGCGCAATCAGCCTGTCTGTTCGTGCGAAAGACGAAGCTGACGAGAAAGATGCAATCGCAACTGTTAACAAACAGGAAGATGCAAACTTCTCCAACAACGCAATGGCTGAAGCTTTCAAAGCAGCTAAAGGCGAGTAA";

        //Random Sequence, has no RBS, should return True.
        String Random = "ATTTTTTTGGGGGGGCCCCCGGGGATG";

        //Random Sequence with Shine-Delgarno Sequence, has RBS, should return False.
        String Random1 = "aaaaaatttttAGGAGGaaaaaATG";

        //5' UTR, have RBS, should return False
        String B3261 = "GCACATTCAACGCCATTGAGGATGCCAGCGAACAGCTGGAGGCGTTGGAGGCATACTTCGAAAATTTTGCGTAAACAGAAATAAAGAGCTGACAGAACTATG";
        String B1237 = "CAATTTTGAATTCCTTACATTCCTGGCTATTGCACAACTGAATTTAAGGCTCTATTATTACCTCAACAAACCACCCCAATATAAGTTTGAGATTACTACAATG";
        String B0907 = "TCACTGAATGATAAAACCGATAGCCACAGGAATAATGTATTACCTGTGGTCGCAATCGATTGACCGCGGGTTAATAGCAACGCAACGTGGTGAGGGGAAATG";
        String B2518 = "TGTTACCTTCATCAATAGTCAACGGCCCTGTTGCTCATTATAATCCGCGCCATCTCGTACGCTGGTACAGACAACAACAGAACAATTTACAGAGGTAAAAATG";

        //rpsA --> should completely lack shine-delgarno sequence (AGGAGG) http://emboj.embopress.org/content/20/15/4222
        //Excel Spreadsheet of RBS 
        //Strong bins: have RBS, should return false
        //String B3 = "TCTAGAGAAAGATTAGACAAACTAGATGAAACCCAAATAATACTAGT";
        //Weak bins: don't have RBS, should return true
        //String G10 = "TCTAGAGAAAGAGTAGATCAACTAGATGAAACCCAAATAATACTAGT";
        //String G12 = "TCTAGAGAAAGATTAGAGTCACTAGATGAAACCCAAATAATACTAGT";
        //Confirmation of Results
        System.out.println("Sequences that contain internal RBS:");

        //Should return False
        boolean result1 = checker.run(BBa_B0029);
        System.out.println("BBa_B0029 (false): " + result1);
        boolean result2 = checker.run(BBa_B0034);
        System.out.println("BBa_B0034 (false): " + result2);
        boolean result3 = checker.run(BBa_B0033);
        System.out.println("BBa_B0033 (false): " + result3);
        boolean result4 = checker.run(BBa_B0031);
        System.out.println("BBa_B0031 (false): " + result4);
        boolean result5 = checker.run(BBa_J61101);
        System.out.println("BBa_J61101 (false): " + result5);
        boolean result6 = checker.run(BBa_J61127);
        System.out.println("BBa_J61127 (false): " + result6);
        boolean result7 = checker.run(BBa_J61100);
        System.out.println("BBa_J61100 (false): " + result7);
        boolean result8 = checker.run(BBa_J61117);
        System.out.println("BBa_J61117 (false): " + result8);
        boolean result9 = checker.run(B2);
        System.out.println("B2 (false): " + result9);
        boolean result9p1 = checker.run(B3261);
        System.out.println("B3261 (false): " + result9p1);
        boolean result9p2 = checker.run(B1237);
        System.out.println("B1237 (false): " + result9p2);
        boolean result9p3 = checker.run(B0907);
        System.out.println("B0907 (false): " + result9p3);
        boolean result9p4 = checker.run(B2518);
        System.out.println("B2518 (false): " + result9p4);

        System.out.println("This algorithm doesn't cover E. coli cases:");

        //Should return True
        boolean result10 = checker.run(Ecoli1);
        System.out.println("Ecoli1 (true): " + result10);
        boolean result11 = checker.run(Ecoli2);
        System.out.println("Ecoli2 (true): " + result11);
        boolean result12 = checker.run(Ecoli3);
        System.out.println("Ecoli3 (true): " + result12);
        boolean result13 = checker.run(Ecoli4);
        System.out.println("Ecoli4 (true): " + result13);
        boolean result14 = checker.run(Ecoli5);
        System.out.println("Ecoli5 (true): " + result14);

        System.out.println("Random Sequences:");

        //Should return True
        boolean result15 = checker.run(Random);
        System.out.println("Random (true): " + result15);

        //Should return False
        boolean result16 = checker.run(Random1);
        System.out.println("Random1 (false): " + result16);

          //boolean result10 = checker.run(B3);
        //System.out.println("B3 (false): " + result10);
        //boolean result9 = checker.run(G10);
        //System.out.println("G10 (true): " + result9);
        //boolean result10 = checker.run(G12);
        //System.out.println("G12 (true): " + result10);
    }
}
