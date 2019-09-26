package org.ucb.c5.composition.checkers;

import org.ucb.c5.sequtils.RevComp;

/**
 * Checks a sequence for internal constitutive promoters
 *
 * Uses a position weight matrix to identify E. coli sigma70 binding sites. Does
 * not look for other sigma factors.
 *
 * https://en.wikipedia.org/wiki/Position_weight_matrix The wikipedia entry does
 * a good job of explaining PWM and PFM
 *
 * @author Joanne Chang student joanne91218
 */
public class PromoterChecker {

    private RevComp revcomp;
    private double[][] pwm;

    public void initiate() {
        revcomp = new RevComp();
        revcomp.initiate();

        //Hard code a PFM matrix assuming 12 input seq with 30 bp with constitutive promoters motifs and rest of the parts perfectly random
        int[][] pfm = new int[][]{
            {0, 0, 0, 12, 0, 12, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 12, 0, 12, 12, 0}, // A
            {0, 0, 0, 0, 12, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0}, // C
            {0, 0, 12, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0}, // G
            {12, 12, 0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 12, 0, 12, 0, 0, 12} // T
        };

        //Create the PWM
        int ncols = pfm[0].length;
        pwm = new double[4][ncols];

        //Calculate the values for the PWM
        for (int x = 0; x < ncols; ++x) {
            double total = 0;
            for (int y = 0; y < 4; ++y) {
                total += pfm[y][x];
            }
            for (int y = 0; y < 4; ++y) {
                double prob_base = 0.25;
                double freq = pfm[y][x];
                double w = (Math.log((freq + Math.sqrt(total) * prob_base) / (total + Math.sqrt(total)) / prob_base)) / (Math.log(2));
                pwm[y][x] = w;
            }
        }
    }

    /**
     * Checks a DNA sequence for constitutive sigma70 promoters
     *
     *
     * @return true if none expected; false if contains a constitutive promoter
     */
    public boolean run(String bob) {
        bob = bob.toUpperCase();
        String rc = revcomp.run(bob);
        String combined = bob + "x" + rc;

        //computing the score by the sliding frame to target the constitutive promoter
        int slidingframe = 29;
        double score = 0.00;
        int pwm_ncol = combined.length();

        //JCA:  Adjusted threshold empirically based on J23119 promoter sequences
        //        double threshold = 20.815721036276983;  // threshold value is a common value defined after passing 12 constitutive promoter containing sequences
        double threshold = 9.134;

        for (int i = 0; i + slidingframe <= pwm_ncol; i++) {
            score = 0.00;
            String partseq = combined.substring(0 + i, slidingframe + i);
            for (int x = 0; x < partseq.length(); x++) {
                int y = -1;
                switch (partseq.charAt(x)) {
                    case 'A':
                        y = 0;
                        break;
                    case 'C':
                        y = 1;
                        break;
                    case 'G':
                        y = 2;
                        break;
                    case 'T':
                        y = 3;
                        break;
                }
                if (y == -1) {
                    continue;
                }
                score += pwm[y][x];
            }
            if (score >= threshold) {           //if the score is greater or equal to the threshold value, pattern recognized! return false
                return false;
            } else {
                continue;
            }
        }

        return true;
    }

    public static void main(String[] args) {
        PromoterChecker checker = new PromoterChecker();
        checker.initiate();

        String ex1 = "CAAGGGTATAATTAGATTTGACAAGAAGGAGTGAATCAATTATAAT";
        String constitutive = "TTGACAATTAATCATCGAACTAGTATAAT";
        String constitutiveBroken = "TTctgAATTAATCATCGAACTAGgcgAAT";

        //High Pcon's, from highest (23119) to still pretty high (J23101)
        String j23119 = "ttgacagctagctcagtcctaggtataatgctagc";
        String J23100 = "ttgacggctagctcagtcctaggtacagtgctagc";
        String J23102 = "ttgacagctagctcagtcctaggtactgtgctagc";
        String J23104 = "ttgacagctagctcagtcctaggtattgtgctagc";
        String j23101 = "tttacagctagctcagtcctaggtattatgctagc";

        //Low Pcon's, from off (j23112) to very low (J23117)
        String J23112 = "ctgatagctagctcagtcctagggattatgctagc";
        String J23103 = "ctgatagctagctcagtcctagggattatgctagc";
        String J23113 = "ctgatggctagctcagtcctagggattatgctagc";
        String J23109 = "tttacagctagctcagtcctagggactgtgctagc";
        String J23117 = "ttgacagctagctcagtcctagggattgtgctagc";

        //returns false if matches constitutive promoter patches
        boolean result = checker.run(constitutive);
        System.out.println("constitutive (false): " + result);

        //returns true if matches the constitutive promoter is broken
        boolean result2 = checker.run(constitutiveBroken);
        System.out.println("constitutiveBroken (true): " + result2);

        boolean result3 = checker.run(ex1);
        System.out.println("ex1: " + result3);

        //String Pcon's (false)
        boolean result4 = checker.run(j23119);
        System.out.println("j23119 (false): " + result4);

        boolean result5 = checker.run(J23100);
        System.out.println("J23100 (false): " + result5);

        boolean result6 = checker.run(J23102);
        System.out.println("J23102 (false): " + result6);

        boolean result7 = checker.run(J23104);
        System.out.println("J23104 (false): " + result7);

        boolean result8 = checker.run(j23101);
        System.out.println("j23101 (false): " + result8);

        //Low PCon's (true)
        boolean result9 = checker.run(J23112);
        System.out.println("J23112 (true): " + result9);

        boolean result10 = checker.run(J23103);
        System.out.println("J23103 (true): " + result10);

        boolean result11 = checker.run(J23113);
        System.out.println("J23113 (true): " + result11);

        boolean result12 = checker.run(J23109);
        System.out.println("J23109 (true): " + result12);

        boolean result13 = checker.run(J23117);
        System.out.println("J23117 (true): " + result13);

    }
}
