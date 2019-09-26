package org.ucb.c5.composition.checkers;

import java.util.HashSet;
import java.util.Set;

/**
 * Check DNA sequence for the presence of 10 base pair repetitive sequence
 *
 * @author Kevin Hong
 */
public class RepeatSequenceChecker {

    //Decide the length of repeating DNA sequence
    private final int length = 10;

    public boolean run(String seq) throws Exception {

        //Uppercase all base pairs
        seq = seq.toUpperCase();

        if (!seq.matches("[ATCG]+")) {
            throw new Exception();
        }

        //Ensure sequence is long enough
        if (seq.length() <= length) {
            throw new Exception();
        }

        //Create Hashset containing unique sequences with decided length
        Set<String> existingSequence = new HashSet<>();

        //Scan through DNA sequence
        for (int i = 0; i < seq.length() - length + 1; i++) {
            String existingSeq = seq.substring(i, i + length);

            // Check if the current sequence is a repeated sequence of any unique sequences in Hashset
            if (existingSequence.contains(existingSeq)) {
                return false;
            }

            //Add this unique sequence to Hashset
            existingSequence.add(existingSeq);
        }
        return true;
    }

    public static void main(String[] args) throws Exception {

        //Create example DNA sequence
        String seqNoRepeat = "ATGCAGTAAAATTTCCG";
        String repeatSeq = "ATGCAGTAAAATTTCCGATGCAGTAAAATTTCCG";

        //Run and output result
        RepeatSequenceChecker RSC = new RepeatSequenceChecker();

        boolean result = RSC.run(seqNoRepeat);
        System.out.println("seq with No Repeat (true):" + result);

        result = RSC.run(repeatSeq);
        System.out.println("Repeatitive sequence (false):" + result);
    }
}
