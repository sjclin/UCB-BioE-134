package org.bioe134.crispr;

import javafx.util.Pair;

import java.sql.SQLOutput;

/**
 *
 * @author J. Christopher Anderson
 */
public class Design_CRISPR_Oligos {
    
    public void initiate() throws Exception {

    }
    
    public Pair<String,String> run(String cds) throws Exception {
        if (!cds.matches("([ATCG])+"))
            throw new Exception();
        String PAMSequence = "GG";
        String searchSequence = cds.substring(21, cds.length() - 1);
        int pamIndex = searchSequence.indexOf(PAMSequence);
        if (pamIndex < 0)
            throw new Exception();
        String oligo1 = "CATAACTAGT" + cds.substring(pamIndex, pamIndex + 20) + "GTTTTAGAGCTAGAAATAGCAAG";
        String oligo2 = "TCAGACTAGTATTATACCTAGGACTGAGCTAG";
        return new Pair<>(oligo1, oligo2);
    }
    
    public static void main(String[] args) throws Exception {
        //Create some example arguments, here the amilGFP coding sequence
        String cds = "ATGTCTTATTCAAAGCATGGCATCGTACAAGAAATGAAGACGAAATACCATATGGAAGGCAGTGTCAATGGCCATGAATTTACGATCGAAGGTGTAGGAACTGGGTACCCTTACGAAGGGAAACAGATGTCCGAATTAGTGATCATCAAGCCTGCGGGAAAACCCCTTCCATTCTCCTTTGACATACTGTCATCAGTCTTTCAATATGGAAACCGTTGCTTCACAAAGTACCCGGCAGACATGCCTGACTATTTCAAGCAAGCATTCCCAGATGGAATGTCATATGAAAGGTCATTTCTATTTGAGGATGGAGCAGTTGCTACAGCCAGCTGGAACATTCGACTCGAAGGAAATTGCTTCATCCACAAATCCATCTTTCATGGCGTAAACTTTCCCGCTGATGGACCCGTAATGAAAAAGAAGACCATTGACTGGGATAAGTCCTTCGAAAAAATGACTGTGTCTAAAGAGGTGCTAAGAGGTGACGTGACTATGTTTCTTATGCTCGAAGGAGGTGGTTCTCACAGATGCCAATTTCACTCCACTTACAAAACAGAGAAGCCGGTCACACTGCCCCCGAATCATGTCGTAGAACATCAAATTGTGAGGACCGACCTTGGCCAAAGTGCAAAAGGCTTTACAGTCAAGCTGGAAGCACATGCCGCGGCTCATGTTAACCCTTTGAAGGTTAAATAA";
        
        //Instantiate and initiate the Function
        Design_CRISPR_Oligos func = new Design_CRISPR_Oligos();
        func.initiate();
        
        //Run the function on the example
        Pair<String,String> oligos = func.run(cds);
        
        //Print out the result
        System.out.println("oligo1: " + oligos.getKey());
        System.out.println("oligo2: " + oligos.getValue());
    }
}