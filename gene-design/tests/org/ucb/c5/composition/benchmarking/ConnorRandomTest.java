package org.ucb.c5.composition.benchmarking;

import org.ucb.c5.sequtils.HairpinCounter;
import org.ucb.c5.composition.checkers.*;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

import org.ucb.c5.composition.model.Composition;
import org.ucb.c5.composition.model.Host;
import org.ucb.c5.composition.model.Transcript;
import org.ucb.c5.composition.model.Construct;
import org.ucb.c5.sequtils.RevComp;

import java.util.ArrayList;
import javafx.util.Pair;
import java.util.List;
import java.util.Random;
import org.ucb.c5.composition.CompositionToDNA;

/**
 * @author Connor Tou Additional Classes Used not originally included:
 * generateRandomProtein, InternalRBSColiGenes, calculate GCContent
 *
 */
public class ConnorRandomTest {

    private static CompositionToDNA c2d;
    private static GenerateRandomProtein randomProtein;
    private static ForbiddenSequenceChecker sequenceChecker;
    private static HairpinCounter HairpinCounter;
    private static RevComp revcomp;
    private double Hairpins_first36;
    private static PromoterChecker promoterChecker;
    private static RepeatSequenceChecker repeatSequenceChecker;
    private static RNAInterferenceChecker rnaInterferenceChecker;

    @BeforeClass
    public static void setUpClass() throws Exception {
        c2d = new CompositionToDNA();
        c2d.initiate();

        randomProtein = new GenerateRandomProtein();
        randomProtein.initiate();

        sequenceChecker = new ForbiddenSequenceChecker();
        sequenceChecker.initiate();

        HairpinCounter = new HairpinCounter();
        HairpinCounter.initiate();

        revcomp = new RevComp();
        revcomp.initiate();

        promoterChecker = new PromoterChecker();
        promoterChecker.initiate();

        repeatSequenceChecker = new RepeatSequenceChecker();

        rnaInterferenceChecker = new RNAInterferenceChecker();
        rnaInterferenceChecker.initiate();

    }

    @Test
    public void test_random2() throws Exception {
        String promoter = "ttatgacaacttgacggctacatcattcactttttcttcacaa";
        String terminator = "TGCCTGGCGGCAGTAGCGCGGTGGTCCCACCTGACCCCATGCC";
        ArrayList<Double> Time = new ArrayList<>();
        ArrayList<Transcript> TranscriptList = new ArrayList<>();

        List<Pair> errorReasons = new ArrayList<>();
        int errors = 0;
        int GC_Content_Errors = 0;
        int RunTime_Errors = 0;
        int SecondaryStruct_Errors = 0;
        int SecondaryStruct36_Errors = 0;
        int ForbiddenSeq_Errors = 0;
        int Internal_RBS_Errors = 0;
        int M_Test_Errors = 0;
        int Proteins_Tested = 0;
        int Constitutive_promoter_Errors = 0;
        int Repeated_Sequence_Errors = 0;
        int RNA_Interference_Sequence_Errors = 0;
        Random rand = new Random();

        //Create 500 Random Proteins OR enter how many
        for (int i = 0; i < 100; i++) {
            ArrayList<String> proteins = new ArrayList<>();
            //Select length of protein sequence (default 6 - 100), and generate/add protein to list
            String protein = randomProtein.run(rand.nextInt(94) + 6);
            Proteins_Tested = Proteins_Tested + 1;

            //generateRandomProtein outputs protein with No starting "M" 20% of the time
            //Test Methionnine First Residue
            if (protein.charAt(0) != 'M') {
                try {
                    //Create a Composition, run Composition2DNA, and fetch run time
                    //Expect Excpetion
                    proteins.add(protein);
                    Composition comp = new Composition(Host.Ecoli, promoter, proteins, terminator);
                    long startTime = System.nanoTime();
                    Construct DNA = c2d.run(comp);
                    long endTime = System.nanoTime();
                    double totalTime = (endTime - startTime) * (.000000001); //nanoseconds to seconds conversion

                    //Failed to Throw Exception --> Continue and Report
                    Time.add(totalTime);
                    TranscriptList.add(DNA.getmRNAs().get(0));

                    //Filaed to Throw Exception --> Add Error Counts
                    errors = errors + 1;
                    M_Test_Errors = M_Test_Errors + 1;
                    continue;

                } catch (Exception e) {
                    continue;
                }
            }

            //Create a Composition, run Composition2DNA, and fetch run time
            proteins.add(protein);
            Composition comp = new Composition(Host.Ecoli, promoter, proteins, terminator);
            long startTime = System.nanoTime();
            Construct DNA = c2d.run(comp);
            long endTime = System.nanoTime();
            double totalTime = (endTime - startTime) * (.000000001); //nanoseconds to seconds conversion

            //Add runtime and Transcript for given protein calculation to respective Lists
            Time.add(totalTime);
            TranscriptList.add(DNA.getmRNAs().get(0));
        }

        List<String> dnaSequences = new ArrayList<>();

        //Build DNA Sequence String for each transcript
        TranscriptList.forEach((transcript) -> {
            StringBuilder dnaSeq = new StringBuilder();
            for (int i = 0; i < transcript.getCodons().length; i++) {
                dnaSeq.append(transcript.getCodons()[i]);
            }
            String addition = new String(dnaSeq);
            dnaSequences.add(addition);
        });

        //________________________________________________________________________________________________________
        //________________________________________________________________________________________________________
        //For each DNA sequence, run tests: (1) ForbiddenSeq/Restriction (2) GC Content (3) Internal RBS
        //(4)RunTime (5) Secondary Structure full sequence (6) Secondary Structure First 36 bp (7) Methionine
        //(8) Constitutive promoter (9) RepeatSequenceChecker (10) RNAInterferenceChecker
        for (int i = 0; i < dnaSequences.size(); i++) {
            String Sequence = dnaSequences.get(i);
            //ForbiddenSeq & Restriction sites
            //modified sequence checker to account for opp strand cut sites
            if (sequenceChecker.run(Sequence) == false) {
                errors = errors + 1;
                ForbiddenSeq_Errors = ForbiddenSeq_Errors + 1;
                Pair<String, String> pair = new Pair("Forbidden Seq or Restriction Site", dnaSequences.get(i));
                errorReasons.add(pair);
            }

            //GC Content
            //Top100 coligenes.txt --> <GCfreq> = 0.5145917192456388
            double GCcontent = 0.0;
            for (int j = 0; j < Sequence.length(); j++) {
                if (Sequence.charAt(j) == 'G' || Sequence.charAt(j) == 'C') {
                    GCcontent += 1;
                }
            }
            double GCfreq = GCcontent / (Sequence.length());
            if (GCfreq > 0.65 || GCfreq < 0.35) {
                errors = errors + 1;
                GC_Content_Errors = GC_Content_Errors + 1;
                Pair<String, String> pair = new Pair("GC Content", Sequence);
                errorReasons.add(pair);
            }

            //Internal RBS sites
            //Top 200 coligenes, 2-3 had 1 internalRBS, frequency = 0.0
            int InternalRBS_Count = 0;
            int index = -1;

            while (true) {
                index = Sequence.indexOf("AGGAGG", index + 1);
                if (index == -1) {
                    break;
                }
                InternalRBS_Count += 1;
            }

            if (InternalRBS_Count > 1) {
                errors = errors + 1;
                Internal_RBS_Errors = Internal_RBS_Errors + 1;
                Pair<String, String> pair = new Pair("Internal RBS", Sequence);
                errorReasons.add(pair);
            }

            //RunTime
            if (Time.get(i) > 15) {
                errors = errors + 1;
                RunTime_Errors = RunTime_Errors + 1;
                Pair<String, String> pair = new Pair("Time_Test", "for " + TranscriptList.size() + " proteins tested");
                errorReasons.add(pair);
            }

            //Secondary Structure for full sequence and first 36 bases
            //<hairpins in full sequence> = 25860.21 from top 100 coligenes.txt
            //<hairpins in first 36 bp> = 136.28
            if (dnaSequences.get(i).length() < 36) {
                Hairpins_first36 = HairpinCounter.run(Sequence.substring(0, Sequence.length() - 1));
            } else {
                Hairpins_first36 = HairpinCounter.run(Sequence.substring(0, 36));
            }
            double Hairpins_fullSeq = HairpinCounter.run(Sequence);
            if (Hairpins_fullSeq > 30000.0) {
                errors = errors + 1;
                SecondaryStruct_Errors = SecondaryStruct_Errors + 1;
                Pair<String, String> pair = new Pair("Secondary Structure full Sequence", Sequence);
                errorReasons.add(pair);
            }
            if (Hairpins_first36 > 200.0) {
                errors = errors + 1;
                SecondaryStruct36_Errors = SecondaryStruct36_Errors + 1;
                Pair<String, String> pair = new Pair("Secondary Structure first 36bp", Sequence);
                errorReasons.add(pair);
            }
            //Methionine First Residue Test done during random protein generation
            //Check if constitutive promoter pattern exists
            if(promoterChecker.run(Sequence) == false){
                errors = errors + 1;
                Constitutive_promoter_Errors = Constitutive_promoter_Errors + 1;
                Pair<String, String> pair = new Pair("Constitutive Promoter", Sequence);
                errorReasons.add(pair);
            }

            //RepearSequenceChecker; Check DNA sequence for the presence of 10 base pair repetitive sequence
            if(repeatSequenceChecker.run(Sequence) == false){
                errors = errors + 1;
                Repeated_Sequence_Errors = Repeated_Sequence_Errors + 1;
                Pair<String, String> pair = new Pair("Repeated Sequence", Sequence);
                errorReasons.add(pair);
            }

            //RNAInterferenceChecker: Check RNA interference from CDS in the composition
            if(rnaInterferenceChecker.run(Sequence) == false){
                errors = errors + 1;
                RNA_Interference_Sequence_Errors = RNA_Interference_Sequence_Errors + 1;
                Pair<String, String> pair = new Pair("Sequence with RNA Interference", Sequence);
                errorReasons.add(pair);
            }
        }

        //Print out Error Totals
        System.out.println();
        System.out.println("Proteins Tested: " + Proteins_Tested);
        System.out.println("Total Number of Errors:" + " " + errors);
        System.out.println("GC Content Errors: " + GC_Content_Errors);
        System.out.println("RunTime Errors: " + RunTime_Errors);
        System.out.println("Secondary Structure Full-Seq Errors: " + SecondaryStruct_Errors);
        System.out.println("SecondaryStruct First-36 Errors: " + SecondaryStruct36_Errors);
        System.out.println("Forbidden Sequence Errors: " + ForbiddenSeq_Errors);
        System.out.println("Internal RBS Errors: " + Internal_RBS_Errors);
        System.out.println("Methinione First Residue Errors: " + M_Test_Errors);
        System.out.println("Constitutive Promoter Errors: " + Constitutive_promoter_Errors);
        System.out.print("Sequence with Repeat Errors: " + Repeated_Sequence_Errors);
        System.out.print("Sequence with RNA Interference Error: " + RNA_Interference_Sequence_Errors);

        //List of Error Type & Sequence which error occured
        System.out.println(errorReasons);

        double threshold_5 = dnaSequences.size() * .05;
        double threshold_10 = dnaSequences.size() * .10;
        //No more than 10% can have any type of error
        assertTrue(errors < threshold_10);
        //No more than 5% can have errors of type:
        assertTrue(GC_Content_Errors < threshold_5);
        assertTrue(SecondaryStruct_Errors < threshold_5);
        assertTrue(SecondaryStruct36_Errors < threshold_5);
        //No Errors of type:
        assertTrue(RunTime_Errors == 0);
        assertTrue(ForbiddenSeq_Errors == 0);
        assertTrue(Internal_RBS_Errors == 0);
        assertTrue(M_Test_Errors == 0);
        assertTrue(Constitutive_promoter_Errors == 0);
        assertTrue(Repeated_Sequence_Errors == 0);
        assertTrue(RNA_Interference_Sequence_Errors ==0);


    }

    private static class GenerateRandomProtein {

        private static final char[] Choices = {'A', 'C', 'D', 'E', 'F', 'H', 'I', 'G', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};

        private Random rnd;
        private Random rnd2;

        public char getAA() {
            return Choices[rnd.nextInt(Choices.length)];
        }

        public void initiate() {
            rnd = new Random();
            rnd2 = new Random();
        }

        public String run(int size) {
            StringBuilder randProtein = new StringBuilder();
            int x = rnd2.nextInt(100);

            //20% of the time, generate a protein that Doesn't start with M
            if (x > 20) {
                randProtein.append("M");
            }
            for (int i = 0; i < size; i++) {
                randProtein.append(getAA());
            }
            return new String(randProtein);
        }
    }
}
