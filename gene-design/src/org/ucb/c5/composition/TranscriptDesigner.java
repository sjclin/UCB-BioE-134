package org.ucb.c5.composition;

import javafx.util.Pair;
import org.ucb.c5.composition.model.RBSOption;
import org.ucb.c5.composition.model.Transcript;
import org.ucb.c5.sequtils.Translate;
import org.ucb.c5.utils.FileUtils;

import java.util.*;

/**
 * This reverse translates a protein sequence to a DNA and chooses an RBS. It
 * uses the highest CAI codon for each amino acid in the specified Host.
 *
 * @author J. Christopher Anderson
 */
public class TranscriptDesigner {

    private Map<Character, String> aminoAcidToCodon;
    private RBSChooser2 rbsChooser;
    private Translate translator;
    private Map<String, List<Pair<String, Double>>> aaMap;
    private Map<String, List<Double>> aaToCumulativeFrequencies;

    public void initiate() throws Exception {
        //Initialize the RBSChooser
        rbsChooser = new RBSChooser2();  //Use new algorithm to choose RBS
        rbsChooser.initiate();
        
        //Construct a map between each amino acid and the highest-CAI codon for E coli
        aminoAcidToCodon = new HashMap<>();

        aminoAcidToCodon.put('A', "GCG");
        aminoAcidToCodon.put('C', "TGC");
        aminoAcidToCodon.put('D', "GAT");
        aminoAcidToCodon.put('E', "GAA");
        aminoAcidToCodon.put('F', "TTC");
        aminoAcidToCodon.put('G', "GGT");
        aminoAcidToCodon.put('H', "CAC");
        aminoAcidToCodon.put('I', "ATC");
        aminoAcidToCodon.put('K', "AAA");
        aminoAcidToCodon.put('L', "CTG");
        aminoAcidToCodon.put('M', "ATG");
        aminoAcidToCodon.put('N', "AAC");
        aminoAcidToCodon.put('P', "CCG");
        aminoAcidToCodon.put('Q', "CAG");
        aminoAcidToCodon.put('R', "CGT");
        aminoAcidToCodon.put('S', "TCT");
        aminoAcidToCodon.put('T', "ACC");
        aminoAcidToCodon.put('V', "GTT");
        aminoAcidToCodon.put('W', "TGG");
        aminoAcidToCodon.put('Y', "TAC");

        String codonFrequencies = FileUtils.readFile("src/org/ucb/c5/composition/data/codon_usage.txt");
        String[] lines = codonFrequencies.split("\\r|\\r?\\n");
        translator = new Translate();
        translator.initiate();
        aaMap = new HashMap<>();
        for (String line: lines) {
            String[] values = line.split("\t");
            String codon = values[0];
            String aa = values[1];
            double frequency = Double.parseDouble(values[2]);
            Pair<String, Double> entry = new Pair<>(codon, frequency);
            if (!aaMap.containsKey(aa)) {
                aaMap.put(aa, new ArrayList<>());
            }
            aaMap.get(aa).add(entry);
        }

        aaToCumulativeFrequencies = new HashMap<>();
        for (String aa: aaMap.keySet()) {
            List<Double> cumulativeSums = new ArrayList<>();
            double currSum = 0;
            for (Pair<String, Double> entry: aaMap.get(aa)) {
                double entryValue = entry.getValue();
                cumulativeSums.add(currSum + entryValue);
                currSum += entryValue;
            }
            cumulativeSums.set(cumulativeSums.size() - 1, 1.0);
            aaToCumulativeFrequencies.put(aa, cumulativeSums);
        }
    }

    private String getCodon(String aa, Random rng) throws IllegalArgumentException {
        if (!aaMap.containsKey(aa)) {
            throw new IllegalArgumentException("Illegal character " + aa + " in peptide string");
        }
        List<Double> cumulativeFrequencies = aaToCumulativeFrequencies.get(aa);
        double randNum = rng.nextDouble();
        int index = 0;
        for (double f: cumulativeFrequencies) {
            if (randNum <= f) {
                break;
            }
            index++;
        }
        String codon = aaMap.get(aa).get(index).getKey();
        return codon;
    }

    public Transcript run(String peptide, Set<RBSOption> ignores) throws Exception {        
        //Choose codons for each amino acid
        String[] codons = new String[peptide.length()];
        Random rng = new Random(23498752);
        for(int i=0; i < peptide.length(); i++) {
            String aa = Character.toString(peptide.charAt(i));
            //Use Monte Carlo method to pick codon at random
            String codon = getCodon(aa, rng);
            codons[i] = codon;
        }
        
        //Choose an RBS
        StringBuilder cds = new StringBuilder();
        for(String codon : codons) {
            cds.append(codon);
        }
        RBSOption selectedRBS = rbsChooser.run(cds.toString(), ignores);
        
        //Construct the Transcript and return it
        Transcript out = new Transcript(selectedRBS, peptide, codons);
        return out;
    }
}
