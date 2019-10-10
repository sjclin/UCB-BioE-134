package org.ucb.c5.composition;

import javafx.util.Pair;
import org.ucb.c5.composition.checkers.ForbiddenSequenceChecker;
import org.ucb.c5.composition.checkers.RepeatSequenceChecker;
import org.ucb.c5.composition.model.RBSOption;
import org.ucb.c5.composition.model.Transcript;
import org.ucb.c5.sequtils.HairpinCounter;
import org.ucb.c5.sequtils.Translate;

import java.util.*;

/**
 * This reverse translates a protein sequence to a DNA and chooses an RBS. It
 * uses the highest CAI codon for each amino acid in the specified Host.
 *
 * @author J. Christopher Anderson
 * @author Stephen Lin
 */
public class TranscriptDesigner {

    private RBSChooser2 rbsChooser;
    private Translate translator;
    private ForbiddenSequenceChecker forbiddenSeqChecker;
    private Map<Character, String[]> aaMap;
    private String aaCharacters;
    private Comparator<String> gcComparator;
    private HairpinCounter hpc;
    private RepeatSequenceChecker repSeqChecker;

    public void initiate() throws Exception {
        //Initialize the RBSChooser
        rbsChooser = new RBSChooser2();  //Use new algorithm to choose RBS
        rbsChooser.initiate();

        translator = new Translate();
        translator.initiate();

        forbiddenSeqChecker = new ForbiddenSequenceChecker();
        forbiddenSeqChecker.initiate();

        aaMap = new HashMap<>();
        Character[] aminoAcids = new Character[]{'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
                                           'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'};
        String[][] codons = new String[][]{{"GCA", "GCT", "GCC", "GCG"}, {"AGA", "AGG", "CGA", "CGT", "CGC", "CGG"},
            {"AAT", "AAC"}, {"GAT", "GAC"}, {"TGT", "TGC"}, {"CAA", "CAG"}, {"GAA", "GAG"}, {"GGA", "GGT", "GGC", "GGG"},
            {"CAT", "CAC"}, {"ATA", "ATT", "ATC"}, {"TTA", "TTG", "CTA", "CTT", "CTC", "CTG"}, {"AAA", "AAG"}, {"ATG"},
            {"TTT", "TTC"}, {"CCA", "CCT", "CCC", "CCG"}, {"AGT", "AGC", "TCA", "TCT", "TCC", "TCG"},
            {"ACA", "ACT", "ACC", "ACG"}, {"TGG"}, {"TAT", "TAC"}, {"GTA", "GTT", "GTC", "GTG"}};
        for (int i = 0; i < 20; i++) {
            aaMap.put(aminoAcids[i], codons[i]);
        }

        aaCharacters = "([ARNDCQEGHILKMFPSTWYV])+";

        gcComparator = Comparator.comparingDouble(this::gcContent);

        hpc = new HairpinCounter();
        hpc.initiate();

        repSeqChecker = new RepeatSequenceChecker();
    }

    public Transcript run(String peptide, Set<RBSOption> ignores) throws Exception {
        if (!peptide.matches(aaCharacters)) {
            throw new IllegalArgumentException("peptide contains non-amino acid character");
        }
        if (peptide.charAt(0) != 'M') {
            throw new IllegalArgumentException("peptide must start with methionine");
        }
        String[] codons = new String[peptide.length()];
        StringBuilder cds = new StringBuilder();
        //Use sliding window approach to choose codons
        for(int i = 0; i < peptide.length(); i++) {
            String aaWindow = peptide.substring(i, Math.min(i + 3, peptide.length()));
            int start = Math.max(0, i * 3 - 9);
            int end = Math.min(i * 3, cds.length());
            String preamble = cds.substring(start, end);
            List<String> synCodonsList = synonymousCodons(aaWindow);
            //Sort synCodons such that codons balancing GC content are prioritized
            if (gcContent(cds.toString()) < 0.5) {
                synCodonsList.sort(gcComparator.reversed());
            } else {
                synCodonsList.sort(gcComparator);
            }
            List<String> validCodons = new ArrayList<>();
            List<Pair<String, Double>> hpScores = new ArrayList<>();
            for (String sCodons: synCodonsList) {
                String testSeq = preamble + sCodons;
                if (forbiddenSeqChecker.run(testSeq)) { // && repSeqChecker.run(cds.toString())) for cds.length() >= 10 {
                    validCodons.add(sCodons);
                    hpScores.add(new Pair<>(sCodons, hpc.run(testSeq))); //TODO: change to hpc.run(cds.toString() + sCodons) or similar for correctness
                }
            }
            if (validCodons.isEmpty()) {
                throw new Exception("no valid codons found for peptide " + peptide);
            }
            //Sort hpCodons such that codons minimizing secondary structure are prioritized
            hpScores.sort(Comparator.comparing(Pair::getValue));
            List<String> hpCodons = new ArrayList<>();
            for (Pair<String, Double> p: hpScores) {
                hpCodons.add(p.getKey());
            }
            int bestRank = Integer.MAX_VALUE;
            String bestCodons = "";
            //Choose the codon that balances GC content and minimizes secondary structure
            for (int j = 0; j < validCodons.size(); j++) {
                String vCodons = validCodons.get(j);
                int rank = j + hpCodons.indexOf(vCodons);
                if (rank < bestRank) {
                    bestRank = rank;
                    bestCodons = vCodons;
                }
            }
            cds.append(bestCodons, 0, 3);
            codons[i] = bestCodons.substring(0, 3);
        }
        //Choose an rbs
        RBSOption selectedRBS = rbsChooser.run(cds.toString(), ignores);
        //Construct the Transcript and return it
        Transcript out = new Transcript(selectedRBS, peptide, codons);
        return out;
    }

    private double gcContent(String seq) {
        double gc = 0;
        for (int i = 0; i < seq.length(); i++) {
            char base = seq.charAt(i);
            if (base == 'G' || base == 'C') {
                gc++;
            }
        }
        return gc / seq.length();
    }

    private List<String> synonymousCodons(String aaWindow) {
        return synonymousCodons(aaWindow, new ArrayList<>());
    }

    /**
     * Given a String aaWindow of amino acids,
     * return all combinations of synonymous
     * codons that can possibly encode it in a List.
     * @param aaWindow Peptide sequence
     * @return synCodons
     */
    private List<String> synonymousCodons(String aaWindow, List<String> synCodons) {
        String[] currCodons = aaMap.get(aaWindow.charAt(0));
        if (aaWindow.length() == 1) {
            Collections.addAll(synCodons, currCodons);
            return synCodons;
        }
        synonymousCodons(aaWindow.substring(1), synCodons);
        String[] synCodonArr = synCodons.toArray(new String[synCodons.size()]);
        synCodons.clear();
        for (String synCodon: synCodonArr) {
            for (String validCodon: currCodons) {
                synCodons.add(validCodon + synCodon);
            }
        }
        return synCodons;
    }
}
