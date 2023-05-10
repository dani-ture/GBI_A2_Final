/**
 * Assignment 02
 * Authors: ANNIKA MAIER & DANIEL TUREGANO
 */
import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

// This is our FastaReader class from Assignment 1
public class FastaReader {

    // Methods


    /**
     *  2.1: Method to read in a fastA file. It generates new Fasta objects for each entry.
     * @param filepath
     * @return List of Fasta objects
     */
    public List<Fasta> readInFasta(String filepath){
        List<Fasta> fastaEntries = new ArrayList<>();

        try {

            //Read one file by specifying its path
            BufferedReader br = new BufferedReader(new FileReader(filepath));

            //First line of a FASTA file should always be a header
            String line = br.readLine();
            String header;

            // Keep reading while the file is non-empty
            while (line != null) {
                header = line;
                StringBuilder sequence = new StringBuilder();
                line = br.readLine();

                // Keep reading while the file is non-empty and the line does not begin with '>'
                // The lines within each sequence will be appended
                while ((line != null) && !(line.charAt(0) == '>')){
                    sequence.append(line);
                    line = br.readLine();
                }
                //Create and add Fasta object
                fastaEntries.add(new Fasta(header,sequence.toString()));
            }

            br.close();
        }
        catch (IOException e) {
            //Print the error code, in case the Input or Output operations failed
            e.printStackTrace();
        }

        return fastaEntries;
    }


    /**
     * 2.2: Returns the sequence length
     * @param fasta (Fasta object)
     * @return (length of the sequence)
     */
    public int calculateSequenceLength(Fasta fasta){
        int length = 0;
        length = fasta.getSequence().length();

        return length;
    }


    /**
     *  2.3 Calculates frequency for each base in percentage
     * @param fasta (Fasta object)
     * @return (a map containing all the bases and their frequencies)
     */

    public Map<Character, Double> calculateBaseFrequency(Fasta fasta) {
        Map<Character, Double> baseFrequencies = new HashMap<>();

        int ACount=0;
        int TCount=0;
        int CCount=0;
        int GCount=0;
        String seq = fasta.getSequence();
        int n = seq.length();

        for (String s : seq.split("")) {
            switch (s) {
                case "A":
                    ACount++; break;
                case "T":
                    TCount++; break;
                case "C":
                    CCount++; break;
                case "G":
                    GCount++; break;
            }
        }

        baseFrequencies.put('A', (double) ACount/n*100);
        baseFrequencies.put('T', (double) TCount/n*100);
        baseFrequencies.put('C', (double) CCount/n*100);
        baseFrequencies.put('G', (double) GCount/n*100);

        return baseFrequencies;

    }

    /**
     * 2.4  Writes a Fasta object to a given file path. Must adhere to fasta format conventions
     * @param fasta (Fasta object)
     * @param filepath (Path to write the file to)
     */
    public void writeOutFasta(Fasta fasta, String filepath){

        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(filepath));

            // Write header, then sequence, repeat

            writer.write(fasta.getHeader()+ "\n");
            writer.write(fasta.getSequence());

            writer.close();
            System.out.println("File has been saved");

        }
        catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * 2.5  Converts the fastA Sequences to the RY-sequences
     * @param fasta (Fasta object)
     * @return (RY-sequence of the provided sequence)
     */

    public String fastaToRYSequence(Fasta fasta){

        return fasta.getSequence().replace("A","R")
                .replace("G","R")
                .replace("C","Y")
                .replace("T", "Y");
    }

    /**
     * 2.6  Saves the RY-Sequence to a file called ry_inputfilename.fasta
     * @param fasta (Fasta object)
     * @param save_location (the desired location where the file will be saved)
     * @param filename (filename of the output file. It will be transformed to ry-filename.fasta).
     */

    public void saveRYSequence(Fasta fasta,String save_location, String filename){
        String rySequence = fastaToRYSequence(fasta);

        try {

            BufferedWriter writer = new BufferedWriter(new FileWriter(save_location + "ry_"+ filename + ".fasta"));

            writer.write(fasta.getHeader() + " ry_version\n");
            writer.write(rySequence);
            writer.close();

            System.out.println("File has been saved");

        }

        catch (IOException e) {
            e.printStackTrace();
        }
    }
}