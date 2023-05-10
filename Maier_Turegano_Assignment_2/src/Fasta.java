/**
 * Class Fasta
 * use to create a Fasta object
 */

// This is the provided Fasta class from Assignment 1

public class Fasta {
    private String header;
    private String sequence;

    public Fasta(String header, String sequence){
        this.header = header;
        this.sequence = sequence;
    }

    public String getSequence(){
        return this.sequence;
    }

    public String getHeader(){
        return this.header;
    }
}