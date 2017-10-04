/**
 * 
 */
package lib.structure.capr;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

/**
 * @author Jan Hoinka
 * Standalone command line interface for CapR
 */
public class CapRCLI {

	private static CapR capr = new CapR();
	
	private static String dir = System.getProperty("user.dir");
	
	/**
	 * Main
	 * @param args list of RNA or DNA strings. Will produce 
	 * results named sequence1, sequence2, ... sequenceN in order
	 * of apperence in args
	 */
	public static void main(String[] args) {

		int sequence_counter = 1;
		
		for (String sequence : args) {
			
			// Make sure we have a valid sequence
			String seq = validateAlphabet(sequence);
			
			// Compute profile
			capr.ComputeStructuralProfile(seq.getBytes(), seq.length());
			
			// Get profile
			double[][] profile =  getMatrix(capr.getStructuralProfile(), sequence.length());
			
			// Output profile
			saveTxt(sequence_counter, String.format("Sequence%s_profile.txt", sequence_counter), sequence, profile);
			
			sequence_counter++;
		}

	}

	
	/**
	 * Beaks if alphabet is invalid, converts Us to Ts 
	 */
	private static String validateAlphabet(String seq) {
		
		//convert to upper
		seq = seq.toUpperCase();
		
		//validate
		StringBuilder sb = new StringBuilder();
		for ( byte c : seq.getBytes()) {
			
			//valid?
			if (c != 'A' || c != 'C' || c != 'G' || c != 'T' || c != 'U') {
				
				throw new RuntimeException(String.format("ERROR: Sequence %s contains invalid character %s", seq, c));
				
			}
			
			//convert (if required) and append
			sb.append( c == 'U' ? 'T' : c);
			
		}
		
		return sb.toString();
		
	}
	
	private static void saveTxt(int number, String filename, String sequence, double[][] probabilities) {
		
		Path path = Paths.get(dir , filename);
		
		//Use try-with-resource to get auto-closeable writer instance
		try (BufferedWriter writer = Files.newBufferedWriter(path))
		{
			// Store sequence ID
		    writer.write(String.format(">Sequence%s\n", number));

		    // Store nucleotides, tab separated
		    writer.write("\t");
		    writer.write(String.join("\t", sequence.split("")));
		    writer.write("\n");
		    
		    // Write probability for each context 
		    String[] hairpin  = new String[sequence.length()];
		    String[] inner    = new String[sequence.length()];
		    String[] bulge    = new String[sequence.length()];
		    String[] multi    = new String[sequence.length()];
		    String[] dangling = new String[sequence.length()];
		    String[] paired   = new String[sequence.length()];
		    
		    for (int pos=0; pos<probabilities.length; pos++) {
		    	
		    	hairpin[pos] = String.format("%.5f", probabilities[pos][0]);
			    inner[pos] = String.format("%.5f", probabilities[pos][1]);
			    bulge[pos] = String.format("%.5f", probabilities[pos][2]);
			    multi[pos] = String.format("%.5f", probabilities[pos][3]);
			    dangling[pos] = String.format("%.5f", probabilities[pos][4]);
			    paired[pos] = String.format("%.5f", probabilities[pos][5]);
		    	
		    }
		    
		    writer.write("\t");
		    writer.write(String.join("\t", hairpin));
		    writer.write("\n");
		    
		    writer.write("\t");
		    writer.write(String.join("\t", inner));
		    writer.write("\n");
		    
		    writer.write("\t");
		    writer.write(String.join("\t", bulge));
		    writer.write("\n");
		    
		    writer.write("\t");
		    writer.write(String.join("\t", multi));
		    writer.write("\n");
		    
		    writer.write("\t");
		    writer.write(String.join("\t", dangling));
		    writer.write("\n");
		    
		    writer.write("\t");
		    writer.write(String.join("\t", paired));
		    
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}

	/**
	 * Converts linear representation of matrix to rectangular
	 * @return [[h1,h2,...,hn],
	 *          [i1,i2,...,in],
	 *          [b1,b2,...,bn],
	 *          [m1,m2,...,mn],
	 *          [d1,d2,...,dn],
	 *          [p1,p2,...,pn]]
	 */
	private static double[][] getMatrix(double[] profile, int length){
		
		double[][] matrix = new double[length][6];
		
		for(int x = 0; x<length; x++) {
			
			matrix[x][0] = profile[x*0]; //H
			matrix[x][1] = profile[x*1]; //I
			matrix[x][2] = profile[x*2]; //B
			matrix[x][3] = profile[x*3]; //M
			matrix[x][4] = profile[x*4]; //D
			matrix[x][5] = 1.0 - (profile[x*0]+profile[x*1]+profile[x*2]+profile[x*3]+profile[x*4]); //P
			
		}
	
		return matrix;
		
	}
	
}
