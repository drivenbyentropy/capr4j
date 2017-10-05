# CapR4J
## CapR4J - A java port of CapR, for calculating probabilities that each RNA base position is located within secondary structural contexts.

CapR4J is a port of the [CapR algorithm](https://github.com/fukunagatsu/CapR) orignially created by Fukunaga et al and implemented in C++. It requires no external depencies is hence platform independent. In addition, an API is available which allows for the integration of CapR into bioinformatics pipelines in a programmatic manner. This port also generates a graphical representation of the structural context probabilities for each input sequence.

**Disclaimer: CapR is implemented and developed by Tsukasa Fukunaga et al. All intellectual credits of this work go to the original authors. My only contribution is the adaptation of the C source code to Java and the addition of the graphical representation.**

## Usage
```
java -jar capr4j INPUT_SEQUENCE
eg. java -jar GGGAGACAAGAATAAACGCTCAACAAACAACAGTACGTAGCATGCATGCTAGCTAGCTACTATGGGATTCGACAGGAGGCTCACAACAGGC
```

## Precomiled Jar File
Please see the [release section](https://github.com/drivenbyentropy/capr4j/releases) for the most recent version of CapR4J.

## Compile From Source
CapR4J is a maven project:
```
git clone https://github.com/drivenbyentropy/capr4j.git
cd capr4j
mvn install
```
The compiled jar file can then be found in `target\capr4j.jar`

## CapR4J API
If you want to integrate `CapR4J` in your own pipeline, you can use the provided API as such
```java

public static void main(String[] args) {
	private static CapR capr = new CapR();
    
	// The sequence to be predicted
	String seq = new String("GGGAGACAAGAATAAACGCTCAACAAACAACAGTACGTAGCATGCATGCTAGCTAGCTACTATGGGATTCGACAGGAGGCTCACAACAGGC");
			
	// Compute profile
	capr.ComputeStructuralProfile(seq.getBytes(), seq.length());

	// Get profile. Returns a double array containing a linear representation of the context probabilities
	// in the form of [h1,h2,...,hn,i1,i2,...,in,b1,b2,...,bn,m1,m2,...,mn,d1,d2,...,dn] where h,i,b,m,d
	// stand for hairpin, inner loop, bulge loop, multiple loop, and dangling end respectively. The probability
	// of a positition being paired is can be computed as 1-(h+i+b+m+d).
	double[] profile =  capr.getStructuralProfile();
}
```
