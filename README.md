# CapR4J
## CapR4J - A java port of CapR, for calculating probabilities that each RNA base position is located within secondary structural contexts.

CapR4J is a port of the CapR program orignially created by Fukunaga et al. It requires no external depencies is hence platform independent. In addition, an API is available which allows for the integration of CapR into bioinformatics pipelines in a programmatic manner. This port also generates a graphical representation of the structural context probabilities for each input sequence.

**Disclaimer: CapR is implemented and developed by Tsukasa Fukunaga et al. All intellectual credits of this work go to the original authors. My only contribution is the adaptation of the C source code to Java and the addition of the graphical representation.**

## Usage
```
java -jar capr4j INPUT_SEQUENCE
```

## Precomiled Jar File
Please see the [release section](https://github.com/drivenbyentropy/capr4j/releases) for the most recent version of CapR4J.

## Compile From Source
CapR4J is a maven project:
```
git clone https://github.com/drivenbyentropy/capr4j.git
cd capr4j
mvn package
```
The compiled jar file can then be found in `target\capr4j.jar`
