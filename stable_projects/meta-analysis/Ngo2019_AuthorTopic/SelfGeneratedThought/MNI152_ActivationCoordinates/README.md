# Data Release

`SelfGeneratedThought_AllCoordinates.txt` contains activation coordinates of 7 tasks of self-generated thought, which were provided by the following studies:

- Spreng, R.N., Mar, R.A. and Kim, A.S., 2009. The common neural basis of autobiographical memory, prospection, navigation, theory of mind, and the default mode: a quantitative meta-analysis. Journal of Cognitive Neuroscience.
- Mar, R.A., 2011. The neural bases of social cognition and story comprehension. Annual Review of Psychology.
- Sevinc, G. and Spreng, R.N., 2014. Contextual and perceptual brain processes underlying moral cognition: a quantitative meta-analysis of moral reasoning and moral emotions. PloS One.

Please cite the 3 papers above when using the activation coordinates.

----

# Input Text File Format
```
//Borg,2006: Nonmoral vs Moral
//Task: MORAL
-12  54 -19
-15  77 -17

//Rosenbaum, 2004: blocked routes vs. vowel counting baseline
//Task: NAV
21  20  49
-36   8  51
28 -37 -12
```

- Each experiment has the following format:
  - First line: start with "//", followed by the identifier of the experiment
  - Second line: start with "//Task: ", followed by one or more task names, delimited by commas.
  - Third line onward: coordinates of activation foci reported in MNI152 space.
- Information of one experiment is separated from another by a single empty line

----

# Pre-process Data
Activation coordinates from 167 studies of 7 tasks of self-generated thought was included in `SelfGeneratedThought_AllCoordinates.txt` following the format described above.

Please see `Ngo2019_AuthorTopic/utilities/preprocessing` for functions to pre-process and convert the raw input text file into Matlab data format used by the Collapsed Variational Bayes algorithm for the author-topic model.
