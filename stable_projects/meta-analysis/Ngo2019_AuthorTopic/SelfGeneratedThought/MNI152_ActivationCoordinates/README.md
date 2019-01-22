# Input file format
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

# Pre-process data
For example, self-generated thought data is saved in: `MNI152_activation_coordinates/SelfGeneratedThought_AllCoordinates.txt`

Please see `Ngo2019_AuthorTopic/utilities/preprocessing` for functions to pre-process and convert the raw input text file into Matlab data format used by the Collapsed Variational Bayese algorithm for the author-topic model.
