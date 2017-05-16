# How to use EDIBLE

EDIBLE ("Experimental Design Information By Likelihood Exploration") is a
computer program written in (hopefully strict) ANSI C to implement experimental
design ideas by Nick Goldman from his 1998 paper "Experimental Design in
Molecular Phylogenetics".

The original code is forked by Richel Bilderbeek in 2017. The original code
can be found on the `original` develop branch.

## To compile

### Contemporary

Using `qmake` and `make`:

```
qmake
make
```

See `.travis.yml` how the environment is set up.

### Classic

```
cd src
make
```

This option failed in the original version and is unmainted. Use the `original`
branch if you intend to use this. Feel encouraged to fix and submit a Pull Request.

## Command line options

  Most of the functionality of EDIBLE is controlled from the command line
(mainly because it seemed like a good idea at the time than any intention
of allowing it to be scripted.)

  The command line dump (from eg: running the program with no options) is:

Usage:

```
  edible 
    <-s [sample_file] sample_size sequence_length percentile> 
    <-m matrix_file> 
    <-p prob_file> 
    <-t> 
    <-b boot_strap_size> 
    <-d>
    <-h pi_a pi_c pi_g pi_t kappa> 
    <-g pi_a pi_c pi_g pi_t GTR_parameters>
    <-r ncat rates probs> 
    <-k> 
    <-c cache_size> 
    <-o> 
   tree_file 
   output_file
```

 * `edible` is the name of the executable to be run. Playing round with 
   different values of this parameter wastes time and has far too many
   interesting effects to be documented here.
 * `tree_file` is the name of the file which should contain a tree in
   bracket-colon format (see below) If not, bad things happen.
 * `output_file` is the name of the output file. 
   EDIBLE creates a file of the given name (overwriting if necessary) 
   and keeps going until it finishes or the hard disk dies.
 * `s` attempts to estimate the "percentile" percentile of the information
   distribution. The distribution in question is the expected amount of 
   information observed from a sequence with "sequence_length" nucleotides
   ('sites'). The estimate is based on "sample_size" samples, with the expected
   implications for variance (pun unintended, but noted). "sample_file" is an
   optional parameter, dumping the sample distribution (sample_size ordered 
   values of observed information per site.)
   * A program "pcalc" is provided to read in the file generated and 
     estimate arbitrary percentile points and sample mean.
 * `m` dumps the intermediate expectation matrices to "matrix_file". This 
   can be a lot of numbers on long runs, numerophiles only. Each element is
   separated by a comma, each row by a new-line.
 * `p` dumps the probabilities calculated for each possible combination
   of nucleotides. If information about only one branch (or node in a 
   clock-like tree) has been specified then the partial derivatives and
   information WRT that branch will also be dumped.
 * `t` dumps each tree used to a file (really only for debugging purposes;
   but may be useful to clarify what the program is doing; not that it isn't
   so obvious even a professor could understand it ;-)
 * `b` turns on sampling behaviour for large (and small) trees. The option 
   name is really a misnomer now, but is does what you'd expect - bases 
   estimates for the expected information on "boot_strap_size" samples. This 
   method shouldn't run into problems with large trees (unlike the exact
   calculation which requires integers of length 4^#leaves and will overflow
   on computers with finite registers for large trees.)
   Only silly people with time to waste on extra key strokes use the `b` and
  `s` options together.
 * `d` return the determinate of the information matrix specified by all
   the parameters chosen. If no parameters are specified, the determinate of 
   the entire matrix is returned. If only one parameter is specified, then
   the information of that parameter about itself is returned.
 * `h` tells EDIBLE to use HKY85 model of nucleotide substitution rather than
   JC. "pi_a", "pi_c", "pi_g", "pi_t" are the base frequencies of certain
   nucleotides. "kappa" is the transition/transversion ratio. Implementation
   of more bases is left to those who feel the need.
 * `g` use the GTR model of nucleotide substitution. The base frequencies are
   given as for HKY85. The five parameters for the substitution matrix should
   be in the same order as PAML: T<->C T<->A T<->G C<->A C<->G normalised so
   A<->G equals 1.0
 * `r` calculate information using rate variation. The number of categories,
   ncat, and the rate and proportion of sites for each category must also be
   specified. Rates must be positive and proportions must sum to 1.
   E.g. `-r 3 0.5 1.0 2.0 0.1 0.7 0.2`
   Three categories with rates 0.5, 1.0 and 2.0. The proportion of sites in
   each category is 0.1, 0.7 and 0.2 respectively.
   Note: this calculation is different from the calculating the information 
   separately for each rate category then forming a weighted sum, which would i
   be the information if we were told which category each site was in.
 * `k` turns _off_ calculations for information about kappa (this way round
   is sensible, regardless of what the boss thinks.)
 * `c` creates a cache of the "cache_size" most probable results 
   generated while sampling. Especially for small trees, this has the 
   potential to reduce sampling times by a considerable factor.

  The experimental version has an extra option <-i iteration_depth>
  "iteration_depth" is the distance from the maximum depth to precalculate.
See implementation.txt for more details. Defaults to zero if no value is
given.

### Examples

#### Start EDIBLE in interactive mode:

```
./edible test/primate.tree out.txt
```

#### Start EDIBLE in non-interactive mode:

```
./edible -o test/primate.tree out.txt
```

This will cause option `1. Calculate the expected information of the tree` being selected from the menu
and print the output file its content to the screen directly

## EDIBLE running

  Once successfully started, EDIBLE presents a small menu of options:

Options:  1. Calculate the expected information of the tree
          2. Scale tree by a given range of factors
          3. Scale a single branch by a given range of factors
          4. Slide a given branch along to others for given lengths
            Please choose:

  Option 1 performs the information calculation on tree provided, with the 
options specified and does very little else of interest. Interesting 
quotes, psychoanalysis and talking paper clips are the reserve of
"fortune", "eliza" and bad flash-backs respectively.

  Option 2 scales the entire tree by a given range of factors (you need to
provide min/max factors and the number of points to use) and does a number
of calculations uniformly distributed between the points (inclusively.)

  Option 3 does the same as option 2, except that only one branch length is
scaled, the rest of the tree remaining constant. For a computational point
of view, this make perfect sense with clock-like trees; people not wishing
to look foolish, beware.

  Option 4 allows one branch to be slid between two others. You need to
specify which branch is to be slid, which branch it is sliding too and 
which branch it is sliding from.
  Anyone commenting about over specification is probably correct.
  eg:
		| A -->
         E    B |    C
	------------------ 
           |
           | D

  Moving A in the direction indicated is sliding "A" from "B" to "C"

  Behaviour on anything other than simple ("binary" - three branches meet at
the node in question) trees has "undefined" effects; rest assured that I
know what happens.
  If the tree is clock-like, EDIBLE will correct the distances involved to
keep the evolutionary time to all other nodes constant.


## Tree specification (bracket,colon format)

  The program accepts trees written in the "bracket" format. 
  A tree is written as branches separated by commas, each branch may
either lead to a leaf or be another (sub)tree. The whole tree is encased in 
brackets.

  Format for a branch leading to a leaf is: "name:length". 'name' is the leaf
name, length is the length of the branch to the current node.
  The name must be less than 16 characters and must not be of the 
form "Node-n", where 'n' is an integer. Names are case sensitive and may be
composed of all standard alpha-numeric characters except white space,
colons and asterisks.

  For a subtree, the format is as if it was a normal tree followed by
":length", where 'length' is the branch length to the subtree.

  Branch lengths may be written as normal decimals, or in exponential format,
the exponent separated from the decimal by 'E' or 'e'.
        eg: 2.5674E-2 = 0.025674

eg: Consider the tree:

```
                      A       B        Length 1 to A   = 0.3
                       \     /                1 to C   = 0.15
                        1---2                 1 to 2   = 0.2
                       /     \                2 to B   = 0.1234
                      C       3-Dog           2 to 3   = 0.9862
                              |               3 to Dog = 2.3
                              E               3 to E   = 0.00231
```

*Note: the expression for the tree is not unique.*

  This would be written as (taking the common node of A and C to be the base
point)

```
(A:0.3,C:0.15,(B:1.234e-1,(Dog:0.23e+01,E:2.31E-3):0.9862):0.2)
```

  The program can convert the unrooted form of the tree to a root form, by 
indicating which branch contains the root node using an 'R'. This
branch _must_ be one the branches leading from the base node. The 'R' is
placed after the length of the branch it is on.
        eg: If root is on the branch from A to 1 in the above example:

```
(A:0.3R,C:0.15,(B:1.234e-1,(Dog:0.23e+01,E:2.31E-3):0.9852):0.2)
```

            If root is on the branch from 1 to 2 in the above example:

```
(A:0.3,C:0.15,(B:1.234e-1,(Dog:0.23e+01.E:2.31E-3):0.9852):0.2R)
```

            If root is on the branch from 2 to 3 in the above example:

```
((A:0.3,C:0.15):0.2,B:0.1234,(Dog:2.3,E:2.31e-3):0.9862R)
```

  If a result from an individual node is required, rather than the determinant
of the expected information matrix, the branch in question may be indicated
by using an asterisk instead of the colon preceding the length.
        eg: Want the expected information about branch from 1 to A:

```
(A*0.3,C:0.15,(B:1.234e-1,(Dog:0.23e+01,E:2.31E-3):0.9852):0.2)
```
            Want the expected information about branch from 3 to E:

```
(A:0.3,C:0.15,(B:1.234e-1,(Dog:0.23e+01,E*2.31E-3):0.9852):0.2)
```

  If the rooted tree is desired, the parameters change and the asterisk above
indicates that the information about the length from the node the branch
leads _from_ to the present day is desired. (From in the sense of away
from the imaginary root of the b-c format, not the specified tree root.)
  If information about the length leading to the root node is desired, this
may be indicated by an asterisk after the external most bracket.
  The behaviour of placing an asterisk here in the unrooted case is undefined.

  Alternatively, a node can be specified as the root instead. This is down by 
replacing the 'R' above by an 'r'. The root node will be the one that the
marked branch leads away from (WRT to the arbitrary 'root' node chosen -
"Node-0", the same as was chosen for putting the tree into brackets format.)

  It is allowable to specify that you require information about more than
one individual nodes (using asterisks) - the result is either the same as
running EDIBLE several times with each parameter individually specified
or the determinate of the matrix form by the information of each parameter
about the others, depending on whether the "d" flag is set or not.
  On the other hand, specifying more than one root is asking to play
result roulette (the table is probably favours those with the patience to
read the source. It isn't worth it.)

## Portability considerations and other features

EDIBLE is bug free, removing unwanted features is an exercise for the user.
There are no integer overflows for large(?) trees in this implementation
(written for a hypothetical Minsky machine). Your mileage on finite
architectures may vary.
  For information, the tree size is limited by the integer size on the
machine. Most current machines have 31bit (signed) integers (lets see how
quickly this looks dated.) We need to differentiate every possible 
nucleotide combination for exact calculation and so are limited to 15 leaf
trees on most machines. This rises to 31 leaves on platforms such as the 
DEC alpha and most of the next generation of micro processors.
  Sampling/"boot strapping" for values has no such limitation.
   
  On most machines memory is not a concern, using less than 1Mb on a 13
leaf tree. An important exception to this is the experimental version,
which has memory usages highly dependent on the iteration depth and 
tree topology. Please see implementation.txt for more details.

## Definitions

  In any case of doubt, my definition of a large tree is:
A maximal acyclic graph of the minimal order to cause the effect in question.
