#  PLUMED Masterclass 21.3: Umbrella Sampling

This lesson was given as part of the PLUMED masterclass series in 2021.  It includes:

* A video that explain the theory covered and a second video which shows you how the exercises should be completed.
* A series of exercises that you should try to complete yourself.
* Some supplementary python notebooks that provide further background information on the exercise.

The flow chart shown below indicates the order in which you should consult the resources.  You can click on the nodes to access the various resources.  Follow the thick black lines for the best results.  The resources that are connected by dashed lines are supplmentary resources that you may find useful when completing the exercise. 

This lesson was the third masterclass in the 2021 series.  You will likely be able to complete the exercise without completing all the exercises in the first two masterclasses in the series.  However, the first masterclass contains instructions for installing PLUMED using conda that you may need to consult if you have not installed it already.  Meanwhile, the second masterclass contains more of the theory behind the block averaging methods that you will be using to complete the exercises here.

```mermaid
flowchart TB;
  A[PLUMED intro] -.-> C[Lecture I] 
  B[Block averaging theory] -.-> C
  C ==> D[Instructions];
  D ==> E[Lecture II];
  D --> F[solution];
  click A "ref1" "This elementary lesson shows you how to install PLUMED and some basic and more advanced syntax. You will only need to complete up to exercise 1 of this earlier lesson to complete these exercises on umbrella sampling.";
  click B "ref2" "This lesson introduces you to more of the theory of block averaging.  You use this method here to calculate the error bars on the averages that you obtain from your umbrella sampling simulations";
  click C "video1" "A lecture that was given on February 15th 2021 as part of the plumed masterclass series that introduces you to the exercises in this lesson"
  click D "INSTRUCTIONS.md" "The instructions for the exercises";
  click E "video2" "A lecture that was given on February 22nd 2021 as part of the plumed masterclass series that goes through the solutions to the exercises in the lesson"
  click F "Solution.ipynb" "A python notebook that contains a full set of solutions to the exercises that are discussed in the masterclass.  This notebook is the one that was editted during the section video lecture of the masterclass"
```
