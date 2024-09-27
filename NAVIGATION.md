# Tutorial - Advanced mechanical pulling

This goes beyond what was presented in PLUMED Masterclass 22.15, and which is now available as a Tutorial on [this page](https://www.plumed-tutorials.org/lessons/22/015/data/NAVIGATION.html)

This tutorial provides theoretical background then gives a set of exercises for you to complete. You may want to watch the videos that accompany the previous tutorial first.

The flow chart shown below indicates the order in which you should consult the resources.  You can click on the nodes to access the various resources.  Follow the thick black lines for the best results.  The resources that are connected by dashed lines are supplmentary resources that you may find useful when completing the exercise.

Some relevant perspective articles to consider are [this one](https://doi.org/10.1021/acs.jpcb.1c06330) from Hocky and [this one](https://doi.org/10.1021/acs.jpcb.1c10715) from Stirnemann.

Solutions are available from [This page](https://github.com/hockyg/plumed-tutorial-force2)


```mermaid
flowchart TB;
  A[syntax] -.-> C[Lecture I];
  B[paper1] -.-> C;
  C ==> D[instructions];
  D ==> E[Lecture II];
  E ==> F[paper2];
  click A "ref1" "A previous tutorial that introduces the basics of PLUMED syntax";
  click B "paper1" "A paper describing the FISST method that is used in this tutorial";
  click C "paper2" "A follow up paper combining FISST with replica exchange";
  click D "video1" "A lecture that was given on October 17th 2022 as part of the plumed masterclass series that introduces you to the exercises in this lesson";
  click E "video2" "A lecture that was given on October 24th 2022 as part of the plumed masterclass series that gives background on modelling mechanobiological processes, and then goes through the solutions to the exercises in the previous lesson";
  click F "INSTRUCTIONS.md" "Instructions for the exercises that you are supposed to complete";
```
