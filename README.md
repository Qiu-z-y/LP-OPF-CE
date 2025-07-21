# LP-OPF-CE

The following points should be noted:

  1、For large-scale cases (IEEE 59/141 system), due to the low solution efficiency of the W matrix, therefore, here we provide a method to accelerate the solution: we adopted dual technology, and the dualization of 
  the problem is  detailed in <<dual-OPF-CE.pdf>>. It should be noted that Reference [Zero Duality Gap in Optimal Power Flow Problem] (doi: 10.1109/TPWRS.2011.2160974) has proven 
  that the relaxation problem and the dual OPF problem in the main text **satisfy strong duality**, that is, its dual gap is zero. **Therefore,the two solving methods are essentially 
  equivalent.**  In addition, there are many other methods to accelerate the solution of SDP, such as partitioning the W matrix and then using the sparsity of the matrix to 
  accelerate the solution [Solution of Optimal Power Flow Problems by Semi-definite Programming] (doi: 10.13334/j.0258-8013.pcsee.2008.19.006).

  2、For Case 1, we stipulate that a solution is deemed successful when the computational accuracy reaches approximately 7e-4. At this accuracy level, the absolute difference between the results and those obtained with an accuracy of 1e-7 is extremely small in the vast majority of cases (mostly ≤ 0.1).

  3、The empirical running times in the text are all based on the times obtained by the second acceleration method.

If you have any questions about the manuscript or code, please contact 1743773910@qq.com or leave a message under this repository.
