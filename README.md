# LP-OPF-CE

The following points should be noted:

  1、For case 2 to case 6, due to the low solution efficiency of the W matrix, therefore, here we provide a method to accelerate the solution: we adopted dual technology, and the dualization of 
  the problem is  detailed in <<dual-OPF-CE.pdf>>. It should be noted that Reference [Zero Duality Gap in Optimal Power Flow Problem] (doi: 10.1109/TPWRS.2011.2160974) has proven 
  that the relaxation problem and the dual OPF problem in the main text **satisfy strong duality**, that is, its dual gap is zero. **Therefore,the two solving methods are essentially 
  equivalent.**  In addition, there are many other methods to accelerate the solution of SDP, such as partitioning the W matrix and then using the sparsity of the matrix to 
  accelerate the solution [Solution of Optimal Power Flow Problems by Semi-definite Programming] (doi: 10.13334/j.0258-8013.pcsee.2008.19.006).

  2、The values in Table 4b do not conform to the rule that "IPM value ≥ global optimal value ≥ SDP value". The reason is as follows: For the OPF problem without deployed control devices, we used the SDP solver built in MATPOWER for solution; however, when using the independently implemented dual-OPF algorithm, the obtained results conform to the above-mentioned rule. Therefore, we speculate that there is a slight deviation in the SDP solver built in MATPOWER for this test case. For the relevant verification process, please refer to the video.

If you have any questions about the manuscript or code, please contact 1743773910@qq.com or leave a message under this repository.
