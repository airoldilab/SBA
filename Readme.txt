================================================
Stochastic Blockmodel Approximation of a Graphon


================================================
This MATLAB package is a supplement to the paper 

E. M. Airoldi, T. B. Costa, and S. H. Chan, "Stochastic blockmodel approximation of a graphon: Theory and consistent estimation", Advances in Neural Information Processing Systems 2013.


================================================
Content:

1. Construct Graphs from a Graphon
	Method 1: [G P u] = construct_a_graph(w,n,T)
		Input:  w - a Graphon 
			n - number of nodes
			T - number of observations
		Output: G - graph (size nxnxT)
			P - probability of each node
			u - label indices

	Method 2: G = construct_a_graph_from_P(P,n,T)
		Input:  P - probability of each node
			n - number of nodes
			T - number of observations
		Output: G - graph (size nxnxT)

2. Stochastic Blockmodel Approximation
	Step 1: B = estimate_blocks_directed(G,Delta)
		Input:  G 	- graph
			Delta 	- the threshold parameter (see Demo_Crossvalidation.m)
		Output: B 	- clusters/blocks
	Step 2: [H,P] = histogram3D(G,B)
		Input:  G 	- graph
			B 	- estimated clusters/blocks
		Output: H	- estimated histogram
			P       - estimated probability of each node (ie graphon)
				(remark, the estimated P is *not* canonical)

3. Cross validation
	Please check Demo_crossvalidation.m


4. Results reported in the paper
	Fig2a.m	Mean Absolute Error vs Number of Nodes
	Fig2b.m	Mean Absolute Error vs Number of Observations
	Fig3a.m Mean Absolute Error vs Number of Blocks
	Fig3b.m Mean Absolute Error vs Percentage of Missing Links
	Fig4    Mean Absolute Error for two types of graphons


5. Compared Methods
	(i)   Largest Gap [1] (estimate_blocks_largest_gap.m)
	(ii)  Universal Singular Value Thresholding [2] (Method_chatterjee.m)
	(iii) Matrix Completion [3] (Method_matrix_completion.m)

	
References
[1] A. Channarond, J. Daudin, and S. Robin. Classification and estimation in the Stochastic Blockmodel based on the empirical degrees. Electronic Journal of Statistics, 6:2574–2601, 2012.
[2] S. Chatterjee. Matrix estimation by universal singular value thresholding. ArXiv:1212.1247. 2012.
[3] R.H. Keshavan, A.Montanari, and S. Oh. Matrix completion from a few entries. IEEE Trans. Information Theory, 56:2980–2998, Jun. 2010.




================================================
COPYRIGHT (C) 2013 Edoardo Airoldi, Thiago Costa, Stanley Chan

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.



================================================
Last update: November 17, 2013



