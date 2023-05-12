There are a few issues with the code in this folder:

1) The analysis matrix is not constructed correctly. This may not end up making a difference in the end because the eigenvalues of the corresponding frame matrix I believe coincide with the eigenvalues of the correct frame matrix

2) The eigenvalue computation is done in an expensive way: all we need to do is form the graph heat kernel and compute the norms of its columns